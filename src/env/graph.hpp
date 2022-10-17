#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "libfly/env/geometry.hpp"
#include "libfly/env/local.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/utility/core.hpp"

//

#include "xxhash.h"

#define MAXN 128

#include "nauty/nauty.h"

namespace fly::env {

  class AdjMat {
  public:
    /**
     * @brief Construct a new adjacency matrix object with ``n`` nodes and no edges.
     *
     * @param n The number of nodes in the graph.
     */
    explicit AdjMat(int n) : m_n(n), m_m(SETWORDSNEEDED(n)) {
      ASSERT(n > 0 && n <= MAXN, "Graph with {} nodes is too big, limit={}", n, MAXN);
      clear();
#ifndef NDEBUG
      nauty_check(WORDSIZE, m_m, m_n, NAUTYVERSIONID);
#endif
    }
    /**
     * @brief Get the number of nodes in the graph.
     */
    int n() const noexcept { return m_n; }

    /**
     * @brief Get the number of rows in used to store the adjacency matrix.
     *
     * Required by Nauty algorithms.
     */
    int m() const noexcept { return m_m; };

    /**
     * @brief Removes all edges from the graph.
     */
    void clear() noexcept { EMPTYGRAPH(m_graph.data(), size_t(m_m), m_n); }

    /**
     * @brief Add an edge between nodes ``v`` and ``w``.
     */
    void make_edge(int v, int w) {
      //
      ASSERT(v >= 0 && v < m_n && v >= 0 && w < m_n, "[{},{}] is an invalid edge with {} nodes", v, w, m_n);
      ADDONEEDGE(m_graph.data(), v, w, size_t(m_m));
    }

    /**
     * @brief
     *
     * @return graph*
     */
    graph* data() noexcept { return m_graph.data(); }

    /**
     * @brief Compute a hash of the graphs adjacency matrix.
     */
    XXH64_hash_t hash() const noexcept { return XXH64(m_graph.data(), size_t(m_m) * size_t(m_n) * sizeof(graph), 0); }

  private:
    int m_n;
    int m_m;

    std::array<graph, MAXN * MAXM> m_graph;
  };

  /**
   * @brief Compute the colour offsets.
   *
   *            data:
   *         indexes: [0, 1, 2, 3, 4, 5]
   *    Atom colours: [2, 0, 2, 0, 2, 2] (total four colours: c_max = 4)
   *
   *          counts: [2, 0, 3, 0, 2]    (we skip the first atom so only 3 of colour #2)
   *         offsets: [1, 3, 3, 6, 2]
   *
   * The offsets array is invariant under the order of geo, encodes the colour of the centre atom and the counts of all the other
   * colours in the geometry therefore, its hash encodes all this information.
   *
   */
  template <typename... T>
  std::vector<std::size_t> offsets(Geometry<T...> const& geo, std::size_t c_max) {
    // Store the counts of each colour in 0,..., c_max -1. Slot c_max stores the central colour for hashing.
    std::vector<std::size_t> offsets(c_max + 1, 0);

    offsets.back() = safe_cast<std::size_t>(geo[0][col_]);

    // Skip first, count the colours.
    for (int i = 1; i < geo.size(); i++) {
      ASSERT(::fly::detail::cmp_less(geo[i][col_], c_max), "Colour #{} is too big!", geo[i][col_]);
      offsets[safe_cast<std::size_t>(geo[i][col_])]++;
    }

    std::size_t sum = 1;

    for (std::size_t i = 0; i < c_max; i++) {
      std::size_t tmp = offsets[i];
      offsets[i] = sum;
      sum += tmp;
    }

    return offsets;
  }

  /**
   * @brief Re-order a geometry into its graph-canonical order and produce a hash of its colour + graph-topology.
   */
  template <typename... T>
  XXH64_hash_t canon_hash(Geometry<T...>& geo, double r_edge, std::size_t c_max) {
    //
    ASSERT(geo.size() <= MAXN, "Graph with {} nodes is too big, limit={}", geo.size(), MAXN);

    // Build the adjacency matrix.

    AdjMat mat(safe_cast<int>(geo.size()));

    for (int i = 0; i < geo.size(); i++) {
      for (int j = 0; j < i; j++) {
        if (gnorm(geo[i][r_] - geo[j][r_]) < r_edge) {
          mat.make_edge(i, j);
        }
      }
    }

    // Build the colour offsets array and compute a colour hash.

    std::vector<std::size_t> off = offsets(geo, c_max);

    XXH64_hash_t c_hash = XXH64(off.data(), off.size() * sizeof(std::size_t), 0);

    //  Compute lab and ptn arrays

    // Want the indexes in colour order (with first atom in a special colour class)

    //            data:
    //         indexes: [0, 1, 2, 3, 4, 5, 6]
    //    Atom colours: [2, 0, 2, 0, 2, 2, 3] (total four colours: c_max = 4)
    //
    //          counts: [2, 0, 3, 1, 2]    (we skip the first atom so only 3 of colour #2)
    //         offsets: [1, 3, 3, 6, 2]

    //            want:
    //             lab: [0, 1, 3, 2, 4, 5, 6]
    //             ptn: [0, 1, 0, 1, 1, 0, 0] Zeros @ ends of colour class: 0, 2, 5, 6

    std::array<int, MAXN> ptn;                    // Partitions/colours
    std::fill_n(ptn.begin(), geo.size() - 1, 1);  // Initialize with contiguous colour class.

    for (std::size_t i = 0; i < c_max; i++) {  // Mark ends of the partitions.
      ptn[off[i] - 1] = 0;
    }

    ptn[safe_cast<std::size_t>(geo.size() - 1)] = 0;  // Last atom must terminate colour class.

    std::array<int, MAXN> lab;

    lab[0] = 0;
    // Counting sort, skip first as it is in position.
    for (int i = 1; i < geo.size(); i++) {
      auto col = safe_cast<std::size_t>(geo[i][col_]);
      lab[off[col]++] = i;
    }

    // Perform graph canonisation.

    // We have to expand a Nauty macro here as we want a constexpr variable.
    static constexpr optionblk options = {
        1, 0, 0, 0, 0, 0, 78, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 100, 0, 1, 0, &dispatch_graph, 0, NULL,
    };

    std::array<graph, MAXN * MAXM> canon;
    statsblk stats;
    std::array<int, MAXN> orbits;

    densenauty(mat.data(),
               lab.data(),
               ptn.data(),
               orbits.data(),
               const_cast<optionblk*>(&options),  // c api does not support const.
               &stats,
               mat.m(),
               mat.n(),
               canon.data());

    // fmt::print("Group size = {} x 10^{}\n", stats.grpsize1, stats.grpsize2);
    // fmt::print("Number of orbits = {}\n", stats.numorbits);
    // fmt::print("Number of generators = {}\n", stats.numgenerators);

    // Now need to re-order the input geometry.

    //          indexes: [0, 1, 2, 3, 4, 5, 6]
    // lab (post canon): [0, 6, 4, 3, 2, 1, 5] means 0->0  6->1 4->2
    //

    //          indexes:  [0, 6, 4, 3, 2, 1, 5]
    //                    [0, 5, 2, 3, 4, 6, 1]

    //              gen:  [0, 1, 2, 3, 4, 5, 6]
    //

    ASSERT(lab[0] == 0, "Colouring has failed lab[0]={}", lab[0]);

    Geometry buff(geo);

    for (int i = 0; i < geo.size(); i++) {
      geo[i] = buff[lab[std::size_t(i)]];
    }

    return XXH64(canon.data(), size_t(mat.m()) * size_t(mat.n()) * sizeof(graph), c_hash);
  }

}  // namespace fly::env