// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/local.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/utility/core.hpp"
#include "xxhash.h"

#define MAXN 128

#include "nauty/nauty.h"

namespace fly::env {

  /**
   * @brief Encapsulate call to the C-lang Nauty library.
   *
   */
  class Graph {
  public:
    /**
     * @brief Construct a new Graph object with ``n`` nodes and no edges.
     *
     * @param n The number of nodes in the graph.
     */
    Graph(int n) : m_n(n), m_m(SETWORDSNEEDED(n)) {
      ASSERT(n > 0 && n <= MAXN, "Graph with {} nodes is too big, limit={}", n, MAXN);
      clear();
#ifndef NDEBUG
      nauty_check(WORDSIZE, m_m, m_n, NAUTYVERSIONID);
#endif
    }
    /**
     * @brief Get the number of nodes in the graph.
     */
    int size() const noexcept { return m_n; }

    /**
     * @brief Removes all edges from the graph.
     */
    void clear() noexcept { EMPTYGRAPH(m_graph.data(), m_m, m_n); }

    /**
     * @brief Add an edge between nodes ``v`` and ``w``.
     */
    void make_edge(int v, int w) {
      //
      ASSERT(v >= 0 && v < m_n && v >= 0 && w < m_n, "[{},{}] is an invalid edge", v, w);

      ADDONEEDGE(m_graph.data(), v, w, m_m);
    }

    /**
     * @brief Compute a hash of the graphs adjacency matrix.
     */
    XXH64_hash_t hash() const noexcept { return XXH64(m_graph.data(), size_t(m_m) * size_t(m_n) * sizeof(graph), 0); }

    /**
     * @brief A type used to communicate the ordering of a canonical permutation.
     */
    using label_t = std::array<int, MAXN>;

    /**
     * @brief Produce a new graph, isomorphic to this graph, in a canonical form.
     *
     * @param label Lists the vertices of this graph in the order in which they need to be relabelled to give the canonical graph, only
     * the first ``size()`` elements are initialized.
     */
    Graph canon(label_t& label) const {
      // We have to expand a Nauty macro here as this is a constexpr variable.
      static constexpr optionblk options = {
          1, 0, 0, 0, 1, 0, 78, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 100, 0, 1, 0, &dispatch_graph, 0, NULL,
      };

      Graph out(size());
      statsblk stats;
      int orbits[MAXN];
      int ptn[MAXN];  // ignored as options specifies default colouring

      // Const casting as C api does not support const.
      densenauty(const_cast<Graph*>(this)->m_graph.data(),
                 label.data(),
                 ptn,
                 orbits,
                 const_cast<optionblk*>(&options),
                 &stats,
                 m_m,
                 m_n,
                 out.m_graph.data());

      //   fmt::print("Group size = {} x 10^{}\n", stats.grpsize1, stats.grpsize2);
      //   fmt::print("Number of orbits = {}\n", stats.numorbits);
      //   fmt::print("Number of generators = {}\n", stats.numgenerators);

      return out;
    }

  private:
    int m_n;
    int m_m;

    std::array<graph, MAXN * MAXM> m_graph;
  };

  // Need a function that accepts a VoS, r_edge and a lambda that computes a graph hash and re-orders the VoS into canonical order.
  // Assumes the colours are 0, 1 ,... , colour_max
  template <typename F, typename... T>
  std::size_t canon_hash(Geometry<T...>& geo, double r_edge, F const& callback) {
    //
    //          indexes: [0, 1, 2, 3, 4, 5]
    //     Atom colours: [2, 0, 2, 0, 2, 2] (total four colours)
    //          offsets: [2, 0, 3, 0]       (we skip the first atom so only 3 of colour #2)
    //      sum offsets: [0, 2, 2, 5]

    // Want the indexes in colour order (with first atom in a special colour class)

    //             lab: [0, 1, 3, 2, 4, 5]
    //             ptn: [0, 1, 0, 1, 1, 0] Zeros @ ends of colour class: 0, 2, 5

    // 1 + sum offsets: [1, 3, 3, 6]

    std::array<int, MAXN> ptn;  // Partitions/colours
    ptn.fill(1);
    ptn[0] = 0;  // First atom is a special colour class.

    int sum = 0;

    // for (int& elem : offsets) {
    //   int tmp = elem;
    //   elem = sum;
    //   sum += tmp;
    // }
  }

}  // namespace fly::env