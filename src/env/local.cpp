// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/local.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

#include "graph.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "xxhash.h"

// static_assert(!is_narrowing_conversion_v<Colour::scalar_t, int>, "Colour's scalar_t is too big for Nauty");

namespace fly::env {

  bool Local::Fingerprint::equiv(Fingerprint const &other, double delta) const {
    if (m_r_0j.size() != other.m_r_0j.size()) {
      return false;
    }

    for (std::size_t i = 0; i < m_r_0j.size(); i++) {
      if (std::abs(m_r_0j[i] - other.m_r_0j[i]) > delta * M_SQRT2) {
        return false;
      }
    }

    ASSERT(m_r_ij.size() == other.m_r_ij.size(), "Secondary distances should match {}!={}", m_r_ij.size(), other.m_r_ij.size());

    for (std::size_t i = 0; i < m_r_ij.size(); i++) {
      if (std::abs(m_r_ij[i] - other.m_r_ij[i]) > delta * M_SQRT2) {
        return false;
      }
    }

    return true;
  }

  // Map id and frozen to a colour.
  int to_colour_idx(TypeID::scalar_t id, Frozen::scalar_t frz) { return safe_cast<int>(2 * id + static_cast<TypeID::scalar_t>(frz)); }

  void Local::rebuild(int ix,
                      system::SoA<TypeID const &, Frozen const &> atoms,
                      neigh::List const &nl,
                      Eigen::Index num_types,
                      double r_env,
                      double r_edge) {
    // // Reuses our storage
    // this->clear();
    // m_fingerprint.m_r_0j.clear();
    // m_fingerprint.m_r_ij.clear();

    // //
    // std::vector<int> offsets(safe_cast<std::size_t>(2 * num_types), 0);

    // // Add central atom
    // this->emplace_back({0, 0, 0}, to_colour_idx(atoms(id_, ix), atoms(fzn_, ix)), ix);

    // // Build the local environment and m_r_0j
    // nl.for_neighbours(ix, r_env, [&](auto n, double r, Vec const &dr) {
    //   //
    //   int col = to_colour_idx(atoms(id_, n), atoms(fzn_, n));

    //   this->emplace_back(dr, col, n);

    //   m_fingerprint.m_r_0j.push_back(r);

    //   ASSERT(col < 2 * num_types, "Colour out of bounds, {}", atoms(id_, n));

    //   offsets[safe_cast<std::size_t>(col)]++;
    // });

    // std::sort(m_fingerprint.m_r_0j.begin(), m_fingerprint.m_r_0j.end());  // Done

    // Graph initial(safe_cast<int>(size()));

    // // Build r_ij part of the fingerprint and fill out the graph

    // for (int i = 0; i < size(); i++) {
    //   for (int j = 0; j < i; j++) {
    //     //
    //     double r = gnorm((*this)[i][r_] - (*this)[j][r_]);

    //     if (r < r_edge) {
    //       initial.make_edge(i, j);
    //     }
    //     m_fingerprint.m_r_ij.push_back(r);
    //   }
    // }

    // std::sort(m_fingerprint.m_r_ij.begin(), m_fingerprint.m_r_ij.end());  // Done

    // // Build the key from the graph.

    // Graph::label_t labels;

    // Graph canon = initial.canon(labels);

    // static_assert(!is_narrowing_conversion_v<XXH64_hash_t, std::size_t>, "Hash is too big");

    // XXH64_hash_t hash = XXH64(colour_counts.data(), colour_counts.size() * sizeof(int), 0);

    // m_key = hash ^ canon.hash();

    // // Set COM == 0,0,0
    // Vec shift = centroid(*this);

    // for (auto &elem : *this) {
    //   elem[r_] -= shift;
    // }
  }

}  // namespace fly::env