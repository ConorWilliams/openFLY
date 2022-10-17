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

  bool Fingerprint::equiv(Fingerprint const &other, double delta) const {
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
    // Reuses our storage
    this->clear();
    m_fingerprint.m_r_0j.clear();
    m_fingerprint.m_r_ij.clear();

    // Add central atom
    this->emplace_back({0, 0, 0}, to_colour_idx(atoms(id_, ix), atoms(fzn_, ix)), ix);

    // Build the local environment and m_r_0j
    nl.for_neighbours(ix, r_env, [&](auto n, double r, Vec const &dr) {
      this->emplace_back(dr, to_colour_idx(atoms(id_, n), atoms(fzn_, n)), n);
      m_fingerprint.m_r_0j.push_back(r);
    });
    std::sort(m_fingerprint.m_r_0j.begin(), m_fingerprint.m_r_0j.end());  // Done

    // Build r_ij part of the fingerprint i != 0, j < i

    for (int i = 1; i < size(); i++) {
      for (int j = 1; j < i; j++) {
        m_fingerprint.m_r_ij.push_back(gnorm((*this)[i][r_] - (*this)[j][r_]));
      }
    }
    std::sort(m_fingerprint.m_r_ij.begin(), m_fingerprint.m_r_ij.end());  // Done

    // // Build the key from the graph.

    m_key = canon_hash(*this, r_edge, 2 * size_t(num_types));

    // Set COM == 0,0,0
    Vec shift = centroid(*this);

    for (auto &elem : *this) {
      elem[r_] -= shift;
    }
  }

  void LocalList::rebuild(system::SoA<Position const &, TypeID const &, Frozen const &> const &info, int num_threads) {
    //
    m_nl.rebuild(info, num_threads);

    m_envs.resize(size_t(info.size()));

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int i = 0; i < info.size(); i++) {
      m_envs[size_t(i)].rebuild(i, info, m_nl, m_num_types, m_opt.r_env, m_opt.r_edge);
    }
  }

}  // namespace fly::env