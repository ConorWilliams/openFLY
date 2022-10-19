#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.centroid>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <functional>
#include <memory>

#include "libfly/env/geometry.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file heuristics.hpp
 *
 * @brief Utilities for quickly comparing local environments without resorting to a full ``Geometry::permute_onto()``.
 */

namespace fly::env {

  /**
   * @brief Reorder a geometry into a canonical order and produce a hash encoding the colours and topology of the geometry.
   *
   * This functions builds a graph based representation of the geometry, encoding atoms as coloured nodes joined with edges if the
   * atoms are closer than ``r_edge``. This graph is then canonised using the nauty library (see: https://pallini.di.uniroma1.it/) and
   * a hash is derived from the canonical adjacency matrix.
   *
   * @param geo The geometry to be reordered.
   * @param r_edge The edge length used to determine if atoms in the geometry are bonded in the graph representation.
   * @param c_max The number of colours in the simulation, equal to ``2 * TypeMap::num_types()``.
   * @param scratch An optional parameter, if provided and ``scratch->size() >= geo.size()`` then no allocation will occur.
   */
  std::size_t canon_hash(Geometry<Index>& geo, double r_edge, std::size_t c_max, Geometry<Index>* scratch = nullptr);

  /**
   * @brief An ordered representation of the intra-atomic distances in a Geometry.
   *
   * The fingerprint is a continuous representation that is cheap(ish)ly comparable.
   */
  class Fingerprint {
  public:
    /**
     * @brief A fast test to see if two local environments **may** be equivalent.
     *
     * Explicitly, for the LE this fingerprint represents, this function must return ``true`` if the LE is able to be
     * ``permute_onto()`` other with the same ``delta``.
     *
     *
     * @param other The other fingerprint.
     * @param delta The tolerance for the equivalence.
     * @return ``true`` If two local environments **may** be equivalent.
     * @return ``false``  If two local environments are not equivalent.
     */
    auto equiv(Fingerprint const& other, double delta) const -> bool;

    /**
     * @brief Get the smallest intra-atomic separation in this environment.
     */
    auto r_min() const -> double {
      ASSERT(!m_r_0j.empty() && !m_r_ij.empty(), "Not enough atoms for r_min!", 0);
      return std::min(m_r_0j.front(), m_r_ij.front());
    }

    /**
     * @brief Rebuild the fingerprint from a Geometry.
     */
    template <typename... T>
    auto rebuild(Geometry<T...> const& geo) -> void {
      // Reuse storage.
      m_r_0j.clear();
      m_r_ij.clear();

      for (int i = 1; i < geo.size(); i++) {
        // Build r_0i part of the fingerprint i != 0
        m_r_0j.push_back(gnorm(geo[i][r_] - geo[0][r_]));

        // Build r_ij part of the fingerprint i != 0, j < i
        for (int j = 1; j < i; j++) {
          m_r_ij.push_back(gnorm(geo[i][r_] - geo[j][r_]));
        }
      }

      std::sort(m_r_0j.begin(), m_r_0j.end());  // Done
      std::sort(m_r_ij.begin(), m_r_ij.end());  // Done
    }

  private:
    std::vector<double> m_r_0j;
    std::vector<double> m_r_ij;
  };

}  // namespace fly::env
