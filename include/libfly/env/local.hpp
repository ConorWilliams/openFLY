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

#include "libfly/env/geometry.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file local.hpp
 *
 * @brief Local environments.
 *
 * The local environment (LE) of an atom is the set of atoms within some neighbourhood. It is assumed in OLKMC that the mechanisms
 * accessible to an atom are completely contained-within and solely a-function-of its LE.
 */

// namespace fly::env {

//   /**
//    * @brief A local environment is a (localised) geometry augmented with a key and a fingerprint.
//    *
//    * Here localised means the first atom is the central (i.e. fixed) atom.
//    *
//    * The key is discrete representation of this topology.
//    * The fingerprint is a continuous representation that is cheap(ish)ly comparable.
//    */
//   class Local : public Geometry<Index> {
//   public:
//     /**
//      * @brief The discrete representation of this environment.
//      *
//      * The first pair is the colour of the centre atom (for .first and .second = 1), the remaining
//      * pairs are .first = colour and .second = count.
//      */
//     using Key = std::vector<std::pair<std::size_t, std::size_t>>;

//     /**
//      * @brief An ordered representation of the intra-atomic distances in this Geometry.
//      */
//     class Fingerprint {
//     public:
//       /**
//        * @brief A fast test to see if two local environments **may** be equivalent.
//        *
//        * Explicitly this tests that the paired ordered intra-atomic distances are within "delta".
//        * If using this as a filter for permute_** algorithms then use "delta" = M_SQRT2 *
//        * delta_for_perm.
//        */
//       friend bool equiv(Fingerprint const& a, Fingerprint const& b, double delta);

//       floating r_min() const {
//         ASSERT(!m_r_0j.empty() && !m_r_ij.empty(), "Not enough atoms for f_min!");
//         return std::min(m_r_0j[0], m_r_ij[0]);
//       }

//     private:
//       friend class Local;

//       std::vector<floating> m_r_0j;
//       std::vector<floating> m_r_ij;
//     };

//     /**
//      * @brief Get a const reference to the discrete key.
//      */
//     Key const& key() const noexcept { return m_key; }

//     /**
//      * @brief Get a const reference to the fingerprint.
//      */
//     Fingerprint const& fingerprint() const noexcept { return m_fingerprint; }

//     /**
//      * @brief Rebuild this local environment to be the env of the idx'th atom.
//      */
//     void rebuild(SimCell const& atoms, neighbour::List const& nl, std::size_t idx);

//   private:
//     Key m_key;

//     Fingerprint m_fingerprint;
//   };

//   /**
//    * @brief This class builds and stores a list of local environments from a SimCell.
//    */
//   class EnvCell {
//   public:
//     struct Options {
//       /** @brief Radius of environment. */
//       floating r_env;
//     };

//     /**
//      * @brief Construct a new Env Cell object.
//      */
//     EnvCell(Options const& opt, OrthoSimBox const& box) : m_nl(box, opt.r_env){};

//     /**
//      * @brief Get the number of Environments in the cell.
//      */
//     std::size_t size() const noexcept { return m_envs.size(); }

//     /**
//      * @brief Get the ith local environment
//      */
//     Local const& operator[](std::size_t i) const noexcept { return m_envs[i]; }

//     /**
//      * @brief Get the ith local environment
//      */
//     Local& operator[](std::size_t i) noexcept { return m_envs[i]; }

//     /**
//      * @brief Construct the local environemnts around all the atom.
//      */
//     void rebuild(SimCell const& atoms, std::size_t num_threads);

//   private:
//     std::vector<Local> m_envs;

//     neighbour::List m_nl;
//   };
// }  // namespace fly::env