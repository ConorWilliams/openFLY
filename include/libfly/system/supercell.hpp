#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <type_traits>

#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
/**
 * \file supercell.hpp
 */

namespace fly::system {

  /**
   * @brief LibFLY's representation of a system of atoms
   *
   * The Supercell is libFLY's amalgamation of all the data required for a simulation. It **is a** SoA containing all the atoms in the
   * system **has a** ``Box`` and a ``TypeMap``. A Supercell always has the ``TypeID`` property for each atom and accepts the rest as
   * template arguments.
   *
   * @tparam Map A specialisation of fly::system::TypeMap.
   * @tparam T Property tags derived from ``Property``.
   */
  template <typename Map, typename... T>
  class Supercell : public SoA<TypeID, T...> {
  private:
    using SOA = SoA<TypeID, T...>;

    static_assert(SOA::owns_all, "Supercells must own all their data");

    static_assert(detail::is_TypeMap<Map>::value, "The Map template param must be a specialisation of system::TypeMap");

  public:
    /**
     * @brief Construct a new Supercell object to store ``num_atoms`` atoms.
     */
    Supercell(Box const& box, Map const& map, Eigen::Index num_atoms) : SOA(num_atoms), m_box(box), m_map(map) {}

    /**
     * @brief Simulation box const-getter.
     */
    Box const& box() const { return m_box; }

    /**
     * @brief Simulation box getter.
     */
    Box& box() { return m_box; }

    /**
     * @brief TypeMap const-getter.
     */
    Map const& map() const { return m_map; }

    /**
     * @brief TypeMap getter.
     */
    Map& map() { return m_map; }

  private:
    Box m_box;
    Map m_map;
  };

  /**
   * @brief Utility to construct a new Supercell to store `num_atoms` atoms.
   *
   * Uses partial function-template argument deduction to deduce Supercell's template parameters.
   *
   * @tparam Property tags derived from ``Property``.
   *
   * @param box Input box to forward to Supercell constructor.
   * @param map Input map to forward to Supercell constructor.
   * @param num_atoms Number of atoms to forward to Supercell constructor.
   *
   * @return A constructed supercell with the Map template-parameter deduced.
   */
  template <typename... T, typename... U>
  auto make_supercell(Box const& box, TypeMap<U...> const& map, Eigen::Index num_atoms) -> Supercell<TypeMap<U...>, T...> {
    return {box, map, num_atoms};
  }

}  // namespace fly::system