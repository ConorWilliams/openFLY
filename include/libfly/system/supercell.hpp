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
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/utility/core.hpp"
/**
 * \file supercell.hpp
 *
 * @brief Classes to represent systems of atoms.
 */

namespace fly::system {

  template <typename...>
  class TypeMap;

  namespace detail {

    template <typename...>
    struct same_properties : std::false_type {};

    template <typename... T>
    struct same_properties<TypeMap<T...>, T...> : std::true_type {};

  }  // namespace detail

  /**
   * @brief Utility for mapping an atoms TypeID to atomic properties.
   *
   * This is used for properties that are the same for many atoms e.g. atomic numbers.
   *
   * @tparam T Tags derived from ``Property``, to describe each property.
   */
  template <typename... T>
  class TypeMap : private SoA<Type, T...> {
  private:
    using SOA = SoA<Type, T...>;

    static_assert(SOA::owns_all, "TypeMap must own all its data,");

  public:
    /**
     * @brief Defaulted copy constructor.
     */
    TypeMap(TypeMap const&) = default;

    /**
     * @brief Defaulted move constructor.
     */
    TypeMap(TypeMap&&) = default;

    /**
     * @brief Construct a TypeMap to hold ``num_types`` types.
     *
     * \rst
     * .. warning::
     *    All the new properties are uninitialised.
     * \endrst
     */
    explicit TypeMap(Eigen::Index num_types) : SOA(num_types) {}

    /**
     * @brief Construct a new TypeMap by slicing a different kind of TypeMap.
     */
    template <typename... U, typename = std::enable_if_t<!detail::same_properties<TypeMap, U...>::value>>
    explicit TypeMap(TypeMap<U...> map) : SOA(static_cast<SoA<Type, U...>>(std::move(map))) {}

    /**
     * @brief Fetch the number of types stored in the TypeMap
     */
    Eigen::Index num_types() const { return SOA::size(); }

    /**
     * @brief Get a property corresponding to the id ``id``.
     */
    template <typename U>
    decltype(auto) get(U, std::uint32_t id) const {
      return SOA::operator()(U{}, safe_cast<Eigen::Index>(id));
    }

    /**
     * @brief Set a property corresponding to the id ``id``.
     */
    template <typename U, typename V = typename remove_cref_t<U>::matrix_t>
    void set(U, std::uint32_t id, V&& value) {
      SOA::operator()(U{}, safe_cast<Eigen::Index>(id)) = std::forward<V>(value);
    }

  private:
    /* clang-format off */ template <typename...>  friend class TypeMap; /* clang-format on */
  };

  /**
   * @brief LibFLY's representation of a system of atoms
   *
   * The Supercell is libFLY's amalgamation of all the data required for a simulation. It **is a** SoA containing all the atoms in the
   * system **has a** ``Box`` and a ``TypeMap``. A Supercell always has the ``TypeID`` property for each atom and accepts the rest as
   * template arguments.
   *
   * @tparam T Tags derived from ``Property``, to describe each member.
   */
  template <typename Map, typename... T>
  class Supercell : public SoA<TypeID, T...> {
  private:
    using SOA = SoA<TypeID, T...>;

    static_assert(SOA::owns_all, "Supercells must own all their data");

  public:
    /**
     * @brief Construct a new Supercell object to store `num_atoms` atoms.
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
   */
  template <typename... T, typename... U>
  Supercell<TypeMap<U...>, T...> make_supercell(Box const& box, TypeMap<U...> const& map, Eigen::Index num_atoms) {
    return {box, map, num_atoms};
  }

}  // namespace fly::system