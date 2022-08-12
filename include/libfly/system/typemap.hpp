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
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
/**
 * \file typemap.hpp
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
   * @brief Utility for mapping an atom's TypeID to a Type and its properties.
   *
   * This is used for properties that are the same for many atoms e.g. atomic numbers, by default each TypeID must always have a Type
   * property.
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
    TypeMap(TypeMap&&) noexcept = default;

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

}  // namespace fly::system