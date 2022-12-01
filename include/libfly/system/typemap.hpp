#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <string_view>
#include <type_traits>

#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file typemap.hpp
 */

namespace fly::io {
  class BinaryFile;  // Forward declaration for friendship.
}

namespace fly::system {

  template <typename...>
  class TypeMap;

  namespace detail {

    template <typename...>
    struct same_properties : std::false_type {};

    template <typename... T>
    struct same_properties<TypeMap<T...>, T...> : std::true_type {};

    template <typename>
    struct is_TypeMap : std::false_type {};

    template <typename... T>
    struct is_TypeMap<TypeMap<T...>> : std::true_type {};

  }  // namespace detail

  /**
   * @brief Utility for mapping an atom's TypeID to a Type and its properties.
   *
   * This is used for properties that are the same for many atoms e.g. atomic numbers, by default each TypeID must always have a Type
   * property.
   *
   * \rst
   * .. note::
   *    TypeMap has special handling for getting and setting ``Type``'s ``Type::matrix_t`` to string-like types for QoL.
   * \endrst
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
     * @brief Defaulted copy-assignment.
     */
    TypeMap& operator=(TypeMap const&) = default;

    /**
     * @brief Defaulted move-assignment.
     */
    TypeMap& operator=(TypeMap&&) = default;

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
    template <typename... U,
              typename = std::enable_if_t<!detail::same_properties<TypeMap, U...>::value
                                          && std::is_constructible_v<SOA, SoA<Type, U...> const&>>>
    explicit TypeMap(TypeMap<U...> const& map) : SOA(static_cast<SoA<Type, U...> const&>(map)) {}

    /**
     * @brief Fetch the number of types stored in the TypeMap
     */
    auto num_types() const -> Eigen::Index { return SOA::size(); }

    /**
     * @brief Get a property form the map.
     *
     * @param id The TypeID of the type who's property you are getting.
     * @param tag Tag to disambiguate property you would like.
     */
    template <typename U>
    auto get(std::uint32_t id, [[maybe_unused]] U tag) const -> decltype(auto) {
      return SOA::operator()(U{}, safe_cast<Eigen::Index>(id));
    }

    /**
     * @brief Overload of get() for ``Type`` that returns a string_view of the type name.
     *
     * \rst
     *
     * .. tip::
     *    With ``truncate = true`` this overload has some runtime-cost to count the char's in the array.
     *
     *    If you would like to call the template overload to get a mapped Eigen::Array use explicit template syntax:
     *
     *    .. code:: cpp
     *
     *       TypeMap map(1);
     *
     *       // set the Type //
     *
     *       map.get<Type>(0, tp_)
     *
     * \endrst
     *
     * @param truncate If ``true`` truncate the returned string_view to the non-null portion.
     * @param id The typeID of the type who's property you are getting.
     */
    auto get(std::uint32_t id, Type, bool truncate = true) const -> std::string_view {
      if (truncate) {
        return {SOA::operator()(Type{}, safe_cast<Eigen::Index>(id)).data()};
      } else {
        return {SOA::operator()(Type{}, safe_cast<Eigen::Index>(id)).data(), Type::size()};
      }
    }

    /**
     * @brief Set a property in the map.
     *
     * @param id The TypeID of the type who's property you are getting.
     * @param tag Tag to disambiguate property you would like.
     * @param value Value to write into the map.
     */
    template <typename U, typename V = typename U::matrix_t>
    auto set(std::uint32_t id, [[maybe_unused]] U tag, V&& value)
        -> std::enable_if_t<std::is_assignable_v<typename U::matrix_t&, V&&>> {
      SOA::operator()(U{}, safe_cast<Eigen::Index>(id)) = std::forward<V>(value);
    }

    /**
     * @brief Overload of set() for Type that accepts string_views.
     *
     * @param id The TypeID of the type who's property you are getting.
     * @param tag Tag to disambiguate property you would like.
     * @param value Value to write into the map.
     */
    auto set(std::uint32_t id, [[maybe_unused]] Type tag, std::string_view value) -> void {
      //
      if (!::fly::detail::cmp_less(value.size(), Type::size())) {
        throw error("Type name length {} but must be < {}", value.size(), Type::size());
      }

      char* buf = (*this)(Type{}, safe_cast<Eigen::Index>(id)).data();

      for (auto&& elem : value) {
        *(buf++) = elem;
      }

      *buf = '\0';
    }

    /**
     * @brief  Set all the properties of a typeID
     *
     * @param id The TypeID of the type who's property you are getting.
     * @param type The value of the Type to write into the map.
     * @param args The rest of values to write into the map.
     */
    auto set(std::uint32_t id, std::string_view type, typename T::matrix_t const&... args) -> void {
      (set(id, Type{}, type));
      ((set<T>(id, T{}, args)), ...);
    }

  private:
    // Friend for converting constructor.
    template <typename...>
    friend class TypeMap;

    // Friend for IO.
    friend class ::fly::io::BinaryFile;
  };

}  // namespace fly::system