
#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <Eigen/Core>
#include <cstddef>
#include <string_view>
#include <type_traits>
#include <utility>

#include "libfly/utility/core.hpp"

/**
 * \file atom.hpp
 *
 * @brief The basic building block for atoms in libFLY.
 */

namespace fly::system {

  namespace detail {

    template <typename Tag>
    struct AtomMem {
    public:
      static_assert(std::is_empty_v<Tag>, "Property tag types are required to be empty");

      AtomMem() = default;

      AtomMem(AtomMem&&) noexcept = default;

      AtomMem(AtomMem const&) = default;

      explicit AtomMem(typename Tag::matrix_t&& data) : m_data(std::move(data)) {}

      explicit AtomMem(typename Tag::matrix_t const& data) : m_data(data) {}

      AtomMem& operator=(AtomMem const&) = default;

      AtomMem& operator=(AtomMem&&) = default;

      template <typename T>
      explicit AtomMem(T&& x) : m_data(std::forward<T>(x)) {}

      typename Tag::matrix_t const& operator[](Tag) const { return m_data; }

      typename Tag::matrix_t& operator[](Tag) { return m_data; }

    private:
      typename Tag::matrix_t m_data;
    };

  }  // namespace detail

  /**
   * @brief The libfly representation of an atom.
   *
   * Libfly uses this type to build atoms to allow integration with libFLY's containers.
   * An Atom behaves as a structure of its properties , which are accessed through ``operator[]``.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/system/atom.cpp
   *    :code:
   *
   * \endrst
   *
   * @tparam T a series of empty types, derived from ``Property``, to describe each member.
   */
  template <typename... T>
  struct Atom : detail::AtomMem<T>... {
    /**
     * @brief Default construct a new Atom object.
     */
    Atom() = default;

    /**
     * @brief Copy construct a new Atom object.
     */
    Atom(Atom&&) noexcept = default;

    /**
     * @brief Move construct a new Atom object.
     */
    Atom(Atom const&) = default;

    /**
     * @brief Default assignment operator.
     */
    Atom& operator=(Atom const&) = default;

    /**
     * @brief Default move-assignment operator.
     */
    Atom& operator=(Atom&&) = default;

    /**
     * @brief Construct a new Atom object, initializing each member.
     */
    explicit Atom(typename T::matrix_t const&... args) : detail::AtomMem<T>(args)... {}

    /**
     * @brief Construct a new Atom object, initializing each member.
     */
    explicit Atom(typename T::matrix_t&&... args) : detail::AtomMem<T>(args)... {}

    using detail::AtomMem<T>::operator[]...;
  };

}  // namespace fly::system
