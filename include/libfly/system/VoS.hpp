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
#include <type_traits>
#include <utility>
#include <vector>

#include "libfly/system/atom.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file VoS.hpp
 *
 * @brief Array of Structures implementation.
 */

namespace fly::system {

  /**
   * @brief A container that models a std::vector of atoms.
   *
   * The members of the "atom" are described through a series of template parameters which should
   * inherit from ``MemTag``
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/system/VoS.cpp
   *    :code:
   *
   * .. todo::
   *    Work out how to make inherited functions display with ``using``.
   *
   * \endrst
   *
   * @tparam Mems a series of empty types, derived from ``MemTag``, to describe each member.
   */
  template <typename... Mems>
  class VoS : private std::vector<Atom<Mems...>> {
  private:
    static_assert(sizeof...(Mems) > 0, "Need at least one member in an VoS");

    using Vector = std::vector<Atom<Mems...>>;

  public:
    // Expose subset of underlying vector API
    using Vector::begin;
    using Vector::clear;
    using Vector::emplace_back;
    using Vector::end;
    using Vector::push_back;
    using Vector::size;
    using Vector::Vector;

    /**
     * @brief Bounds checked version of operator [].
     */
    decltype(auto) operator[](std::size_t i) {
      ASSERT(i < size(), "Out of bounds: {} !< {}", i, size());
      return Vector::operator[](i);
    }

    /**
     * @brief Bounds checked version of operator [] const.
     */
    decltype(auto) operator[](std::size_t i) const {
      ASSERT(i < size(), "Out of bounds: {} !< {}", i, size());
      return Vector::operator[](i);
    }

    /**
     * @brief Provides an emplace_back with explicit ``matrix_t`` parameters.
     */
    decltype(auto) emplace_back(typename Mems::matrix_t const &...args) { return Vector::emplace_back(args...); }

    /**
     * @brief Provides an emplace_back with explicit ``matrix_t`` parameters.
     */
    decltype(auto) emplace_back(typename Mems::matrix_t &&...args) { return Vector::emplace_back(std::move(args)...); }
  };

}  // namespace fly::system