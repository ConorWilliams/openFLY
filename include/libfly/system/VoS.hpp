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

#include "libfly/system/atom.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file VoS.hpp
 */

namespace fly::system {

  /**
   * @brief A container that models a std::vector of atoms.
   *
   * The properties  of the "atom" are described through a series of template parameters which should
   * inherit from ``Property``
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
   * @tparam T a series of empty types, derived from ``Property``, to describe each member.
   */
  template <typename... T>
  class VoS : public Vector<Atom<T...>> {
  private:
    static_assert(sizeof...(T) > 0, "Need at least one property in an VoS");

    using vector = Vector<Atom<T...>>;

  public:
    /**
     * @brief The underlying atom type.
     */
    using atom_t = Atom<T...>;

    using vector::emplace_back;

    /**
     * @brief Provides an emplace_back with explicit ``matrix_t const&`` parameters.
     */
    decltype(auto) emplace_back(typename T::matrix_t const &...args) { return vector::emplace_back(args...); }

    /**
     * @brief Provides an emplace_back with explicit ``matrix_t &&`` parameters.
     */
    decltype(auto) emplace_back(typename T::matrix_t &&...args) { return vector::emplace_back(std::move(args)...); }

    /**
     * @brief Lib cereal serialization support.
     */
    template <class Archive>
    void serialize(Archive &archive) {
      archive(cereal::base_class<Vector<Atom<T...>>>(this));
    }
  };

}  // namespace fly::system