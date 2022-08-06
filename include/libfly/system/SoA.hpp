#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <Eigen/Core>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "libfly/system/adaptor.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file SoA.hpp
 *
 * @brief Struct of Arrays implementation.
 */

namespace fly::system {

  template <typename... Ms> class SoA;  // Required for detail below.

  namespace detail {

    template <typename, typename> struct different_SoA : std::false_type {};

    template <typename... Ms> struct different_SoA<SoA<Ms...>, SoA<Ms...>> : std::false_type {};

    template <typename... Ms, typename... Mx> struct different_SoA<SoA<Ms...>, SoA<Mx...>> : std::true_type {};

  }  // namespace detail

  /**
   * @brief A container for atoms that stores each member in its own array.
   *
   * The default container type used in the libfly; an ``SoA`` models an array of ``Atom`` types
   * but decomposes the atom type and stores each member in a separate array. This enables efficient
   * cache use. Like ``Atom``, the members of the "atom" are described through a series of template parameters
   * which should inherit from ``MemTag``. The members of each atom can be accessed by the index of the
   * atom or as an ``Eigen::Array`` to enable collective operations.
   *
   * SoA also supports slicing and reference members which transform that member into a view, this enables SoA to act as a concrete
   * type in interfaces.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/system/SoA.cpp
   *    :code:
   *
   * \endrst
   *
   * @tparam Ms a series of empty types, derived from ``MemTag``, to describe each member.
   */
  template <typename... Ms> class SoA : private detail::Adaptor<Ms>... {
  public:
    /**
     * @brief True if this SoA is a pure view i.e. all its members are reference types.
     */
    static constexpr bool owns_none = (std::is_reference_v<Ms> && ...);

    /**
     * @brief True if this SoA is purely owning i.e. non of its members are reference types.
     */
    static constexpr bool owns_all = (!std::is_reference_v<Ms> && ...);

    /**
     * @brief Detect if a type is a specialization of a SoA but different from this specialization.
     */
    template <typename T> static constexpr bool different_SoA_v = detail::different_SoA<SoA<Ms...>, remove_cref_t<T>>::value;

    /**
     * @brief Construct a new empty SoA.
     */
    SoA() = default;

    /**
     * @brief Defaulted move constructor.
     */
    SoA(SoA&&) = default;

    /**
     * @brief Defaulted copy constructor.
     */
    SoA(SoA const&) = default;

    /**
     * @brief Construct a new SoA conatining 'size' atoms.
     *
     * \rst
     *
     * Only SFINE enabled if this SoA owns all its arrays.
     *
     * .. warning::
     *    New atoms are uninitialized.
     * \endrst
     */
    template <bool OwnsAll = owns_all> explicit SoA(int size, std::enable_if_t<OwnsAll>* = 0)
        : detail::Adaptor<Ms>(size)..., m_size(size) {}

    /**
     * @brief Implicitly construct a new SoA object from SoA 'other' with different members.
     *
     * Only SFINE enabled if this SoA owns non of its arrays.
     *
     */
    template <typename T, typename = std::enable_if_t<different_SoA_v<T> && owns_none>> SoA(T&& other)
        : detail::Adaptor<Ms>(std::forward<T>(other))..., m_size(other.size()) {
      // Ok to ``std::move`` ``other`` multiple times but this is ok as the detail::Adaptor constructor will only move its
      // corresponding base slice.
    }

    /**
     * @brief Explicitly construct a new SoA object from SoA 'other' with different members.
     *
     * SFINE enabled if this SoA owns some of its arrays.
     *
     */
    template <typename T, typename = std::enable_if_t<different_SoA_v<T> && !owns_none>, typename = void> explicit SoA(T&& other)
        : detail::Adaptor<Ms>(std::forward<T>(other))..., m_size(other.size()) {
      // Ok to ``std::move`` ``other`` multiple times but this is ok as the detail::Adaptor constructor will only move its
      // corresponding base slice.
    }

    /**
     * @brief Defaulted copy assignment operator.
     *
     * \rst
     *
     * See the `note on assignment to references`_.
     *
     * \endrst
     */
    SoA& operator=(SoA const&) = default;

    /**
     * @brief Defaulted move assignment operator.
     *
     * \rst
     *
     * See the `note on assignment to references`_.
     *
     * \endrst
     */
    SoA& operator=(SoA&&) = default;

    /**
     * @brief Assign to a SoA with with different members.
     *
     * \rst
     *
     * .. _`note on assignment to references`:
     *
     * .. note::
     *    Assignement to reference members follows pointer-semantics:
     *
     *    .. include:: ../../examples/system/SoA_assign.cpp
     *       :code:
     *
     * \endrst
     */
    template <typename T, typename = std::enable_if_t<different_SoA_v<T>>> SoA& operator=(T&& other) {
      (static_cast<void>(static_cast<detail::Adaptor<Ms>&>(*this) = std::forward<T>(other)), ...);
      return *this;
    }

    // Inherit operators

    using detail::Adaptor<Ms>::operator()...;

    using detail::Adaptor<Ms>::operator[]...;

    /**
     * @brief Get the number of atoms in the SoA.
     */
    int size() const noexcept { return m_size; }

    /**
     * @brief Resize the SoA.
     *
     * \rst
     * .. warning::
     *    Destroys all the data in the SoA, new atoms are all uninitialized.
     * \endrst
     *
     */
    template <bool OwnsAll = owns_all> void destructive_resize(int new_size, std::enable_if_t<OwnsAll>* = 0) {
      if (std::exchange(m_size, new_size) != new_size) {
        (static_cast<void>(get(Ms{}).resize(new_size * Ms::size(), Eigen::NoChange)), ...);
      }
    }

  private:
    int m_size = 0;

    /* clang-format off */ template <typename...>  friend class SoA; /* clang-format on */
  };

}  // namespace fly::system