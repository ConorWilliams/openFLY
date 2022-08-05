#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: MPL-2.0

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Core>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "libfly/system/adaptor.hpp"
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file SoA.hpp
 *
 * @brief Struct of Arrays implementation
 */

namespace fly::system {

  /**
   * @brief Forward declaration
   */
  template <typename... Ms> class SoA;

  namespace detail {

    template <typename, typename> struct different_SoA : std::false_type {};

    template <typename... Ms> struct different_SoA<SoA<Ms...>, SoA<Ms...>> : std::false_type {};

    template <typename... Ms, typename... Mx> struct different_SoA<SoA<Ms...>, SoA<Mx...>> : std::true_type {};

  }  // namespace detail

  /**
   * @brief A container for atoms that stores each member in its own array.
   *
   * \rst
   *
   * The default container type used in the libfly; an ``SoA`` models an array of ``Atom`` types
   * but decomposes the atom type and stores each member in a separate array. This enables efficient
   * cache use. The members of the atom" are described through a series of template parameters
   * which should inherit from ``fly::MemTag``. A selection of canonical members are provided in the
   * namespace ``builtin_members``. The members of each atom can be accessed by the index of the
   * atom or as an Eigen array to enable collective operations.
   *
   * \endrst
   *
   * Example of use:
   *
   * @code{.cpp}
   *
   * #include "libatom/atom.hpp"
   *
   * using namespace otf;
   *
   * AtomArray<Position, AtomicNum> atoms{10}; // Initialise an array of 10 atoms.
   *
   * // Add a hydrogen atom at the origin.
   *
   * atoms(Position{}, 0) =  Vec3<double>{0, 0, 0};
   * atoms(AtomicNum{}, 0) = 1;
   *
   * Vec3 xyz = atoms(Position{}, 0); // Get the position of the zeroth atom.
   *
   * std::size_t n = = atoms(AtomicNum{}, 0); // Get the atomic number of the zeroth atom.
   *
   * atoms(Position{}) += 1; // Add 1 to each of every atoms coordinates.
   *
   * atoms(AtomicNum{}) = 6; // Set all the atoms to carbon atoms.
   *
   * @endcode
   */
  template <typename... Ms> class SoA : private detail::Adaptor<Ms>... {
  private:
    static constexpr bool owns_none = (std::is_reference_v<Ms> && ...);
    static constexpr bool owns_all = (!std::is_reference_v<Ms> && ...);

    template <typename T> static constexpr bool different_SoA_v = detail::different_SoA<SoA<Ms...>, remove_cref_t<T>>::value;

  public:
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
     * @brief Construct a new SoA conatining 'size' uninitialized atoms.
     *
     * Only SFINE enabled if this SoA owns all its arrays.
     */
    template <bool OwnsAll = owns_all> explicit SoA(int size, std::enable_if_t<OwnsAll>* = 0)
        : detail::Adaptor<Ms>(size)..., m_size(size) {}

    /**
     * @brief Implicitly construct a new SoA object from SoA 'other' with different members.
     *
     * Only SFINE enabled if this SoA owns non of its arrays.
     *
     * .. note::
     *    The implementation may ``std::move`` ``other`` multiple times but this is ok as the detail::Adaptor constructor will only
     *    move its corresponding base slice.
     *
     */
    template <typename T, typename = std::enable_if_t<different_SoA_v<T> && owns_none>> SoA(T&& other)
        : detail::Adaptor<Ms>(std::forward<T>(other))..., m_size(other.size()) {}

    /**
     * @brief Explicitly construct a new SoA object from SoA 'other' with different members.
     *
     * SFINE enabled if this SoA owns some of its arrays.
     *
     * .. note::
     *    The implementation may ``std::move`` ``other`` multiple times but this is ok as the detail::Adaptor constructor will only
     *    move its corresponding base slice.
     */
    template <typename T, typename = std::enable_if_t<different_SoA_v<T> && !owns_none>, typename = void> explicit SoA(T&& other)
        : detail::Adaptor<Ms>(std::forward<T>(other))..., m_size(other.size()) {}

    /**
     * @brief Defaulted copy assignment operator.
     *
     * \rst
     * .. warning::
     *    Assignement to reference members follows pointer-semantics.
     * \endrst
     */
    SoA& operator=(SoA const&) = default;

    /**
     * @brief Defaulted move assignment operator.
     *
     * \rst
     * .. warning::
     *    Assignement to reference members follows pointer-semantics.
     * \endrst
     */
    SoA& operator=(SoA&&) = default;

    /**
     * @brief Assign to a SoA with with different members.
     *
     * Assignement to reference members follows pointer-semantics.
     */
    template <typename T, typename = std::enable_if_t<different_SoA_v<T>>> SoA& operator=(T&& other) {
      (static_cast<void>(static_cast<detail::Adaptor<Ms>&>(*this) = std::forward<T>(other)), ...);
      return *this;
    }

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