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
   * @brief A base type to derive from for defining members of an Atom type.
   *
   * Members must be matrices of arithmetic types or default constructible 1x1 matricies.
   *
   * 1x1 matricies are unwrapped into scalars
   *
   * @tparam Scalar This member represents a matrix of Scalar elements.
   * @tparam Rows Number of rows in this member.
   * @tparam Cols Number of colums in this member.
   * @tparam Rep The Eigen3 template, Eigen::[matrix||array], to use for this member.
   */
  template <typename Scalar, int Rows = 1, int Cols = 1, template <typename, auto...> typename Rep = Eigen::Matrix> struct MemTag {
    /** @brief True if this member represents a 1x1 matrix. */
    static constexpr bool is_1x1 = Rows == 1 && Cols == 1;

    static_assert(std::is_arithmetic_v<Scalar> || (Rows == 1 && Cols == 1), "Non-scalar members must be arithmetic.");

    static_assert(std::is_default_constructible_v<Scalar>, "Scalar members must be default constructable.");

    static_assert(Rows > 0 && Cols > 0, "Invalid member extents.");

    /** @brief This member represents a matrix of elements of scalar_t. */
    using scalar_t = Scalar;

    /** @brief The matrix type that this member represents (1x1 matrices are unwrapped to scalars). */
    using matrix_t = std::conditional_t<is_1x1, Scalar, Rep<Scalar, Rows, Cols>>;
    /** @brief A reference-like type to a matrix_t. */
    using matrix_ref_t = std::conditional_t<is_1x1, Scalar&, Eigen::Map<matrix_t>>;
    /** @brief A const-reference-like type to a const matrix_t. */
    using matrix_cref_t = std::conditional_t<is_1x1, Scalar const&, Eigen::Map<matrix_t const>>;

    /** @brief The Eigen type used to store a dynamic collection of contiguous matrix_t. */
    using array_t = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
    /** @brief A reference-like type to the underlying array_t. */
    using array_ref_t = Eigen::ArrayBase<array_t>&;
    /** @brief A const-reference-like type to the underlying array_t. */
    using array_cref_t = Eigen::ArrayBase<array_t> const&;

    /**
     * @brief Get the number of elements in the matrix_t.
     */
    static constexpr int size() { return Rows * Cols; }
  };

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
   * @brief A container that stores a Atom types decomposed by member.
   *
   * The default container type used in the libatom; an AtomArray models an array of "atom" types
   * but decomposes the atom type and stores each member in a separate array. This enables efficient
   * cache use. The members of the "atom" are described through a series of template parameters
   * which should inherit from otf::MemTag. A selection of canonical members are provided in the
   * namespace builtin_members. The members of each atom can be accessed either by the index of the
   * atom or as an Eigen array to enable collective operations.
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

    //

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
     *    Destroys all the data in the SoA and
     * \endrst
     *
     * @tparam OwnsAll
     * @param new_size
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