
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
 * @brief Atom type and basic building blocks.
 */

namespace fly::system {

  /**
   * @brief A base type to derive from for defining members of an Atom type.
   *
   * Members must be matrices of arithmetic types or default constructible 1x1 matricies.
   *
   * 1x1 matricies are unwrapped into scalars.
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

    /**
     * @brief Get the number of rows in the matrix_t.
     */
    static constexpr int rows() { return Rows; }

    /**
     * @brief Get the number of cols in the matrix_t.
     */
    static constexpr int cols() { return Cols; }
  };

  namespace detail {

    template <typename Tag> struct AtomMem {
    public:
      static_assert(std::is_empty_v<Tag>, "Member tag types are required to be empty");

      AtomMem(AtomMem&&) = default;

      AtomMem(AtomMem const&) = default;

      explicit AtomMem(typename Tag::matrix_t&& data) : m_data(data) {}

      explicit AtomMem(typename Tag::matrix_t const& data) : m_data(data) {}

      template <typename T> explicit AtomMem(T&& x) : m_data(std::forward<T>(x)) {}

      typename Tag::matrix_t const& operator[](Tag) const { return m_data; }

      typename Tag::matrix_t& operator[](Tag) { return m_data; }

    private:
      typename Tag::matrix_t m_data;
    };

  }  // namespace detail

  /**
   * @brief The libfly representation of an atom.
   *
   * @tparam Mems a series of empty types, derived from otf::MemTag, to describe each member.
   */
  template <typename... Mems> struct Atom : detail::AtomMem<Mems>... {
    /**
     * @brief Copy construct a new Atom object
     */
    Atom(Atom&&) = default;

    /**
     * @brief Move construct a new Atom object
     */
    Atom(Atom const&) = default;

    /**
     * @brief Construct a new Atom object, initializing each member.
     */
    explicit Atom(typename Mems::matrix_t const&... args) : detail::AtomMem<Mems>(args)... {}

    /**
     * @brief Construct a new Atom object, initializing each member.
     */
    explicit Atom(typename Mems::matrix_t&&... args) : detail::AtomMem<Mems>(args)... {}

    // Expose tagged dispatch.
    using detail::AtomMem<Mems>::operator[]...;
  };

}  // namespace fly::system