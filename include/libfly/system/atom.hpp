
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

  /**
   * @brief A base type to derive from for defining members of an Atom type.
   *
   * Members must be matrices of arithmetic types or default constructible 1x1 matricies.
   *
   * 1x1 matricies are unwrapped into scalars.
   *
   * A selection of canonical members, deriving from this type, are provided in the namespace ``builtin_m``.
   *
   * @tparam Scalar This member represents a matrix of Scalar elements.
   * @tparam Rows Number of rows in this member.
   * @tparam Cols Number of colums in this member.
   * @tparam Rep The Eigen3 template, Eigen::[matrix||array], to use for this member.
   */
  template <typename Scalar, int Rows = 1, int Cols = 1, template <typename, auto...> typename Rep = Eigen::Matrix>
  struct MemTag {
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

    template <typename Tag>
    struct AtomMem {
    public:
      static_assert(std::is_empty_v<Tag>, "Member tag types are required to be empty");

      AtomMem(AtomMem&&) = default;

      AtomMem(AtomMem const&) = default;

      explicit AtomMem(typename Tag::matrix_t&& data) : m_data(std::move(data)) {}

      explicit AtomMem(typename Tag::matrix_t const& data) : m_data(data) {}

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
   * Libfly uses this type to build atoms to allow integration with libfly's containers.
   * An Atom behaves as a struct of its members, which are accessed through ``operator[]``.
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
   * @tparam Mems a series of empty types, derived from ``MemTag``, to describe each member.
   */
  template <typename... Mems>
  struct Atom : detail::AtomMem<Mems>... {
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

    using detail::AtomMem<Mems>::operator[]...;
  };

}  // namespace fly::system

namespace fly {

  /**
   * @brief An inline namespace providing a selection of canonical members for Atom.
   */
  inline namespace builtin_m {

    // GSD Schema

    /**
     * @brief Tag type for atom's type.
     *
     * Types are abstract labels assigned to atoms, if two atoms have the same type they are completely interchangeable in every way.
     */
    struct Type : system::MemTag<std::string> {
      static constexpr char const* tag = "particles/types";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's id.
     *
     * Each atom has a Type, an atoms TypeID is a number between 0 and N - 1 with N the number of types in the simulation. A TypeMap
     * object performs the transformation from TypeId -> Type.
     */
    struct TypeID : system::MemTag<std::uint32_t> {
      static constexpr char const* tag = "particles/typeid";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's image count, i.e. how many times atom has wrapped.
     *
     * Not used by libFLY, for consistency with HOOMD schema.
     */
    struct Image : system::MemTag<std::int32_t, spatial_dims> {
      static constexpr char const* tag = "particles/image";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's mass.
     */
    struct Mass : system::MemTag<double> {
      static constexpr char const* tag = "particles/mass";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's charge.
     */
    struct Charge : system::MemTag<double> {
      static constexpr char const* tag = "particles/charge";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's effective diameter.
     */
    struct Diameter : system::MemTag<double> {
      static constexpr char const* tag = "particles/diameter";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's moment of inertia in an atoms body frame.
     */
    struct MomentInertia : system::MemTag<double> {
      static constexpr char const* tag = "particles/moment_inertia";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's position in real space.
     */
    struct Position : system::MemTag<double, spatial_dims> {
      static constexpr char const* tag = "particles/position";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's velocity in real space.
     */
    struct Velocity : system::MemTag<double, spatial_dims> {
      static constexpr char const* tag = "particles/velocity";  ///< GSD chunk label.
    };

    // Custom

    /**
     * @brief Tag type for atom's index i.e. position in some array.
     */
    struct Index : system::MemTag<std::uint32_t> {
      static constexpr char const* tag = "log/particles/index";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's status as a frozen atom.
     *
     * Note this cannot be written to a GSD file as it a boolean.
     */
    struct Frozen : system::MemTag<bool> {};

    /**
     * @brief Tag type for atom's contribution to a dimer axis.
     */
    struct Axis : system::MemTag<double, spatial_dims> {
      static constexpr char const* tag = "log/particles/axis";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's potential gradient acting on an atom.
     */
    struct PotentialGradient : system::MemTag<double, spatial_dims> {
      static constexpr char const* tag = "log/particles/potential_gradient";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's acceleration.
     */
    struct Acceleration : system::MemTag<double, spatial_dims> {
      static constexpr char const* tag = "log/particles/acceleration";  ///< GSD chunk label.
    };

    /**
     * @brief Type literal.
     */
    inline constexpr Type tp_;

    /**
     * @brief TypeID literal.
     */
    inline constexpr TypeID id_;

    /**
     * @brief Image literal.
     */
    inline constexpr Image img_;

    /**
     * @brief Mass literal.
     */
    inline constexpr Mass m_;

    /**
     * @brief Charge literal.
     */
    inline constexpr Charge q_;

    /**
     * @brief Diameter literal.
     */
    inline constexpr Diameter d_;
    /**
     * @brief MomentInertia literal.
     */
    inline constexpr MomentInertia I_;
    /**
     * @brief Position literal.
     */
    inline constexpr Position r_;

    /**
     * @brief Velocity literal.
     */
    inline constexpr Velocity v_;

    /**
     * @brief Index literal.
     */
    inline constexpr Index i_;

    /**
     * @brief Frozen literal.
     */
    inline constexpr Frozen fzn_;

    /**
     * @brief Axis literal.
     */
    inline constexpr Axis ax_;

    /**
     * @brief PotentialGradient literal.
     */
    inline constexpr PotentialGradient g_;

    /**
     * @brief Acceleration literal.
     */
    inline constexpr Acceleration a_;

  }  // namespace builtin_m

}  // namespace fly