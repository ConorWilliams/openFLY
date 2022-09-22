
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
#include <string>
#include <type_traits>
#include <utility>

#include "libfly/utility/core.hpp"

/**
 * \file property.hpp
 *
 * @brief The basic building block for atoms in libFLY.
 *
 * Properties in libFLY describe an atomic quantity i.e. mass, position, etc. They take the form of fixed-size matrices and must
 * always be derived from the Property class. The most fundamental property that every atom must have is a fly::TypeID. Properties can
 * either be associated per-atom using the fly::system::Atom, fly::system::VoS and fly::system::SoA data-structures or per-type using
 * the fly::system::TypeMap.
 */

namespace fly::system {

  /**
   * @brief A base type to derive from for defining properties.
   *
   * Properties may be fixed-size matrices/arrays of arithmetic types or scalar (1x1) default-constructible types.
   *
   * 1x1 matrices are unwrapped into scalars.
   *
   * A selection of canonical properties , deriving from this type, are provided in the inline namespace ``builtin_properties``.
   *
   * @tparam Scalar This property represents a matrix of Scalar elements.
   * @tparam Rows Number of rows in this property
   * @tparam Cols Number of colums in this property
   * @tparam Rep The Eigen3 template, Eigen::[matrix||array], to use for this property.
   */
  template <typename Scalar, int Rows = 1, int Cols = 1, template <typename, auto...> typename Rep = Eigen::Matrix>
  struct Property {
    /** @brief True if this property represents a 1x1 matrix. */
    static constexpr bool is_1x1 = Rows == 1 && Cols == 1;

    static_assert(std::is_arithmetic_v<Scalar> || (Rows == 1 && Cols == 1), "Non-scalar properties must be arithmetic.");

    static_assert(std::is_default_constructible_v<Scalar>, "Scalar properties must be default constructable.");

    static_assert(Rows > 0 && Cols > 0, "Invalid property extents.");

    /** @brief This property represents a matrix of elements of scalar_t. */
    using scalar_t = Scalar;

    /** @brief The matrix type that this property represents (1x1 matrices are unwrapped to scalars). */
    using matrix_t = std::conditional_t<is_1x1, Scalar, Rep<Scalar, Rows, Cols>>;
    /** @brief A reference-like type to a matrix_t. */
    using matrix_ref_t = std::conditional_t<is_1x1, Scalar&, Eigen::Map<matrix_t>>;
    /** @brief A const-reference-like type to a matrix_t. */
    using matrix_cref_t = std::conditional_t<is_1x1, Scalar const&, Eigen::Map<matrix_t const>>;

    /** @brief The Eigen type used to store a dynamic collection of contiguous matrix_t. */
    using array_t = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
    /** @brief A reference-like type to the underlying array_t. */
    using array_ref_t = Eigen::Map<array_t>;
    /** @brief A const-reference-like type to the underlying array_t. */
    using array_cref_t = array_t const&;

    /**
     * @brief Get the number of elements in the matrix_t.
     */
    static constexpr int size() { return Rows * Cols; }
  };

}  // namespace fly::system

namespace fly {

  /**
   * @brief An inline namespace providing a selection of canonical properties  for Atom.
   */
  inline namespace builtin_properties {

    // GSD Schema

    /**
     * @brief Tag type for atom's type.
     *
     * Types are abstract labels assigned to atoms, if two atoms have the same type they are completely interchangeable in every way.
     */
    struct Type : system::Property<char, 32, 1, Eigen::Array> {
      static constexpr char const* tag = "particles/types";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's id.
     *
     * Each atom has a Type, an atoms TypeID is a number between 0 and N - 1 with N the number of types in the simulation. A TypeMap
     * object performs the transformation from TypeId -> Type.
     */
    struct TypeID : system::Property<std::uint32_t> {
      static constexpr char const* tag = "particles/typeid";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's image count, i.e. how many times atom has wrapped.
     *
     * Not used by libFLY, for consistency with HOOMD schema.
     */
    struct Image : system::Property<std::int32_t, spatial_dims> {
      static constexpr char const* tag = "particles/image";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's mass.
     */
    struct Mass : system::Property<double> {
      static constexpr char const* tag = "particles/mass";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's charge.
     */
    struct Charge : system::Property<double> {
      static constexpr char const* tag = "particles/charge";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's effective diameter.
     */
    struct Diameter : system::Property<double> {
      static constexpr char const* tag = "particles/diameter";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's moment of inertia in an atoms body frame.
     */
    struct MomentInertia : system::Property<double> {
      static constexpr char const* tag = "particles/moment_inertia";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's position in real space.
     */
    struct Position : system::Property<double, spatial_dims> {
      static constexpr char const* tag = "particles/position";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's velocity in real space.
     */
    struct Velocity : system::Property<double, spatial_dims> {
      static constexpr char const* tag = "particles/velocity";  ///< GSD chunk label.
    };

    // Custom

    /**
     * @brief Tag type for atom's index i.e. position in some array.
     */
    struct Index : system::Property<Eigen::Index> {
      static constexpr char const* tag = "log/particles/index";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's status as a frozen atom.
     *
     * Note this cannot be written to a GSD file as it a boolean.
     */
    struct Frozen : system::Property<bool> {
      static constexpr char const* tag = "log/particles/frozen";  ///< Custom chunk label.
    };

    /**
     * @brief Tag type for atom's contribution to a dimer axis.
     */
    struct Axis : system::Property<double, spatial_dims> {
      static constexpr char const* tag = "log/particles/axis";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's potential gradient acting on an atom.
     */
    struct PotentialGradient : system::Property<double, spatial_dims> {
      static constexpr char const* tag = "log/particles/potential_gradient";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for a changes in some atomic property.
     */
    struct Delta : system::Property<double, spatial_dims> {
      static constexpr char const* tag = "log/particles/delta";  ///< GSD chunk label.
    };

    /**
     * @brief Tag type for atom's acceleration.
     */
    struct Acceleration : system::Property<double, spatial_dims> {
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
     * @brief Delta literal.
     */
    inline constexpr Delta del_;

    /**
     * @brief Acceleration literal.
     */
    inline constexpr Acceleration a_;

  }  // namespace builtin_properties

}  // namespace fly