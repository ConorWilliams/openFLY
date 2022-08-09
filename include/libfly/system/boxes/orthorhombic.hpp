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
#include <optional>

#include "libfly/system/atom.hpp"
#include "libfly/system/boxes/hypergrid.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file orthorhombic.hpp
 *
 * @brief Specialised simulation box for orthorhombic supercells.
 */

namespace fly::system {

  /**
   * @brief Maps position vector -> ND integer-tuples -> 1D index.
   *
   * OrthoGrid is constructable using the Orthorhombic::make_grid factory.
   */
  class OrthoGrid : public HyperGrid {
  public:
    /**
     * @brief Generate the periodic image of an atom along a particular axis ``ax``.
     *
     * \rst
     * .. warning::
     *    Only supports atoms in the canonical cell.
     *
     * .. note::
     *    Not a template for implementation reasons.
     * \endrst
     *
     * @tparam S Direction along axis which to generate image.
     * @param x Position of the atom who's image we are computing.
     * @param ax Index of axis along which to generate image.
     * @return std::optional<Vec<double>> If the atom's image is beyond the cut-off (``r_cut``) for atomic interactions then
     * the image's position otherwise ``std::nullopt``.
     */
    template <Sign S>
    std::optional<Vec<double>> gen_image(Vec<double> x, int ax) {
      //
      static_assert(S == Sign::plus || S == Sign::minus, "Unreachable");

      if constexpr (S == Sign::plus) {
        if (x[ax] < HyperGrid::r_cut()) {
          x[ax] += m_extents[ax];
          return x;
        }
      } else {
        if (x[ax] > m_extents[ax] - HyperGrid::r_cut()) {
          x[ax] -= m_extents[ax];
          return x;
        }
      }

      return std::nullopt;
    }

  private:
    friend class Orthorhombic;

    Arr<double> m_extents;

    /**
     * @brief Construct a new OrthoGrid object.
     *
     * @param extents The size of the box along each axis.
     * @param r_cut The cut-off radius for atomic interactions.
     */
    OrthoGrid(Arr<double> const& extents, double r_cut) : HyperGrid(extents, r_cut), m_extents(extents) {}
  };

  /**
   * @brief Generalised orthogonal simulation box.
   */
  class Orthorhombic {
  public:
    /**
     * @brief Construct a new Orthorhombic box object.
     *
     * @param ex Length of simulation box along each axis.
     * @param pd True for each periodic axis.
     */
    Orthorhombic(Arr<double> const& ex, Arr<bool> const& pd) : m_extents{ex}, m_periodic{pd}, m_inv_extents(1.0 / ex) {
      VERIFY((m_extents > 0).all(), "Orthorhombic extents are negative");
    }

    /**
     * @brief Fetch the basis vectors of this cell.
     *
     * @return Eigen::DiagonalMatrix<double, spatial_dims> A matrix with each column corresponding to a basis vector.
     */
    Eigen::DiagonalMatrix<double, spatial_dims> basis() const noexcept {
      return Eigen::DiagonalMatrix<double, spatial_dims>{m_extents.matrix()};
    }

    /**
     * @brief Query if the ``i``th axis is periodic.
     */
    bool periodic(int i) const noexcept { return m_periodic[i]; }

    /**
     * @brief Maps an atom into the canonical cell.
     *
     * \rst
     *
     * .. warning::
     *
     *    Assumes atoms are within the extents of the non-periodic axes!
     *
     * Guarantees :math:`0 \le x_i < \text{extent}_i` for all :math:`i`.
     *
     * \endrst
     */
    template <typename E>
    Vec<double> canon_image(Eigen::MatrixBase<E> const& x) const {
      // Non-periodic atoms are within the simulation box extents so x[i] * inv_extents less than 1 and x[i]
      // remains unaffected, hence no non-periodic switch/select.
      ASSERT((m_periodic || (x.array() >= Arr<double>::Zero() && x.array() < m_extents)).all(), "Out of box");
      return x.array() - m_extents * (x.array() * m_inv_extents).floor();
    }

    /**
     * @brief Compute shortest vector connecting ``a`` to a periodic image of ``b``.
     *
     * This function is branchy and should be avoided in hot code.
     *
     * See: https://doi.org/10.1524/zpch.2013.0311
     *
     * @deprecated No triclinic generalization.
     */
    template <typename A, typename B>
    [[deprecated("No triclinic generalization")]] Vec<double> min_image(Eigen::MatrixBase<A> const& a,
                                                                        Eigen::MatrixBase<B> const& b) const noexcept {
      Arr<double> dr = b - a;
      return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
    }

    /**
     * @brief Comparison operator, no surprises.
     */
    friend bool operator==(Orthorhombic const& a, Orthorhombic const& b) noexcept {
      return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
    }

    /**
     * @brief Make an OrthoGrid object.
     *
     * @param r_cut The cut-off radius for atomic interactions.
     */
    OrthoGrid make_grid(double r_cut) const { return {m_extents, r_cut}; }

  private:
    Arr<double> m_extents = Arr<double>::Zero();
    Arr<bool> m_periodic = Arr<bool>::Zero();
    Arr<double> m_inv_extents = Arr<double>::Zero();
  };

}  // namespace fly::system