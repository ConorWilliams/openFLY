#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

/**
 * \file orthorhombic.hpp
 *
 * @brief orthorhombic simulation box implementation.
 */

#include <Eigen/Core>
#include <optional>

#include "libfly/system/atom.hpp"
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

namespace fly::system {

  /**
   * @brief Maps position vector tuples to a 1D grid.
   *
   * The OrthorhombicGrid is responsible for managing the mapping from ND->1D for the neighbour::List. The supercell is cut into a grid
   * of cells and each cell has a unique index. Atoms are assigned to a cell (at least as large in each dimension as r_cut).
   *
   */
  class OrthorhombicGrid {
  public:
    /**
     * @brief Get the total number of cells in the OrthorhombicGrid along each axis.
     */
    Arr<int> shape() const noexcept { return m_shape; }

    /**
     * @brief Compute the cell index from an atoms canonical position.
     *
     * This function only supports atoms in the canonical cell or atoms in the surrounding one grid-cell buffer i.e. ghost atoms.
     */
    template <typename E>
    int cell_idx(Eigen::MatrixBase<E> const& x) const noexcept {
      return to_1D(clamp_to_grid_idxs(x + m_cell.matrix()));
    }

  private:
    friend class Orthorhombic;

    Arr<int> m_shape = Arr<int>::Zero();
    Arr<int> m_prod_shape = Arr<int>::Zero();

    Arr<Position::scalar_t> m_cell = Arr<Position::scalar_t>::Zero();
    Arr<Position::scalar_t> m_inv_cell = Arr<Position::scalar_t>::Zero();

    OrthorhombicGrid(Arr<Position::scalar_t> const& extents, Position::scalar_t r_cut) {
      // Sanity checks
      VERIFY(r_cut > 0, "r_cut is negative");
      VERIFY((extents > r_cut).all(), "r_cut is too big");
      //
      m_shape = 2 + (extents / r_cut).cast<int>();
      m_cell = extents / (extents / r_cut).floor();
      m_inv_cell = 1.0 / m_cell;

      // Cumulative product of m_shape

      m_prod_shape = Arr<int>::Ones();

      for (int i = 1; i < spatial_dims; i++) {
        m_prod_shape[i] = m_prod_shape[i - 1] * m_shape[i - 1];
      }
    }

    /**
     * @brief Cast a position to its grid indexes and clamp on interval [0, m_shape -1].
     */
    template <typename E>
    Arr<int> clamp_to_grid_idxs(Eigen::MatrixBase<E> const& x) const noexcept {
      return (x.array() * m_inv_cell).template cast<int>().cwiseMax(0).cwiseMin(m_shape - 1);
    }

    /**
     * @brief Get the 1D index from the nD index.
     */
    int to_1D(Arr<int> const& indexes) const noexcept { return (indexes * m_prod_shape).sum(); }
  };

  /**
   * @brief Provides details of the simulations orthorhombic supercell geometry.
   *
   * All queries of the space in which the atoms exist are provided by this class. It is assumed
   * (and must be ensured) all non-periodic atoms are within the Orthorhombic extents.
   *
   * \rst
   * .. todo::
   *    Make test for this class work in ND.
   * \endrst
   */
  class Orthorhombic {
  public:
    /**
     * @brief Construct a new Orthorhombic box object.
     *
     * @param ex Length of simulation box along each axis.
     * @param pd True for each periodic axis.
     */
    Orthorhombic(Arr<Position::scalar_t> const& ex, Arr<bool> const& pd) : m_extents{ex}, m_periodic{pd}, m_inv_extents(1.0 / ex) {
      VERIFY((m_extents > 0).all(), "Orthorhombic extents are negative");
    }

    /**
     * @brief Fetch the basis vectors of this cell.
     *
     * @return Eigen::DiagonalMatrix<Position::scalar_t, spatial_dims> A matrix with each column corresponding to a basis vector.
     */
    Eigen::DiagonalMatrix<Position::scalar_t, spatial_dims> basis() const noexcept {
      return Eigen::DiagonalMatrix<Position::scalar_t, spatial_dims>{m_extents.matrix()};
    }

    /**
     * @brief Get an array detailing the periodicity along each axis.
     */
    Arr<bool> const& periodic() const noexcept { return m_periodic; }

    /**
     * @brief Generate the periodic image of an atom.
     *
     * @param x Position of the atom who's image we are computing.
     * @param ax Axis along which to generate image (magnitude = index, sign = direction).
     * @param r_cut Cut-off radius for atomic interactions.
     * @return std::optional<Position::matrix_t> If the atom's image is beyond the cut-off (``r_cut``) for atomic interactions then
     * the image's position otherwise, ``std::nullopt``.
     */
    template <typename A, typename B>
    std::optional<Position::matrix_t> gen_image(Eigen::MatrixBase<B> const& x, int ax, Position::scalar_t r_cut) {
      x + ax + r_cut;
      return {};
    }

    /**
     * @brief Maps atom into the canonical cell.
     *
     * \rst
     *
     * Guarantees:
     *
     * .. math::
     *    0 \le x_i < \text{extent}_i
     *
     * for all :math:`i`.
     *
     * \endrst
     */
    template <typename E>
    Position::matrix_t canon_image(Eigen::MatrixBase<E> const& x) const {
      // Non-periodic atoms are within the simulation box extents so x[i] * inv_extents less than 1 and x[i]
      // remains unaffected, hence no non-periodic switch/select.
      ASSERT((m_periodic || (x.array() >= Arr<Position::scalar_t>::Zero() && x.array() < m_extents)).all(), "Out of box");
      return x.array() - m_extents * (x.array() * m_inv_extents).floor();
    }

    /**
     * @brief Compute shortest vector connecting ``a`` to a periodic image of ``b``.
     *
     * This function is branchy and should be avoided in hot code.
     *
     * See: https://doi.org/10.1524/zpch.2013.0311
     */
    template <typename A, typename B>
    Position::matrix_t min_image(Eigen::MatrixBase<A> const& a, Eigen::MatrixBase<B> const& b) const noexcept {
      Arr<Position::scalar_t> dr = b - a;
      return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
    }

    /**
     * @brief Comparison operator, no surprises.
     */
    friend bool operator==(Orthorhombic const& a, Orthorhombic const& b) noexcept {
      return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
    }

    /**
     * @brief Make an OrthorhombicGrid object.
     *
     * @param r_cut The cut-off radius for atomic interactions.
     */
    OrthorhombicGrid make_grid(Position::scalar_t r_cut) const { return {m_extents, r_cut}; }

  private:
    Arr<Position::scalar_t> m_extents = Arr<Position::scalar_t>::Zero();
    Arr<bool> m_periodic = Arr<bool>::Zero();
    Arr<Position::scalar_t> m_inv_extents = Arr<Position::scalar_t>::Zero();
  };

}  // namespace fly::system