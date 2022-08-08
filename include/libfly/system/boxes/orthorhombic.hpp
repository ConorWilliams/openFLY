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
#include "libfly/utility/core.hpp"

/**
 * \file orthorhombic.hpp
 *
 * @brief Specialised simulation box for orthorhombic supercells.
 */

namespace fly::system {

  /**
   * @brief Generalised orthogonal simulation box.
   */
  class Orthorhombic {
  public:
    class Grid;

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
     *
     * @deprecated No triclinic generalization.
     */
    template <typename A, typename B>
    [[deprecated("No triclinic generalization")]] Position::matrix_t min_image(Eigen::MatrixBase<A> const& a,
                                                                               Eigen::MatrixBase<B> const& b) const noexcept {
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
     * @brief Make a Grid object.
     *
     * @param r_cut The cut-off radius for atomic interactions.
     */
    Grid make_grid(Position::scalar_t r_cut) const;  // Implementation after Grid at EoF

  private:
    Arr<Position::scalar_t> m_extents = Arr<Position::scalar_t>::Zero();
    Arr<bool> m_periodic = Arr<bool>::Zero();
    Arr<Position::scalar_t> m_inv_extents = Arr<Position::scalar_t>::Zero();
  };

  /**
   * @brief Maps position vector -> ND integer-tuples -> 1D index.
   *
   * Grid is constructable using the Orthorhombic::make_grid factory.
   */
  class Orthorhombic::Grid {
  public:
    /**
     * @brief Get the total number of cells in the Grid along each axis.
     */
    Arr<int> shape() const noexcept { return m_shape; }

    /**
     * @brief Compute the cell index from an atoms canonical position.
     *
     * \rst
     * .. warning::
     *
     *    Assumes atoms are in the canonical cell or in the surrounding one Grid-cell buffer i.e. ghost atoms.
     * \endrst
     */
    template <typename E>
    int cell_idx(Eigen::MatrixBase<E> const& x) const noexcept {
      return to_1D(clamp_to_grid_idxs(x + m_cell.matrix()));
    }

    /**
     * @brief Generate the image of an atom along a particular axis ``ax``.
     *
     * \rst
     * .. warning::
     *    Only supports atoms in the canonical cell.
     * \endrst
     *
     * @tparam S Direction along axis which to generate image.
     * @param x Position of the atom who's image we are computing.
     * @param ax Index of axis along which to generate image.
     * @return std::optional<Position::matrix_t> If the atom's image is beyond the cut-off (``r_cut``) for atomic interactions then
     * the image's position otherwise ``std::nullopt``.
     */
    template <Sign S>
    std::optional<Position::matrix_t> gen_image(Position::matrix_t x, int ax) {
      //
      static_assert(S == Sign::plus || S == Sign::minus, "Unreachable");

      if constexpr (S == Sign::plus) {
        if (x[ax] < m_r_cut) {
          x[ax] += m_extents[ax];
          return x;
        }
      } else {
        if (x[ax] > m_extents[ax] - m_r_cut) {
          x[ax] -= m_extents[ax];
          return x;
        }
      }

      return std::nullopt;
    }

  private:
    friend class Orthorhombic;

    Arr<Position::scalar_t> m_extents = Arr<Position::scalar_t>::Zero();

    Position::scalar_t m_r_cut = 0;

    Arr<int> m_shape = Arr<int>::Zero();
    Arr<int> m_prod_shape = Arr<int>::Zero();

    Arr<Position::scalar_t> m_cell = Arr<Position::scalar_t>::Zero();
    Arr<Position::scalar_t> m_inv_cell = Arr<Position::scalar_t>::Zero();

    // Private constructor
    Grid(Arr<Position::scalar_t> const& extents, Position::scalar_t r_cut) : m_extents(extents), m_r_cut(r_cut) {
      // Sanity checks
      VERIFY(r_cut > 0, "r_cut is negative");
      VERIFY((extents > r_cut).all(), "r_cut is too big");
      //
      m_shape = 2 + (extents / r_cut).cast<int>();
      m_cell = extents / (extents / r_cut).floor();
      m_inv_cell = 1.0 / m_cell;
      m_prod_shape = product_scan(m_shape);
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

  inline Orthorhombic::Grid Orthorhombic::make_grid(Position::scalar_t r_cut) const { return {m_extents, r_cut}; }

}  // namespace fly::system