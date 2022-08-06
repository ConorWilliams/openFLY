#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * \file orthorhombic.hpp
 *
 * @brief orthorhombic simulation box implementation.
 */

#include <Eigen/Core>

#include "libfly/system/atom.hpp"
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

namespace fly::system {
  /**
   * @brief Provides details of the simulations orthorhombic supercell geometry.
   *
   * All queries of the space in which the atoms exist are provided by this class. It is assumed
   * (and must be ensured) all non-periodic atoms are within the Orthorhombic extents.
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
     * The result is in the form of an ``Eigen::DiagonalMatrix`` with each column corresponding to a basis vector.
     */
    Eigen::DiagonalMatrix<floating, spatial_dims> basis() const noexcept { return {m_extents[0], m_extents[1], m_extents[2]}; }

    /**
     * @brief Get an array detailing the periodicity along each axis.
     */
    Arr<bool> const& periodic() const noexcept { return m_periodic; }

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
    Vec<floating> canon_image(Vec<floating> const& x) const noexcept {
      // Non-periodic atoms are within the simulation box extents so x[i] * inv_extents less than 1 and x[i]
      // remains unaffected, hence no non-periodic switch/select.
      ASSERT((m_periodic || (x.array() >= Arr<floating>::Zero() && x.array() < m_extents)).all(), "Out of box");
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
    [[deprecated]] Vec<floating> min_image(Eigen::MatrixBase<A> const& a, Eigen::MatrixBase<B> const& b) const noexcept {
      Arr<floating> dr = b - a;
      return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
    }

    /**
     * @brief Comparison operator, no surprises.
     */
    friend bool operator==(Orthorhombic const& a, Orthorhombic const& b) noexcept {
      return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
    }

    /**
     * @brief Maps position vector tuples to a 1D grid.
     *
     * The Grid is responsible for managing the mapping from ND->1D for the neighbour::List. The supercell is cut into a grid of cells
     * and each cell has a unique index. Atoms are assigned to a cell (at least as large in each dimension as r_cut).
     *
     */
    class Grid {
    public:
      /**
       * @brief Get the total number of cells in the Grid along each axis.
       */
      Arr<int> shape() const noexcept { return m_shape; }

      /**
       * @brief Get a vector which displaces an atom by one grid cell along each axis.
       */
      Vec<floating> cell_offset() const noexcept { return m_cell; }

      /**
       * @brief Compute the cell index from an atoms position.
       *
       * Assumes atom is in the canonical cell + displaced by cell_offset().
       */
      int cell_idx(Vec<floating> const& x) const { return to_1D(clamp_to_grid_idxs(x)); }

    private:
      friend Orthorhombic;

      Arr<int> m_shape = Vec<int>::Zero();
      Arr<int> m_prod_shape = Vec<int>::Zero();

      Arr<floating> m_cell = Vec<floating>::Zero();
      Arr<floating> m_inv_cell = Vec<floating>::Zero();

      Grid(Arr<floating> const& extents, floating rcut) {
        m_shape = 2 + (extents / rcut).cast<int>();
        m_cell = extents / (extents / rcut).floor();
        m_inv_cell = 1.0 / m_cell;

        // Sanity checks
        VERIFY(rcut > 0, "r_cut is negative");
        VERIFY((extents > rcut).all(), "r_cut is too big");

        // Cumulative product of m_shape

        m_prod_shape = Arr<int>::Ones();

        for (int i = 1; i < spatial_dims; i++) {
          m_prod_shape[i] = m_prod_shape[i - 1] * m_shape[i - 1];
        }
      }

      /**
       * @brief Cast a Vec<floating> to Vec<int> and clamp on interval [0, m_shape -1].
       */
      Arr<int> clamp_to_grid_idxs(Vec<floating> const& x) const& {
        return (x.array() * m_inv_cell).cast<int>().cwiseMax(0).cwiseMin(m_shape - 1);
      }

      /**
       * @brief Get the 1D index from the nD index.
       */
      int to_1D(Arr<int> const& indexes) const& { return (indexes * m_prod_shape).sum(); }
    };

    /**
     * @brief Make a grid
     *
     * @param rcut
     * @return Grid
     */
    Grid make_grid(floating rcut) const { return {m_extents, rcut}; }

  private:
    Arr<floating> m_extents = Arr<floating>::Zero();
    Arr<bool> m_periodic = Arr<bool>::Zero();
    Arr<floating> m_inv_extents = Arr<floating>::Zero();
  };

}  // namespace fly::system