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
 * \file hypergrid.hpp
 *
 * @brief Generalised utilities for sub-dividing hyperrectangles.
 */

namespace fly::system {

  /**
   * @brief Maps a position-vector -> ND-integer-tuples -> 1D-index.
   *
   * Dives a hyperrectangle into smaller hyperrectangle called cells.
   */
  class HyperGrid {
  public:
    /**
     * @brief Construct a new Grid object with cells at least `r_cut` along each axis.
     *
     * @param extents The size of the box along each axis.
     * @param r_cut The cut-off radius for atomic interactions.
     */
    HyperGrid(Arr<double> const& extents, double r_cut) : m_r_cut(r_cut) {
      // Sanity checks
      verify(r_cut > 0, "r_cut={} is negative", r_cut);
      verify((extents > r_cut).all(), "r_cut={} is too big for {}", r_cut, extents);

      //
      m_shape = 2 + (extents / r_cut).cast<int>();
      m_cell = extents / (extents / r_cut).floor();
      m_inv_cell = 1.0 / m_cell;
      m_prod_shape = product_scan(m_shape);
    }

    /**
     * @brief Fetch the cut-off radius for atomic interactions.
     */
    double r_cut() const noexcept { return m_r_cut; }

    /**
     * @brief Get the total number of cells in the HyperGrid along each axis.
     */
    Arr<int> shape() const noexcept { return m_shape; }

    /**
     * @brief Compute the cell index from an atoms canonical position.
     *
     * \rst
     * .. warning::
     *
     *    Assumes atoms are in the extents of the HyperGrid or in the surrounding one-cell buffer i.e. ghost atoms.
     * \endrst
     */
    template <typename E>
    int cell_idx(Eigen::MatrixBase<E> const& x) const noexcept {
      return to_1D(clamp_to_grid_idxs(x + m_cell.matrix()));
    }

  private:
    double m_r_cut = 0;

    Arr<int> m_shape = Arr<int>::Zero();
    Arr<int> m_prod_shape = Arr<int>::Zero();

    Arr<double> m_cell = Arr<double>::Zero();
    Arr<double> m_inv_cell = Arr<double>::Zero();

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

}  // namespace fly::system