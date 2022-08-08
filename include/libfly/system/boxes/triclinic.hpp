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
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file triclinic.hpp
 *
 * @brief Specialised simulation box for Triclinic supercells.
 */

namespace fly::system {

  /**
   * @brief Generalized triclinic simulation box.
   */
  class Triclinic {
  public:
    class Grid;

    /**
     * @brief Construct a new Triclinic box object.
     *
     * @param basis Matrix with columns equal to the basis vectors of the parallelotope.
     * @param pd True for each periodic axis.
     */
    Triclinic(Mat<Position::scalar_t> const& basis, Arr<bool> const& pd)
        : m_basis{basis.triangularView<Eigen::Upper>()}, m_basis_inv{m_basis.inverse()}, m_periodic{pd} {
      // Verify basis
      Position::scalar_t eps = 1e-5;

      VERIFY(((basis - m_basis).array().abs() < eps).all(), "Basis must be an upper triangular matrix.");

      //
      VERIFY((m_basis.array() >= 0).all(), "Basis elements must be positive");
      VERIFY((m_basis.array().pow(2).colwise().sum().sqrt() > eps).all(), "Basis vectors too small");
      VERIFY((m_basis.diagonal().array() > eps).all(), "Diagonal elements must be non-zero");

      ASSERT(m_basis.determinant() > eps, "Should be invertible by above.");

      // Compute hyperplane normals
      for (int i = 0; i < spatial_dims; i++) {
        //
        Mat<Position::scalar_t> points = m_basis;
        points.col(i).array() = 0;

        // Want normal in same direction as edge.
        if (auto n = hyperplane_normal(points); gdot(n, m_basis.col(i)) < 0) {
          //   m_hyper_normals.col(i) = -n;
        } else {
          //   m_hyper_normals.col(i) = n;
        }
      }
    }

    /**
     * @brief Fetch the basis vectors of this cell.
     *
     * @return  Mat<Position::scalar_t> A matrix with each column corresponding to a basis vector.
     */
    Mat<Position::scalar_t> basis() const noexcept { return m_basis; }

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
     * \endrst
     */
    template <typename E>
    Position::matrix_t canon_image(Eigen::MatrixBase<E> const& x) const {
      // Convert to fractional coordinates
      Position::matrix_t f = m_basis_inv * x;

      ASSERT((m_periodic || (f.array() >= 0 && f.array() < 1)).all(), "Out of box A");

      // Do the canonizing in fractional basis.
      Position::matrix_t f_canon = f.array() - f.array().floor();

      return m_basis * f_canon;  // Transform back to real coordinates
    }

    /**
     * @brief Comparison operator, no surprises.
     */
    friend bool operator==(Triclinic const& a, Triclinic const& b) noexcept {
      return a.m_basis == b.m_basis && (a.m_periodic == b.m_periodic).all();
    }

    /**
     * @brief Make a Grid object.
     *
     * @param r_cut The cut-off radius for atomic interactions.
     */
    Grid make_grid(Position::scalar_t r_cut) const;  // Implementation after Grid at EoF

  private:
    Mat<Position::scalar_t> m_basis = Mat<Position::scalar_t>::Zero();
    Mat<Position::scalar_t> m_basis_inv = Mat<Position::scalar_t>::Zero();

    Arr<bool> m_periodic = Arr<bool>::Zero();
  };

  //   /**
  //    * @brief Maps position vector -> ND integer-tuples -> 1D index.
  //    *
  //    * Grid is constructable using the Triclinic::make_grid factory.
  //    */
  //   class Triclinic::Grid {
  //   public:
  //     /**
  //      * @brief Get the total number of cells in the Grid along each axis.
  //      */
  //     Arr<int> shape() const noexcept { return m_shape; }

  //   private:
  //     friend class Triclinic;

  //     Arr<Position::scalar_t> m_extents = Arr<Position::scalar_t>::Zero();

  //     Position::scalar_t m_r_cut = 0;

  //     Mat<Position::scalar_t> m_hyper_normals = Mat<Position::scalar_t>::Zero();

  //     Arr<int> m_shape = Arr<int>::Zero();
  //     Arr<int> m_prod_shape = Arr<int>::Zero();

  //     // Arr<Position::scalar_t> m_cell = Arr<Position::scalar_t>::Zero();
  //     // Arr<Position::scalar_t> m_inv_cell = Arr<Position::scalar_t>::Zero();

  //     // Private constructor
  //     Grid(Triclinic const& box, Position::scalar_t r_cut) : m_extents(box.m_basis), m_r_cut(r_cut) {
  //       // Sanity checks
  //       VERIFY(r_cut > 0, "r_cut is negative");
  //       VERIFY((extents > r_cut).all(), "r_cut is too big");
  //       //
  //       m_shape = 2 + (extents / r_cut).cast<int>();
  //       m_cell = extents / (extents / r_cut).floor();
  //       m_inv_cell = 1.0 / m_cell;

  //       // Cumulative product of m_shape

  //       m_prod_shape = Arr<int>::Ones();

  //       for (int i = 1; i < spatial_dims; i++) {
  //         m_prod_shape[i] = m_prod_shape[i - 1] * m_shape[i - 1];
  //       }
  //     }

  //     /**
  //      * @brief Cast a position to its grid indexes and clamp on interval [0, m_shape -1].
  //      */
  //     template <typename E>
  //     Arr<int> clamp_to_grid_idxs(Eigen::MatrixBase<E> const& x) const noexcept {
  //       return (x.array() * m_inv_cell).template cast<int>().cwiseMax(0).cwiseMin(m_shape - 1);
  //     }

  //     /**
  //      * @brief Get the 1D index from the nD index.
  //      */
  //     int to_1D(Arr<int> const& indexes) const noexcept { return (indexes * m_prod_shape).sum(); }
  //   };

}  // namespace fly::system