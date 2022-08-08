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
#include <cmath>
#include <complex>
#include <optional>

#include "libfly/system/atom.hpp"
#include "libfly/system/boxes/hypergrid.hpp"
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file triclinic.hpp
 *
 * @brief Specialised simulation box for Triclinic supercells.
 */

namespace fly::system {

  /**
   * @brief Maps position vector -> ND integer-tuples -> 1D index.
   *
   * TriGrid is constructable using the Triclinic::make_grid factory.
   *
   * Although the grid is for an triclinic cell it still uses a HyperGrid as triclinic cells are worse approximations to
   * spheres so the neighbour list construction would become more expensive.
   *
   */
  class TriGrid : public HyperGrid {
  public:
    /**
     * \copydoc OrthoGrid::gen_image
     */
    template <Sign S>
    std::optional<Position::matrix_t> gen_image(Position::matrix_t x, int ax) {
      if constexpr (S == Sign::plus) {
        // Shortest distance point to hyperplane
        auto dx = gdot(x, m_hyper.col(ax));
        ASSERT(dx > 0, "sign error");
        if (dx < HyperGrid::r_cut()) {
          return x + m_basis.col(ax);
        }
      } else {
        auto dx = gdot(m_basis.col(ax) - x, m_hyper.col(ax));
        ASSERT(dx > 0, "sign error");
        if (dx < HyperGrid::r_cut()) {
          return x - m_basis.col(ax);
        }
      }

      return std::nullopt;
    }

  private:
    friend class Triclinic;

    Mat<Position::scalar_t> m_basis = Mat<Position::scalar_t>::Zero();
    Mat<Position::scalar_t> m_hyper = Mat<Position::scalar_t>::Zero();

    // Private constructor
    TriGrid(Mat<Position::scalar_t> const& bs, Position::scalar_t r_cut) : HyperGrid(bs.colwise().sum(), r_cut), m_basis(bs) {
      // Compute hyperplane normals and hyper plane spacing.
      for (int i = 0; i < spatial_dims; i++) {
        //
        Mat<Position::scalar_t> points = m_basis;

        points.col(i).array() = 0;  // ith hyperplane through all axes EXCEPT ith

        auto n = hyperplane_normal(points);
        auto w = gdot(n, m_basis.col(i));

        VERIFY(std::abs(w) > r_cut, "Box too small for this r_cut");

        // Want normal in same direction as basis/edge.
        if (w < 0) {
          m_hyper.col(i) = -n;
        } else {
          m_hyper.col(i) = n;
        }
      }
    }
  };

  /**
   * @brief Generalized triclinic simulation box.
   */
  class Triclinic {
  public:
    /**
     * @brief Construct a new Triclinic box object.
     *
     * @param basis Matrix with columns equal to the basis vectors of the parallelotope.
     * @param pd True for each periodic axis.
     */
    Triclinic(Mat<Position::scalar_t> const& basis, Arr<bool> const& pd)
        : m_basis{basis.triangularView<Eigen::Upper>()}, m_basis_inv{m_basis.inverse()}, m_periodic{pd} {
      [[maybe_unused]] Position::scalar_t eps = 1e-5;

      VERIFY(((basis - m_basis).array().abs() < eps).all(), "Basis must be an upper triangular matrix.");

      VERIFY((m_basis.array() >= 0).all(), "Basis elements must be positive");
      VERIFY((m_basis.array().pow(2).colwise().sum().sqrt() > eps).all(), "Basis vectors too small");
      VERIFY((m_basis.diagonal().array() > eps).all(), "Diagonal elements must be non-zero");

      ASSERT(m_basis.determinant() > eps, "Should be invertible by above.");
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

      // Do the canonizing in the fractional basis.
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
     * @brief Make a TriGrid object.
     *
     * @param r_cut The cut-off radius for atomic interactions.
     */
    TriGrid make_grid(Position::scalar_t r_cut) const { return {m_basis, r_cut}; }

  private:
    Mat<Position::scalar_t> m_basis = Mat<Position::scalar_t>::Zero();
    Mat<Position::scalar_t> m_basis_inv = Mat<Position::scalar_t>::Zero();

    Arr<bool> m_periodic = Arr<bool>::Zero();
  };

}  // namespace fly::system