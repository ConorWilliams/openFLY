#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <type_traits>

#include "libfly/utility/core.hpp"

/**
 * \file hessian.hpp
 *
 * @brief Represent and manipulate Hessians.
 */

namespace fly::system {

  /**
   * @brief A class to represent the blocked hessian of a system.
   *
   * The hessian is a symmetric nm by nm matrix with m the number of active atoms and n the number of spatial
   * dimensions. Each n by n sub-matrix is a "block" of the hessian.
   */
  class Hessian {
  public:
    /**
     * @brief The vector types used to store the eigen values of the hessian.
     */
    using Vector = Eigen::Vector<double, Eigen::Dynamic>;

    /**
     * @brief The Matrix types used to store the eigen vectors of the hessian.
     */
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    /**
     * @brief Allocates and zeros Hessian large enough for ``num_active`` atoms.
     *
     * @param num_active Number of active atoms that hessian will contain.
     */
    auto zero_for(Eigen::Index num_active) -> void {
      verify(num_active >= 0, "{} is not a valid number of active atoms", num_active);
      m_hess = Matrix::Zero(spatial_dims * num_active, spatial_dims * num_active);
    }

    /**
     * @brief Get a view of the block corresponding to the ``ij`` atom pair.
     *
     * @param i The index (in the hessian matrix) of the first atom in the pair.
     * @param j The index (in the hessian matrix) of the second atom in the pair.
     */
    auto operator()(Eigen::Index i, Eigen::Index j) {
      ASSERT(i >= 0 && i < m_hess.rows(), "{} is not a valid index with {} rows", i, m_hess.rows());
      ASSERT(j >= 0 && j < m_hess.cols(), "{} is not a valid index with {} cols", i, m_hess.cols());
      return m_hess.block<spatial_dims, spatial_dims>(spatial_dims * i, spatial_dims * j);
    }

    /**
     * @brief Compute and return a reference to an ordered vector of Eigen Values.
     *
     *
     * See: https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html only reads the lower
     * diagonal portion of the matrix.
     *
     * The eigenvalues are repeated according to their algebraic multiplicity, so there are as many
     * eigenvalues as rows in the matrix. The eigenvalues are sorted in increasing order.
     */
    Vector const& eigenvalues() {
      //
      static_assert(std::is_same_v<decltype(m_solver.eigenvalues()), Vector const&>, "Eigen allocation");

      m_solver.compute(m_hess, Eigen::EigenvaluesOnly);

      return m_solver.eigenvalues();
    }

    /**
     * @brief A struct that contains references to the eigen values/vectors.
     */
    struct eigen_t {
      Vector const& values;   ///< The n (maybe repeated) eigen values of the (n x n) matrix.
      Matrix const& vectors;  ///< The ``k``th column stores the ``k`` eigenvector corresponding to the
                              ///< ``k``th eigenvalue. The eigenvectors are normalized to have (Euclidean)
                              ///< norm equal to one.
    };

    /**
     * @brief Computes the eigenvalues and eigenvectors.
     *
     * See: https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html only reads the lower
     * diagonal portion of the matrix.
     */
    eigen_t eigen() {
      //
      static_assert(std::is_same_v<decltype(m_solver.eigenvalues()), Vector const&>, "Eigen allocation");
      static_assert(std::is_same_v<decltype(m_solver.eigenvectors()), Matrix const&>, "Eigen allocation");

      m_solver.compute(m_hess);

      return {m_solver.eigenvalues(), m_solver.eigenvectors()};
    }
    /**
     * @brief Fetch a reference to the underlying matrix.
     */
    Matrix const& get() const noexcept { return m_hess; }

  private:
    Matrix m_hess;
    Eigen::SelfAdjointEigenSolver<Matrix> m_solver;
  };

}  // namespace fly::system
