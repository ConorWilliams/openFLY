#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: MPL-2.0

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
   * @brief Provides details of the simulations Orthogonal supercell geometry,
   *
   * All queries of the space in which the atoms exist are provided by this class. It is assumed
   * (and must be ensured) all non-periodic atoms are within the Orthorhombic extents.
   */
  class Orthorhombic {
  public:
    /**
     * @brief Construct a new Orthorhombic object.
     *
     * @param extents Length of simulation box along each axis.
     * @param periodic True for each periodic axis.
     */
    Orthorhombic(Arr<Position::scalar_t> const &extents, Arr<bool> const &periodic)
        : m_extents{extents}, m_periodic{periodic}, m_inv_extents(1.0 / extents) {
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
    Arr<bool> const &periodic() const noexcept { return m_periodic; }

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
    template <typename T> Arr<floating> canon_image(Eigen::ArrayBase<T> const &x) const noexcept {
      // Non-periodic atoms are within the simulation box extents so x[i] * inv_extents less than 1 and x[i]
      // remains unaffected, hence no non-periodic switch/select.
      ASSERT((m_periodic || (x >= Arr<floating>::Zero() && x < m_extents)).all(), "Out of box");
      return x - m_extents * (x * m_inv_extents).floor();
    }

    /**
     * @brief Compute shortest vector connecting ``a`` to a periodic image of ``b``.
     *
     * This function is branchy and should be avoided in hot code.
     *
     * See: https://doi.org/10.1524/zpch.2013.0311
     */
    template <typename A, typename B>
    Vec<floating> min_image(Eigen::ArrayBase<A> const &a, Eigen::ArrayBase<B> const &b) const noexcept {
      Arr<floating> dr = b - a;
      return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
    }

    /**
     * @brief Comparison operator, no surprises.
     */
    friend bool operator==(Orthorhombic const &a, Orthorhombic const &b) noexcept {
      return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
    }

  private:
    Arr<floating> m_extents = Arr<floating>::Zero();
    Arr<bool> m_periodic = Arr<bool>::Zero();
    Arr<floating> m_inv_extents = Arr<floating>::Zero();
  };

}  // namespace fly::system