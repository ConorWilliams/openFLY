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
     * @brief Construct a new Ortho Sim Box object.
     *
     * @param extents Length of simulation box along each axis.
     * @param periodic True for each periodic axis.
     */
    Orthorhombic(Vec3<floating> const &extents, Vec3<bool> const &periodic)
        : m_extents{extents}, m_periodic{periodic}, m_inv_extents(1.0 / extents) {
      VERIFY((m_extents > 0).all(), "Orthorhombic extents are negative");
    }

    /**
     * @brief Extents getter (const).
     */
    Vec3<floating> const &extents() const noexcept { return m_extents; }

    /**
     * @brief Periodicity getter (const).
     */
    Vec3<bool> const &periodic() const noexcept { return m_periodic; }

    /**
     * @brief Maps atom into canonical cell, 0 <= r_i < extent_i for all i which are periodic.
     * Non-periodic atoms are within the simbox extents so x[i] * inv_extents less than 1 and x[i]
     * remains unaffected, hence no non-periodic switch/select.
     */
    template <typename T> Vec3<floating> canon_image(Eigen::ArrayBase<T> const &x) const noexcept {
      ASSERT((m_periodic || (x >= Vec3<floating>::Zero() && x < m_extents)).all(), "Out of box");
      return x - m_extents * (x * m_inv_extents).floor();
    }

    /**
     * @brief Compute the shortest Vector connecting a to a periodic image of b. This function is
     * branchy and should be avoided in hot code.
     */
    template <typename A, typename B>
    Vec3<floating> min_image(Eigen::ArrayBase<A> const &a, Eigen::ArrayBase<B> const &b) const noexcept {
      Vec3<floating> dr = b - a;
      return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
    }

    /**
     * @brief Comparison operator, no surprises.
     */
    friend bool operator==(Orthorhombic const &a, Orthorhombic const &b) noexcept {
      return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
    }

  private:
    Vec3<floating> m_extents = Vec3<floating>::Zero();
    Vec3<bool> m_periodic = Vec3<bool>::Zero();
    Vec3<floating> m_inv_extents = Vec3<floating>::Zero();
  };

}  // namespace fly::system