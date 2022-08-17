#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <variant>

#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/boxes/triclinic.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file box.hpp
 *
 * @brief Classes for describing the simulation space.
 *
 * \rst
 *
 * Atoms in libFLY exist within a simulation space (or box) and this space is often periodic. Periodic boxes must tessellate in
 * N-dimensional space hence, libFLY's boxes are parallelotopes. The parallelotopes is described by a series of *basis* vectors. These
 * are assembled into an upper-triangular matrix of the form:
 *
 * .. math::
 *
 *     \begin{bmatrix}
 *     a_1 & b_1 & c_1 & \cdots\\
 *     0 & b_2 & c_2 & \\
 *     0 & 0 & c_3 & \\
 *     \vdots &  & & \ddots
 *     \end{bmatrix}
 *
 * with: each column corresponding to a basis vector, all non-zero entries positive and all diagonal elements non-zero. This
 * constrains some of the rotational degrees of freedom of the parallelotope and ensures the existence of an inverse.
 *
 * The *canonical* cell/box/simulation-space is the volume of space inside the parallelotope with edges formed from the basis
 * vectors rooted at the origin. Along periodic axes atoms are allowed to have projected coordinates outside the canonical-cell
 * however, along non-periodic axes atoms must be inside the canonical cell.
 *
 * LibFLY resolved periodicity using *ghost atoms*, these are the periodic images of *real atom*. To do this efficiently the canonical
 * box is split into a *grid* of smaller parallelotopes. These *grid cells* are determined by the maximum interatomic interaction
 * distance.
 *
 * \endrst
 */

namespace fly::system {

  namespace detail {
    inline std::variant<Orthorhombic, Triclinic> build_varient(Mat const& ex, Arr<bool> const& pd) {
      //
      Mat skew = ex;

      for (int i = 0; i < spatial_dims; i++) {
        skew(i, i) = 0;
      }

      if (gnorm(skew) < 0.0001) {
        return Orthorhombic(ex.diagonal(), pd);
      } else {
        return Triclinic(ex, pd);
      }
    }
  }  // namespace detail

  /**
   * @brief Generalised simulation box.
   *
   * Essentially a ``std::variant`` of all supported crystal-systems /boxes.
   */
  class Box {
  public:
    /**
     * @brief The result of calling Box::make_grid()
     */
    using Grid = std::variant<OrthoGrid, TriGrid>;

    /**
     * @brief Construct a new Box box object.
     *
     * Automatically select the optimal underlying crystal system.
     *
     * @param basis Matrix with columns equal to the basis vectors of the box (parallelotope).
     * @param periodic True for each periodic axis.
     */
    Box(Mat const& basis, Arr<bool> const& periodic) : m_sys(detail::build_varient(basis, periodic)) {}

    /**
     * \copydoc Triclinic::basis
     */
    Mat basis() const {
      return std::visit([](auto const& m_box) -> Mat { return m_box.basis(); }, m_sys);
    }

    /**
     * \copydoc Triclinic::periodic
     */
    bool periodic(Eigen::Index i) const {
      return std::visit([i](auto const& m_box) -> bool { return m_box.periodic(i); }, m_sys);
    }

    /**
     * \copydoc Triclinic::canon_image
     */
    template <typename E>
    Vec canon_image(Eigen::MatrixBase<E> const& x) const {
      return std::visit([&x](auto const& m_box) -> Vec { return m_box.canon_image(x); }, m_sys);
    }

    /**
     * @brief Make an ``std::variant`` of grid objects.
     *
     * Grids have a common API, unwrap with ``std::visit``.
     *
     * @param r_cut The cut-off radius for atomic interactions.
     */
    std::variant<OrthoGrid, TriGrid> make_grid(double r_cut) const {
      return std::visit([r_cut](auto const& m_box) -> std::variant<OrthoGrid, TriGrid> { return m_box.make_grid(r_cut); }, m_sys);
    }

    /**
     * @brief Check if Box is currently using ``T`` as its crystal system.
     */
    template <typename T>
    bool holding() const noexcept {
      return std::holds_alternative<T>(m_sys);
    }

    /**
     * @brief Get the variant underlying this box.
     */
    auto get() const noexcept -> std::variant<Orthorhombic, Triclinic> const& { return m_sys; }

    /**
     * @brief Comparison operator, no surprises.
     */
    friend bool operator==(Box const& a, Box const& b) noexcept { return a.m_sys == b.m_sys; }

  private:
    std::variant<Orthorhombic, Triclinic> m_sys;
  };

}  // namespace fly::system

//    * The canonical grid cell is the parallelepiped of space spanned by the extents of the canonical cell, canon_grid_pos() maps a
//    * point into the canonical cell and then nudges it by one grid cell in the (1,1,1) direction. This leaves room for a layer of
//    * ghost atom cells around the canonical cell.