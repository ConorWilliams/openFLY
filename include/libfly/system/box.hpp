#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <array>
#include <cstddef>
#include <nonstd/span.hpp>
#include <vector>

#include "libfly/utility/asserts.hpp"
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
 * are assembled into a matrix of the form:
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
 * with each column corresponding to a basis vector and all non-zero entries positive. This constrains some of the rotational degrees
 * of freedom of the parallelotope.
 *
 * The *canonical* cell/box/simulation-space is the volume of space inside the parallelotope with edges formed from these basis
 * vectors rooted at the origin. Along periodic axes atoms are allowed to have projected coordinates outside the canonical cell
 * however, along non-periodic axes atoms must be inside the canonical cell.
 *
 * LibFLY resolved periodicity using *ghost atoms*, these are the periodic images of *real atom*. To do this efficiently the canonical
 * box is split into a *grid* of smaller parallelotopes. These *grid cells* are determined by the maximum interatomic interaction
 * distance.
 *
 * \endrst
 */

namespace fly::system {}  // namespace fly::system

//    * The canonical grid cell is the parallelepiped of space spanned by the extents of the canonical cell, canon_grid_pos() maps a
//    * point into the canonical cell and then nudges it by one grid cell in the (1,1,1) direction. This leaves room for a layer of
//    * ghost atom cells around the canonical cell.