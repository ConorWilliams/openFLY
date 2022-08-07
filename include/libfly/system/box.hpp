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
 * @brief Static polymorphism for box.
 */

namespace fly::system {}  // namespace fly::system

//    * The canonical grid cell is the parallelepiped of space spanned by the extents of the canonical cell, canon_grid_pos() maps a
//    * point into the canonical cell and then nudges it by one grid cell in the (1,1,1) direction. This leaves room for a layer of
//    * ghost atom cells around the canonical cell.