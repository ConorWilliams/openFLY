#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file perturb.hpp
 *
 * @brief Utility for perturbing a supercell.
 */

namespace fly::saddle {

  /**
   * @brief Provide a random pertubation to every atom within rcut of centre.
   *
   * Uses the minimum image distance to determine distance from centre. The pertubation is gaussian
   * in each coordinate axis and has standard deviation stddev. An envelope function will linearly
   * decrease the size of each atoms pertubation based off its distance to the centre.
   *
   * @param out
   * @param urbg
   * @param box
   * @param centre
   * @param cell
   * @param rcut
   * @param stddev
   */
  void perturb(system::SoA<Position&, Axis&> out,
               Xoshiro& urbg,
               system::Box const& box,
               Position::matrix_t const& centre,
               system::SoA<Position const&, Frozen const&> cell,
               double rcut,
               double stddev);

}  // namespace fly::saddle