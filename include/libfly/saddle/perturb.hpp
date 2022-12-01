#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/potential/generic.hpp"
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
   * @brief Provide a random perturbation to every atom within r_cut of a point.
   *
   * Uses the minimum image distance to determine distance from centre. The perturbation is Gaussian in each coordinate axis and has
   * standard deviation ``stddev``. An envelope function will linearly decrease the size of each atoms perturbation based off its
   * distance to the centre. The axis will be normalised as appropriate. Frozen atoms will not be perturbed.
   *
   * \rst
   * .. tip::
   *    This function uses an expensive method for computing the minimum image convention. If calling this many times on the same
   *    initial condition it may be more efficient to compute and re-use a neighbour list.
   *
   * \endrst
   *
   * @param out Output the perturbed state and randomised axis here.
   * @param urbg Random number generator.
   * @param box Simulation space.
   * @param centre Origin of perturbation.
   * @param cell Input parameters.
   * @param r_cut All atoms within ``r_cut`` of ``centre`` will be perturbed.
   * @param stddev Standard deviation of Gaussian perturbations.
   */
  void perturb(system::SoA<Position&, Axis&> out,
               Xoshiro& urbg,
               system::Box const& box,
               Position::matrix_t const& centre,
               system::SoA<Position const&, Frozen const&> cell,
               double r_cut,
               double stddev);

}  // namespace fly::saddle