#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

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
 * \file find_sp.hpp
 *
 * @brief Utility for perturbing a supercell.
 */

namespace fly::saddle {

  /**
   * @brief
   *
   * @tparam F
   * @param urbg
   * @param in
   * @param minimiser
   * @param pot
   */
  template <typename F>
  void find_connected_min(Xoshiro& urbg, system::SoA<Position const&, Frozen const&> in, F& minimiser, potential::Generic& pot) {
    //
    // Store working
    system::SoA<Position, Frozen const&, Axis, PotentialGradient> probe(in.size());

    dcell[r_] = cell[r_];
    dcell[fzn_] = cell[fzn_];
    dcell[id_] = cell[id_];
    dcell[ax_] = 1;
    dcell[ax_] /= gnorm(dcell[ax_]);

    //
    potential::Generic dimer{potential::Dimer{
        {},
        pot,
    }};
  }

}  // namespace fly::saddle