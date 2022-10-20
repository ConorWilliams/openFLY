#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.centroid>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/SoA.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file mechanisms.hpp
 *
 * @brief Representations of mechanisms.
 */

namespace fly::system {

  /**
   * @brief A representation of a local mechanism containing the the displacements of the atoms in a LE.
   *
   * Frozen atoms must have a displacement equal to ``Vec::Zero()``.
   */
  class LocalMech {
  public:
    double barrier;       ///< Energy barrier (eSp - e0).
    double delta;         ///< Energy change (eF - e0).
    double kinetic_pre;   ///< Arrhenius pre factor.
    double capture_frac;  ///< norm(delta_fwd in local) / norm(delta_fwd in global)

    system::VoS<Delta> delta_sp;   ///< Displacement vectors from initial to saddle-point.
    system::VoS<Delta> delta_fwd;  ///< Displacement vectors from initial to final.
  };

}  // namespace fly::system