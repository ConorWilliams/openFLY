#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cmath>
#include <cstddef>
#include <memory>

#include "libfly/minimise/LBFGS/lbfgs_core.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file rotor.hpp
 *
 * @brief Class to orient a dimer.
 */

namespace fly::saddle {

  /**
   * @brief Orients a dimer to align it with the minimum eigen-mode.
   */
  class Rotor {
  public:
    /**
     * @brief Configure the dimers rotation minimisation pass.
     */
    struct Options {
      /** @brief L-BFGS core rotation-history size. */
      int n = 6;
      /** @brief Maximum number of iterations during rotation minimization. */
      int iter_max_rot = 20;
      /** @brief Half dimer length. */
      double delta_r = 0.001;
      /** @brief (Rad) rotation convergence criterion. */
      double theta_tol = 1 * M_PI / 360.;
      /** @brief If true when convex we return only the component parallel to the min mode. */
      bool relax_in_convex = true;
      /** @brief Print out debug info. */
      bool debug = false;
    };

    /**
     * @brief Construct a new Rotor object.
     *
     * @param opt Options for rotor minimisation pass.
     */
    explicit Rotor(Options const &opt) : m_opt(opt), m_core(opt.n) {}

    /**
     * @brief Copy construct a new Rotor object.
     */
    Rotor(Rotor const &) = default;

    /**
     * @brief Move construct a new Rotor object.
     */
    Rotor(Rotor &&) = default;

    /**
     * @brief Compute the effective potential's gradient.
     *
     * This is achieved by inverting the component of the gradient parallel to the minimum eigen-mode of the potential. The dimer is
     * rotated using the LBFGS algorithm to make it align with the minimum eigen-mode.
     *
     * Assumes the neighbour list is ready, force on frozen atoms will be zero.
     *
     * \rst
     *
     * .. note::
     *    This function does actually modify neighbour but returns it to an identical state after it is done. It guarantees not to call
     *    ``neigh::List::rebuild()``.
     *
     * \endrst
     *
     * @param in Per-atom input properties
     * @param out Write effective gradient here.
     * @param inout The axis input and output.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called) configured with a cut-off
     * at least ``r_cut()``.
     * @param pot The potential energy function.
     * @param threads Number of openMP threads to use.
     * @return The approximate curvature of the wrapped potential along the output Axis.
     */
    auto eff_gradient(system::SoA<PotentialGradient &> out,
                      system::SoA<Axis &> inout,
                      system::SoA<TypeID const &, Frozen const &> in,
                      potential::Generic &pot,
                      neigh::List &nl,
                      int threads = 1) -> double;

    /**
     * @brief Compute the cut-off radius.
     *
     * This is the cut-off that ``nl`` in call to ``eff_grad()`` with potential ``pot`` must be configured with.
     *
     * @param pot Potential that ``eff_grad()`` will be called with.
     * @return auto A slightly larger cut-off that enables ``eff_grad`` to displace atoms a little without rebuilding the nl.
     */
    auto r_cut(potential::Generic const &pot) const noexcept { return pot.r_cut() + m_opt.delta_r * 2; }

  private:
    Options m_opt;
    minimise::StepLBFGS m_core;

    system::SoA<Delta> m_delta;       // Store displacement for updates
    system::SoA<Delta> m_delta_prev;  // Store previous displacement for updates

    system::SoA<Axis> m_axisp;  // Temporary axis

    system::SoA<PotentialGradient> m_g0;   // Central grad
    system::SoA<PotentialGradient> m_g1;   // Temp end grad
    system::SoA<PotentialGradient> m_g1p;  // Temp primed end grad

    system::SoA<Delta> m_delta_g;  // Gradient difference
  };

}  // namespace fly::saddle