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
#include "libfly/saddle/rotor.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file dimer.hpp
 *
 * @brief Saddle-point finding with the dimer method.
 */

namespace fly::saddle {

  /**
   * @brief Saddle point locator.
   *
   * A dimer is two images of a system. These images have a centre point and a unit axis. The images are (conceptually) located
   * ``delta_r`` in the plus/minus directions along the axis.
   *
   * This implementation uses the superlinear dimer method to find saddle-points. This alternates optimizing the orientation of the
   * dimer along the minimum eigen-mode then translating the dimer along an effective potential to find the SP.
   */
  class Dimer {
  public:
    /**
     * @brief Configure the dimers internal translation minimisation pass.
     */
    struct Options {
      /** @brief L-BFGS core translation-history size. */
      int n = 10;
      /** @brief Maximum number of steps before failing. */
      int max_steps = 1000;
      /** @brief Every ``hist_check_freq`` number of steps will check not converging to known SP. */
      int hist_check_freq = 5;
      /** @brief Abort if the cosine of the angle between the dimer and a known SP is greater than this. */
      double cos_theta_tol = std::cos(10. / 360. * 2 * M_PI);
      /** @brief Force convergence criterion (eV/Angstroms). */
      double f2norm = 1e-5;
      /**
       * @brief Used to determine the skin size.
       *
       * The skin thickness is determined such that the expected number of neighbours is skin_frac *
       * (num_neighbours if no skin used). Hence skin_frac must be larger than one. Neighbour lists
       * are built less often when skin_frac is higher but there will be more more non-neighbours in
       * neighbour lists.
       */
      double skin_frac = 1.1;
      /** @brief Trust tolerance, set larger to reduce trust radius change. */
      double proj_tol = 0;
      /** @brief Maximum trust radius e.g max steps size (Angstroms). */
      double max_trust = 0.5;
      /** @brief Minimum trust radius e.g initial step size (Angstroms). */
      double min_trust = 0.05;
      /** @brief Trust radius expansion rate. */
      double grow_trust = 1.5;
      /** @brief Trust radius contraction rate. */
      double shrink_trust = 0.5;
      /** @brief If  ``convex_max`` steps with +Ve curvature then exit early. */
      double convex_max = 3;
      /** @brief Print out debug info. */
      bool debug = false;
      /**
       * @brief If provided at each frame of the translation ``out`` will be written to this file.
       *
       * It is the users responsibility to ensure the lifetime of ``fout`` is at least as long as the lifetime of the ``LBFGS`` object
       * and ensure only a single ``LBFGS`` object writes to ``fout`` at any one time. This really only exist for debugging purposes...
       */
      io::BinaryFile *fout = nullptr;
    };

    /**
     * @brief Construct a new Dimer object.
     *
     * @param opt Options for translations pass.
     * @param r_opt Options for rotation pass.
     * @param box Description of simulation space.
     */
    Dimer(Options const &opt, Rotor::Options const &r_opt, system::Box const &box)
        : m_opt(opt), m_core(opt.n), m_rotor(r_opt), m_box(box) {}

    /**
     * @brief Copy construct a new Dimer object.
     */
    Dimer(Dimer const &) = default;

    /**
     * @brief Move construct a new Dimer object.
     */
    Dimer(Dimer &&) = default;

    /**
     * @brief Return codes for `step()``.
     */
    enum Exit : int {
      success = 0,     ///< Found saddle-point.
      convex = -1,     ///< Failed: stuck in +curvature region.
      iter_max = -2,   ///< Exceeded the maximum number of iterations but in -curvature region.
      collision = -3,  ///< Approached a known SP during search.
    };

    /**
     * @brief Advance the dimer towards a saddle-point.
     *
     * Performs at most ``max_step`` translation steps. The translation uses the LBFGS algorithm with a trust-radius limited step size
     * and an early-exit condition if the curvature is positive for too long.
     *
     * @param out Final position and axis orientation written here.
     * @param in Initial position, axis and per-particle data forwarded to potential.
     * @param pot Potential energy function.
     * @param in_min initial minima that dimer is climbing FROM.
     * @param hist_sp Previously discovered saddle points.
     * @param num_threads Number of openMP threads to use.
     * @return Exit-code.
     */
    auto step(system::SoA<Position &, Axis &> out,
              system::SoA<Position const &, Axis const &, TypeID const &, Frozen const &> in,
              system::SoA<Position const &> in_min,
              potential::Generic &pot,
              std::vector<system::SoA<Position>> const &hist_sp,
              int num_threads = 1) -> Exit;

  private:
    Options m_opt;

    minimise::StepLBFGS m_core;

    system::SoA<PotentialGradient> m_eff_grad;

    Rotor m_rotor;

    system::Box m_box;
    double m_r_cut = -1;

    std::optional<neigh::List> m_nl;
  };

}  // namespace fly::saddle
