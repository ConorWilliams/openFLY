#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <fmt/core.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <utility>

#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs_core.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"

/**
 * \file lbfgs.hpp
 *
 * @brief Find local minima of a potential energy surface using the Limited-memory BFGS method (LBFGS).
 *
 * The LBFGS algorithm stores a portion of the history of a pathway through phase-space as it move towards a local minimum. This
 * history is used to approximate the inverse hessian which is then used to compute the approximate newton step towards the minimum.
 */

namespace fly::minimise {

  /**
   * @brief A minimiser that uses the LBFGS algorithm.
   */
  class LBFGS {
  public:
    /**
     * @brief Used to configure the minimiser.
     */
    struct Options {
      /** @brief Number of previous steps held in memory. */
      int n = 10;
      /** @brief Number of steps before exit with failure. */
      int iter_max = 2000;
      /** @brief Force convergence criterion (eV/Angstroms). */
      double f2norm = 1e-6;
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
      double max_trust = 0.50;
      /** @brief Minimum trust radius e.g initial step size (Angstroms). */
      double min_trust = 0.01;
      /** @brief Trust radius expansion rate. */
      double grow_trust = 1.5;
      /** @brief Trust radius contraction rate. */
      double shrink_trust = 0.5;
      /** @brief Print out debug info. */
      bool debug = false;
      /**
       * @brief If provided at each frame of the minimisation ``out`` will be written to this file.
       *
       * It is the users responsibility to ensure the lifetime of ``fout`` is at least as long as the lifetime of the ``LBFGS`` object
       * and ensure only a single ``LBFGS`` object writes to ``fout`` at any one time. This really only exist for debugging purposes...
       */
      io::BinaryFile *fout = nullptr;
    };

    /**
     * @brief Construct a new LBFGS object
     *
     * @param opt Options strut.
     * @param box Description of simulation space.
     */
    LBFGS(Options const &opt, system::Box const &box) : m_opt{opt}, m_core{opt.n}, m_box{box} {}

    /**
     * @brief Find the nearest local minimum of a potential.
     *
     * \rst
     *
     * .. todo::
     *    Currently there is a small but non negligible non-parallelised part of this workload in the LBFGS step function. If this is a
     *    bottleneck (require profiling in real workload) consider working on this.
     *
     * \endrst
     *
     * @param out Write final position/gradient here (forwarded to potential).
     * @param in Inputs required by potential and at-least Position.
     * @param pot Potential to minimise.
     * @param num_threads Number of openMP threads to use.
     * @return False if converged to a minima.
     */
    auto minimise(system::SoA<Position &, PotentialGradient &> out,
                  system::SoA<Position const &, TypeID const &, Frozen const &> in,
                  potential::Generic &pot,
                  int num_threads = 1) -> bool;

  private:
    Options m_opt;
    StepLBFGS m_core;

    system::Box m_box;
    double m_r_cut = -1;

    std::optional<neigh::List> m_nl;
  };

}  // namespace fly::minimise
