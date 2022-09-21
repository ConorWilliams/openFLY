#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <memory>
#include <utility>

#include "libfly/neigh/list.hpp"
#include "libfly/potential/property.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"

/**
 * \file lbfgs.hpp
 *
 * @brief ...
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
      std::size_t n = 10;
      /** @brief Number of steps before exit with failure. */
      std::size_t iter_max = 2000;
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
      /** @brief If minimising a dimer and convex_max steps with +Ve curvature exit. */
      double convex_max = 3;
      /** @brief Print out debug info and dumps minimisation trace to "lbfgs_debug.xyz" */
      bool debug = false;
    };

    /**
     * @brief Construct a new LBFGS object
     *
     * @param opt
     */
    explicit LBFGS(Options const &opt) : m_opt{opt}, m_core{opt.n} {}

    /**
     * @brief Move the atoms in the SimCell to a local minimum of the potential.
     *
     * @param in
     * @param threads
     * @return true
     * @return false
     */
    bool minimise(system::SoA<TypeID const &, Frozen const &> in, potential::Base &, int threads = 1);

  private:
    Options m_opt;
    StepLBFGS m_core;
    std::optional<neigh::List> m_nl;
  };

}  // namespace fly::minimise