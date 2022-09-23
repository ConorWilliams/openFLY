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
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file dimer.hpp
 *
 * @brief A second-order potential.
 */

namespace fly::potential {

  /**
   * @brief Forward declare.
   */
  class Generic;

  /**
   * @brief Dimer is a potential adaptor.
   *
   * Wraps a potential and inverts the component of the gradient parallel to the minimum mode.
   */
  class Dimer {
  public:
    /**
     * @brief Configure the dimers internal minimisation pass.
     */
    struct Options {
      /** @brief L-BFGS core rotation-history size. */
      int n = 6;
      /** @brief Maximum number of iterations during rotation minimization. */
      int iter_max_rot = 20;
      /** @brief Half dimer length. */
      double delta_r = 0.001;
      /** @brief (Rad) rotation convergence criterion. */
      double theta_tol = 2 * M_PI / 360.;
      /** @brief If true when convex we return only the component parallel to the min mode. */
      bool relax_in_convex = true;
      /** @brief Print out debug info. */
      bool debug = false;
    };

    /**
     * @brief Construct a new Dimer object.
     *
     * @param opt Options for minimisation pass.
     * @param to_wrap Potential to wrap.
     */
    Dimer(Options const& opt, potential::Generic const& to_wrap);

    /**
     * @brief Copy construct a new Dimer object.
     */
    Dimer(Dimer const&);

    /**
     * @brief Move construct a new Dimer object.
     */
    Dimer(Dimer&&) noexcept;

    /**
     * @brief Get this potentials cut-off radius that the neighbour lists should be configured for.
     *
     * This is the wrapped potentials cut-off plus two times delta_r, this means we can displace every atom by delta_r during the
     * minimisation pass and the not need to rebuild the neighbour lists.
     */
    double r_cut() const noexcept;

    /**
     * @brief Compute the effective potential's gradient.
     *
     * This is achieved by inverting the component of the gradient parallel to the minimum eigen-mode of the wrapped potential. The
     *
     * Assumes the neighbour list are ready, force on frozen atoms will be zero.
     *
     * @param in Per-atom input properties
     * @param out Write effective gradient and final axis orientation here.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     * @return The approximate curvature of the wrapped potential along the output Axis.
     */
    auto eff_gradient(system::SoA<PotentialGradient&, Axis&> out,
                      system::SoA<TypeID const&, Frozen const&, Axis const&> in,
                      neigh::List& nl,
                      int threads = 1) -> double;

  private:
    Options m_opt;
    minimise::StepLBFGS m_core;
    std::unique_ptr<potential::Generic> m_wrapped;

    system::SoA<Delta> m_delta;       // Store displacement for updates
    system::SoA<Delta> m_delta_prev;  // Store previous displacement for updates

    system::SoA<Axis> m_axisp;  // Temporary axis

    system::SoA<PotentialGradient> m_g0;   // Central grad
    system::SoA<PotentialGradient> m_g1;   // Temp end grad
    system::SoA<PotentialGradient> m_g1p;  // Temp primed end grad

    system::SoA<Delta> m_delta_g;  // Gradient difference
  };

}  // namespace fly::potential