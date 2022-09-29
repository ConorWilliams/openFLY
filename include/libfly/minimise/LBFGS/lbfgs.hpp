#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

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
      int n = 10;
      /** @brief Number of steps before exit with failure. */
      int iter_max = 2000;
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
     * @brief Return codes for ``minimise()``.
     *
     */
    enum Exit : int {
      success = 0,         ///< Found minimum.
      positive_curv = -1,  ///< In positive curvature region for too long (must be using a Dimer as the potential).
      iter_max = -2,       ///< Exceeded the maximum number of iterations.
    };

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
     * @return Exit Status code detailing why minimisation stopped.
     */

    template <typename... P1, typename... P2>
    auto minimise(system::SoA<P1...> &out, system::SoA<P2...> const &in, potential::Generic &pot, int num_threads = 1) -> Exit {
      // Check inputs

      verify(in.size() == out.size(), "LBFGS minimizer inputs size mismatch, in={} out={}", in.size(), out.size());

      m_core.clear();  // Clear history from previous runs.

      double skin = std::max(std::pow(m_opt.skin_frac, 1. / 3.) - 1, 0.0) * pot.r_cut();

      if (m_opt.debug) {
        fmt::print("LBFGS: threads = {}\n", num_threads);
        fmt::print("LBFGS: Skin = {}\n", skin);
      }

      // Avoid floating point comparison warnings.
      constexpr std::equal_to<> eq;

      if (double r_cut = pot.r_cut() + skin; !eq(std::exchange(m_r_cut, r_cut), r_cut) || !m_nl) {
        m_nl = neigh::List(m_box, r_cut);

        if (m_opt.debug) {
          fmt::print("LBFGS: Reallocating neigh list\n");
        }
      }

      out[r_] = in[r_];  // Initialise out = in;

      m_nl->rebuild(out, num_threads);

      pot.gradient(out, in, *m_nl, num_threads);

      double trust = m_opt.min_trust;

      double acc = 0;
      double convex_count = 0;

      for (int i = 0; i < m_opt.iter_max; ++i) {
        //
        double mag_g = gnorm_sq(out[g_]);

        if (m_opt.debug) {
          fmt::print("LBFGS: i={:<4} trust={:f} acc={:f} norm(g)={:e} rebuild={}\n", i, trust, acc, std::sqrt(mag_g), eq(acc, 0.0));
        }

        if (m_opt.fout) {
          m_opt.fout->commit([&] { (m_opt.fout->write(P1{}, out), ...); });
        }

        if (mag_g < m_opt.f2norm * m_opt.f2norm) {
          return success;
        } else if (convex_count >= m_opt.convex_max) {
          if (m_opt.debug) {
            fmt::print("LBFGS: Early exit - curvature\n");
          }
          return positive_curv;
        }

        auto &Hg = m_core.newton_step<Position, PotentialGradient>(out, out);

        ASSERT(gdot(out[g_], Hg[del_]) > 0, "Ascent direction in lbfgs", 0);

        // Limit step size.
        Hg[del_] *= std::min(1.0, trust / gnorm(Hg[del_]));

        // Add distance of most displaced atom
        acc += std::sqrt((Hg[del_] * Hg[del_]).colwise().sum().maxCoeff());

        // Update positions in real space
        out[r_] -= Hg[del_];

        if (acc > 0.5 * skin) {
          m_nl->rebuild(out, num_threads);
          acc = 0;
        } else {
          m_nl->update(Hg);
        }

        pot.gradient(out, in, *m_nl, num_threads);

        pot.visit([&](auto const &underlying) {
          if constexpr (std::is_same_v<remove_cref_t<decltype(underlying)>, potential::Dimer>) {
            if (underlying.curv() > 0) {
              convex_count += 1;
            } else {
              convex_count = 0;
            }
          }
        });

        double proj = gdot(out[g_], Hg[del_]);

        if (proj < -m_opt.proj_tol) {
          trust = std::max(m_opt.min_trust, m_opt.shrink_trust * trust);
        } else if (proj > m_opt.proj_tol) {
          trust = std::min(m_opt.max_trust, m_opt.grow_trust * trust);
        }
      }

      return iter_max;
    }

  private:
    Options m_opt;
    StepLBFGS m_core;

    system::Box m_box;
    double m_r_cut = -1;

    std::optional<neigh::List> m_nl;
  };

}  // namespace fly::minimise