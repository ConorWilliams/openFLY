// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

#include "libfly/saddle/rotor.hpp"

#include <cmath>
#include <cstddef>
#include <memory>

#include "libfly/minimise/LBFGS/lbfgs_core.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::saddle {

  template <typename T>
  static void trans_proj(system::SoA<T &> x) {
    //
    Vec proj = Vec::Zero();

    for (Eigen::Index i = 0; i < x.size(); ++i) {
      proj += x(T{}, i);
    }

    proj *= 1.0 / static_cast<typename T::scalar_t>(x.size());

    for (Eigen::Index i = 0; i < x.size(); ++i) {
      x(T{}, i) -= proj;
    }
  }

  auto Rotor::eff_gradient(system::SoA<PotentialGradient &> out,
                           system::SoA<Axis &> inout,
                           system::SoA<TypeID const &, Frozen const &> in,
                           potential::Generic &pot,
                           neigh::List &nl,
                           int count_frozen,
                           int num_threads) -> double {
    //
    verify(in.size() == out.size(), "Effective gradient size mismatch in={} out={}", in.size(), out.size());
    verify(
        inout.size() == out.size(), "Axis gradient size mismatch inout={} out={}", inout.size(), out.size());

    m_delta.destructive_resize(in.size());       // Store displacement for updates
    m_delta_prev.destructive_resize(in.size());  // Store previous displacement for updates

    m_axisp.destructive_resize(in.size());  // Temporary axis

    m_g0.destructive_resize(in.size());   // Central grad
    m_g1.destructive_resize(in.size());   // Temp end grad
    m_g1p.destructive_resize(in.size());  // Temp primed end grad

    m_delta_g.destructive_resize(in.size());  // Gradient difference

    // Reset LBFGS history
    m_core.clear();

    pot.gradient(m_g0, in, nl, num_threads);  // Gradient at centre (g0)

    m_delta[del_] = -m_opt.delta_r * inout[ax_];

    nl.update(m_delta);

    using std::swap;  // ADL

    swap(m_delta, m_delta_prev);

    pot.gradient(m_g1, in, nl, num_threads);  // Gradient at end (g1)

    double m_curv = [&] {
      for (int i = 0;; i++) {
        m_delta_g[del_] = m_g1[g_] - m_g0[g_];
        m_delta_g[del_] -= gdot(m_delta_g[del_], inout[ax_]) * inout[ax_];  // Torque

        if (count_frozen == 0) {
          trans_proj<Delta>(m_delta_g);
        }

        // Use lbfgs to find rotation plane
        auto &theta = m_core.newton_step<Axis, Delta>(inout, m_delta_g);

        theta[del_] -= gdot(theta[del_], inout[ax_]) * inout[ax_];  // Ortho
        theta[del_] *= 1 / gnorm(theta[del_]);                      //      normalization

        double b_1 = gdot(m_g1[g_] - m_g0[g_], theta[del_]) / m_opt.delta_r;
        double c_x0 = gdot(m_g1[g_] - m_g0[g_], inout[ax_]) / m_opt.delta_r;
        double theta_1 = -0.5 * std::atan(b_1 / std::abs(c_x0));  // Trial rotation angle

        if (std::abs(theta_1) < m_opt.theta_tol || (i >= m_opt.iter_max_rot && c_x0 < 0)
            || i > 10 * m_opt.iter_max_rot) {
          return c_x0;
        } else {
          // Trial rotation
          m_axisp[ax_] = inout[ax_] * std::cos(theta_1) + theta[del_] * std::sin(theta_1);

          m_delta[del_] = -m_opt.delta_r * m_axisp[ax_];  // Temporarily store next m_delta_prev into m_delta
          swap(m_delta, m_delta_prev);                    // Now put it in the correct place
          m_delta[del_] = m_delta_prev[del_] - m_delta[del_];
          nl.update(m_delta);

          pot.gradient(m_g1p, in, nl, num_threads);  // Gradient at primed end (g1p)

          double c_x1 = gdot(m_g1p[g_] - m_g0[g_], m_axisp[ax_]) / m_opt.delta_r;
          double a_1 = (c_x0 - c_x1 + b_1 * sin(2 * theta_1)) / (1 - std::cos(2 * theta_1));
          double theta_min = 0.5 * std::atan(b_1 / a_1);  // Optimal rotation

          // Flip if extrema is maxima
          if (a_1 * std::cos(2 * theta_min) - a_1 + b_1 * std::sin(2 * theta_min) > 0) {
            theta_min += M_PI / 2;
          }

          inout[ax_] = inout[ax_] * std::cos(theta_min) + theta[del_] * std::sin(theta_min);

          // Interpolate force at new rotation
          m_g1[g_] = (std::sin(theta_1 - theta_min) / std::sin(theta_1)) * m_g1[g_]
                     + (std::sin(theta_min) / std::sin(theta_1)) * m_g1p[g_]
                     + (1 - std::cos(theta_min) - std::sin(theta_min) * std::tan(0.5 * theta_1)) * m_g0[g_];

          if (m_opt.debug) {
            fmt::print("\tDimer: i={:<4} theta={:f} curv={:f}\n",
                       i,
                       std::abs(theta_min),
                       gdot(m_g1[g_] - m_g0[g_], inout[ax_]) / m_opt.delta_r);
          }

          if (std::abs(theta_min) < m_opt.theta_tol) {
            return gdot(m_g1[g_] - m_g0[g_], inout[ax_]) / m_opt.delta_r;
          }
        }
      }
    }();

    m_delta[del_] = -m_delta_prev[del_];  // Flip to reverse last.
    nl.update(m_delta);                   // Reset neighbour lists to input state.

    if (!m_opt.relax_in_convex && m_curv > 0) {
      out[g_] = -gdot(m_g0[g_], inout[ax_]) * inout[ax_];
    } else {
      out[g_] = m_g0[g_] - 2 * gdot(m_g0[g_], inout[ax_]) * inout[ax_];
    }

    if (count_frozen == 0) {
      trans_proj<PotentialGradient>(out);
    }

    return m_curv;
  }

}  // namespace fly::saddle