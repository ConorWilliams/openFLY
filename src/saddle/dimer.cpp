// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/saddle/dimer.hpp"

#include "libfly/io/gsd.hpp"

namespace fly::saddle {

  auto Dimer::step(system::SoA<Position &, Axis &> out,
                   system::SoA<Position const &, Axis const &, TypeID const &, Frozen const &> in,
                   potential::Generic &pot,
                   int max_steps,
                   int num_threads) -> Exit {
    // Check inputs

    verify(in.size() == out.size(), "Dimer stepper inputs size mismatch, in={} out={}", in.size(), out.size());

    verify(max_steps > 0, "Doing {} steps with Dimer is probably a mistake?", max_steps);

    m_eff_grad.destructive_resize(in.size());

    m_core.clear();  // Clear history from previous runs.

    // Avoid floating point comparison warnings.
    constexpr std::equal_to<> eq;

    double skin = std::max(std::pow(m_opt.skin_frac, 1. / 3.) - 1, 0.0) * m_rotor.r_cut(pot);

    if (m_opt.debug) {
      fmt::print("Dimer: Skin = {}\n", skin);
    }

    if (double r_cut = m_rotor.r_cut(pot) + skin; !eq(std::exchange(m_r_cut, r_cut), r_cut) || !m_nl) {
      m_nl = neigh::List(m_box, r_cut);

      if (m_opt.debug) {
        fmt::print("Dimer: Reallocating neigh list\n");
      }
    }

    out[r_] = in[r_];    // Initialise out = in
    out[ax_] = in[ax_];  // Initialise out = in

    m_nl->rebuild(out, num_threads);

    double curv = m_rotor.eff_gradient(m_eff_grad, out, in, pot, *m_nl, num_threads);

    double trust = m_opt.min_trust;

    double acc = 0;
    double convex_count = 0;

    for (int i = 0; i < max_steps; ++i) {
      //
      double mag_g = gnorm_sq(m_eff_grad[g_]);

      if (m_opt.debug) {
        fmt::print("Dimer: i={:<4} trust={:f} acc={:f} norm(g)={:e} curv={:f} rebuild={}\n",
                   i,
                   trust,
                   acc,
                   std::sqrt(mag_g),
                   curv,
                   eq(acc, 0.0));
      }

      if (m_opt.fout) {
        m_opt.fout->commit([&] {
          m_opt.fout->write(r_, out);
          m_opt.fout->write(ax_, out);
        });
      }

      if (mag_g < m_opt.f2norm * m_opt.f2norm) {
        return success;
      } else if (convex_count >= m_opt.convex_max) {
        if (m_opt.debug) {
          fmt::print("Dimer: Early exit - curvature\n");
        }
        return fail;
      }

      auto &Hg = m_core.newton_step<Position, PotentialGradient>(out, m_eff_grad);

      ASSERT(gdot(m_eff_grad[g_], Hg[del_]) > 0, "Ascent direction in Dimer", 0);

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

      curv = m_rotor.eff_gradient(m_eff_grad, out, in, pot, *m_nl, num_threads);

      if (curv > 0) {
        convex_count += 1;
      } else {
        convex_count = 0;
      }

      double proj = gdot(m_eff_grad[g_], Hg[del_]);

      if (proj < -m_opt.proj_tol) {
        trust = std::max(m_opt.min_trust, m_opt.shrink_trust * trust);
      } else if (proj > m_opt.proj_tol) {
        trust = std::min(m_opt.max_trust, m_opt.grow_trust * trust);
      }
    }

    return iter_max;
  }

}  // namespace fly::saddle