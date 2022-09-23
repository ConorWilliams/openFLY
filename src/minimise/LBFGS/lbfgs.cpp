// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/minimise/LBFGS/lbfgs.hpp"

#include <fmt/core.h>
#include <fmt/os.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <utility>

#include "libfly/minimise/LBFGS/lbfgs_core.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::minimise {

  bool LBFGS::minimise(system::SoA<Position &> out,
                       system::SoA<Position const &, TypeID const &, Frozen const &> in,
                       potential::Generic &pot,
                       int num_threads) {
    // Check inputs

    verify(in.size() == out.size(), "LBFGS minimizer inputs size mismatch, in={} out={}", in.size(), out.size());

    m_grad.destructive_resize(in.size());

    // Clear history from previous runs;

    m_core.clear();

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

    out[r_] = in[r_];  // Initialise out - in;

    m_nl->rebuild(out, num_threads);

    pot.gradient(m_grad, in, *m_nl, num_threads);

    double trust = m_opt.min_trust;

    double acc = 0;
    double convex_count = 0;

    for (int i = 0; i < m_opt.iter_max; ++i) {
      //
      double mag_g = gnorm_sq(m_grad[g_]);

      if (m_opt.debug) {
        fmt::print("LBFGS: i={:<4} trust={:f} acc={:f} norm(g)={:e} rebuild={}\n", i, trust, acc, std::sqrt(mag_g), eq(acc, 0.0));
      }

      if (m_opt.fout) {
        m_opt.fout->commit([&] {
          m_opt.fout->write(fly::r_, out);  //< Write the positions of the atoms.
        });
      }

      if (mag_g < m_opt.f2norm * m_opt.f2norm) {
        return true;
      } else if (convex_count >= m_opt.convex_max) {
        return false;
      }

      auto &Hg = m_core.newton_step<Position, PotentialGradient>(out, m_grad);

      ASSERT(gdot(m_grad[g_], Hg[del_]) > 0, "Ascent direction in lbfgs", 0);

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

      pot.gradient(m_grad, in, *m_nl, num_threads);

      //   if (std::optional curv = ; curv > 0) {
      //     convex_count += 1;
      //   } else {
      //     convex_count = 0;
      //   }

      double proj = gdot(m_grad[g_], Hg[del_]);

      if (proj < -m_opt.proj_tol) {
        trust = std::max(m_opt.min_trust, m_opt.shrink_trust * trust);
      } else if (proj > m_opt.proj_tol) {
        trust = std::min(m_opt.max_trust, m_opt.grow_trust * trust);
      }
    }

    return false;
  }

}  // namespace fly::minimise
