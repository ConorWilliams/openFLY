// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/minimise/LBFGS/lbfgs_core.hpp"

#include <fmt/core.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::minimise {

  system::SoA<Delta> &StepLBFGS::newton_step_impl(GenericArr const &r, GenericArr const &g) {
    //
    verify(r.size() == g.size(), "LBFGS core inputs size mismatch, r={} g={}", r.size(), g.size());
    //
    int prev = (m_k - 1) % m_n;

    m_r.destructive_resize(r.size() / spatial_dims);

    // Compute the k-1 th y, s and rho
    if (m_k > 0) {
      //
      if (r.size() != m_prev_x.size()) {
        throw error("Call to newton_step() with {} atoms but history has {} atoms", r.size(), m_prev_x.size());
      }
      //
      m_hist[prev].s = r - m_prev_x;
      m_hist[prev].y = g - m_prev_g;

      // If Wolfie conditions fulfilled during the line search then dot(y, s) > 0. Otherwise we take
      // absolute value to prevent ascent direction.
      m_hist[prev].rho = 1.0 / std::abs(gdot(m_hist[prev].s, m_hist[prev].y));
    }

    m_prev_x = r;
    m_prev_g = g;

    m_q = g;

    int incur = m_k <= m_n ? 0 : m_k - m_n;
    int bound = m_k <= m_n ? m_k : m_n;

    // Loop 1: over most recent steps
    for (int i = bound - 1; i >= 0; --i) {
      int j = (i + incur) % m_n;

      m_hist[j].alpha = m_hist[j].rho * gdot(m_hist[j].s, m_q);
      m_q -= m_hist[j].alpha * m_hist[j].y;
    }

    // Scaling Hessian_0.
    if (m_k > 0) {
      m_r[del_] = m_q * (1.0 / m_hist[prev].rho / gdot(m_hist[prev].y, m_hist[prev].y));
    } else {
      m_r[del_] = m_q;  // Start with identity hessian.
    }

    // Loop 2:
    for (int i = 0; i <= bound - 1; ++i) {
      int j = (i + incur) % m_n;

      double beta = m_hist[j].rho * gdot(m_hist[j].y, std::as_const(m_r)[del_]);

      m_r[del_] += (m_hist[j].alpha - beta) * m_hist[j].s;
    }

    ++m_k;

    return m_r;
  }

}  // namespace fly::minimise