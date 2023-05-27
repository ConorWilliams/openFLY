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

#include "libfly/utility/spline.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <cstddef>

namespace fly {

  Spline::Spline(std::vector<double> y, double dx) : m_dx(dx), m_inv_dx(1 / dx) {
    // Soft pad end of spline
    double last = y.back();

    for (size_t i = 0; i < 5; i++) {
      y.push_back(last);
    }

    std::size_t n = y.size() - 1;

    // 1
    std::vector<double> a(n + 1);
    // 2
    std::vector<double> b(n);
    std::vector<double> d(n);
    // 4
    std::vector<double> alpha(n);
    // 5
    std::vector<double> c(n + 1);
    std::vector<double> l(n + 1);
    std::vector<double> mu(n + 1);
    std::vector<double> z(n + 1);

    // 1
    for (std::size_t i = 0; i <= n; ++i) {
      a[i] = y[i];
    }

    // 3
    // h_i = _dx

    // 4
    for (std::size_t i = 1; i <= n - 1; ++i) {
      alpha[i] = 3 * (a[i + 1] - 2 * a[i] + a[i - 1]) / dx;
    }

    // 6
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    // 7
    for (std::size_t i = 1; i <= n - 1; ++i) {
      l[i] = dx * (4 - mu[i - 1]);
      mu[i] = dx / l[i];
      z[i] = (alpha[i] - dx * z[i - 1]) / l[i];
    }

    // 8
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    // 9
    for (std::size_t i = 1; i <= n; ++i) {
      std::size_t j = n - i;

      c[j] = z[j] - mu[j] * c[j + 1];
      b[j] = (a[j + 1] - a[j]) / dx - dx * (c[j + 1] + 2 * c[j]) / 3;
      d[j] = (c[j + 1] - c[j]) / (3 * dx);
    }

    // 11
    for (std::size_t i = 0; i <= n - 1; ++i) {
      m_spines.push_back({a[i], b[i], c[i], d[i]});
    }

    // fmt::print("final spine: a={}, b={}, c={}, d={}\n",
    //            m_spines.back().a,
    //            m_spines.back().b,
    //            m_spines.back().c,
    //            m_spines.back().d);

    // ASSERT(m_spines.back().b >= 0, "Spline must be monotonic increasing but b={}.", m_spines.back().b);
    // m_spines.back().c = std::max(m_spines.back().b, 1.0);  //
    // m_spines.back().d = 0;
  }

}  // namespace fly