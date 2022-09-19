
#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <string>
#include <vector>

#include "libfly/utility/core.hpp"

/**
 * \file spline.hpp
 *
 * @brief Utility for working with natural cubic splines.
 */

namespace fly {

  /**
   * @brief Computes a set of natural cubic spline coefficients for a uniformly tabulated function.
   *
   * See Wikipedia algorithm: https://en.wikipedia.org/wiki/Spline_(mathematics)
   */
  class Spline {
  public:
    /**
     * @brief Construct an empty Spline
     */
    Spline() = default;

    /**
     * @brief Construct a new Spline object.
     *
     * From (n + 1) evenly spaced tabulated values on the interval {0, dx, ..., ndx}.
     *
     * @param y Tabulated values of the function ``f``.
     * @param dx Tabulation interval.
     */
    Spline(std::vector<double> y, double dx);

    /**
     * @brief Interpolate tabulated function.
     *
     * @param x Point to interpolate ``f``.
     * @return ``f(x)`` The interpolated value of the function at ``x``.
     */
    auto f(double x) const -> double {
      auto [x0, spine] = fetch(x);
      return spine.a + x0 * (spine.b + x0 * (spine.c + x0 * spine.d));
    }

    /**
     * @brief Interpolate tabulated function's gradient.
     *
     * @param x Point to interpolate ``f'``.
     * @return ``f'(x)`` The interpolated gradient of the function at ``x``.
     */
    double fp(double x) const {
      auto [x0, spine] = fetch(x);
      return spine.b + x0 * (2.0 * spine.c + x0 * 3.0 * spine.d);
    }

    /**
     * @brief Interpolate tabulated function's second derivative.
     *
     * \rst
     *
     * .. note::
     *    This may not be continuous or smooth.
     *
     * \endrst
     *
     * @param x Point to interpolate ``f''``.
     * @return ``f''(x)`` The interpolated second derivative of the function at ``x``.
     */
    auto fpp(double x) const -> double {
      auto [x0, spine] = fetch(x);
      return 2 * spine.c + 6.0 * x0 * spine.d;
    }

  private:
    struct Spine {
      double a, b, c, d;
    };

    std::vector<Spine> m_spines;

    double m_dx = 0;
    double m_inv_dx = 0;

    std::pair<double, Spine const&> fetch(double x) const {
      //
      ASSERT(x >= 0, "x={}, not less than zero", x);

      auto i = static_cast<std::size_t>(x * m_inv_dx);

      // Could clamp this value:  i = std::min(i, m_spines.size() - 1)

      ASSERT(i < m_spines.size(), "x={} is outside tabulated region with i={}, len={}", x, i, m_spines.size());

      return {x - static_cast<double>(i) * m_dx, m_spines[i]};
    }
  };

}  // namespace fly