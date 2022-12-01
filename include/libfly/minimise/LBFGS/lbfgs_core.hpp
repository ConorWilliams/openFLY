#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file lbfgs_core.hpp
 *
 * @brief Exposes the Newton-step calculation part of the L-BFGS algorithm for re-use.
 */

namespace fly::minimise {

  /**
   * @brief Class for storing L-BFGS position/gradient history.
   */
  class StepLBFGS {
  public:
    /**
     * @brief Construct a new StepLBFGS object
     *
     * @param n The number of steps to hold in the history.
     */
    explicit StepLBFGS(int n) : m_n{n}, m_hist(n){};

    /**
     * @brief Reset the history.
     */
    auto clear() -> void { m_k = 0; }

  private:
    using GenericArr = system::Property<double>::array_t;

  public:
    /**
     * @brief Compute the approximate Newton-step.
     *
     * \rst
     *
     * .. important::
     *
     *    It is expected that the user of this method will set :math:`x \gets x - \alpha q` after calling this method and before the
     *    next call.
     *
     *
     * The Newton-step is:
     *
     * .. math::
     *
     *    q = \hat{H} g
     *
     * with  :math:`g` the gradient of the potential and :math:`\hat{H}` the approximate inverse Hessian matrix of the potential
     * calculated using the L-BFGS two-loop recursion. The Newton-step, :math:`q`, is such-that the optimal step towards the minimum
     * from the current position :math:`x` is:
     *
     * .. math::
     *
     *    x \gets x - \alpha q
     *
     * with :math:`\alpha = 1` the best guess in the absence of a line-search.
     *
     * This function accepts generic ``SoA``'s to allow exotic notions of position/gradient.
     *
     * \endrst
     *
     * @param g The current 'potential gradient'.
     * @param x The current 'position' of the atoms.
     * @return A view of the approximate Newton-step (see above) (it is OK to modify this view it will be overwritten upon next call).
     */
    template <typename T1, typename T2>
    auto newton_step(system::SoA<T1 const &> x, system::SoA<T2 const &> g) -> system::SoA<Delta> & {
      static_assert(std::is_same_v<typename T1::matrix_t, typename T2::matrix_t>, ".newton_step()'s r and g must be compatible");
      static_assert(std::is_same_v<decltype(x[T1{}]), GenericArr const &>, ".newton_step() will allocate");
      static_assert(std::is_same_v<decltype(g[T2{}]), GenericArr const &>, ".newton_step() will allocate");
      return newton_step_impl(x[T1{}], g[T2{}]);
    }

  private:
    auto newton_step_impl(GenericArr const &r, GenericArr const &g) -> system::SoA<Delta> &;

    int m_n;
    int m_k = 0;

    struct Elem {
      GenericArr s;
      GenericArr y;
      double rho;
      double alpha;
    };

    Vector<Elem> m_hist;

    GenericArr m_prev_x;
    GenericArr m_prev_g;

    GenericArr m_q;

    system::SoA<Delta> m_r;
  };

}  // namespace fly::minimise