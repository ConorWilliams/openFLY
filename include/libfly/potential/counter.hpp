#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <atomic>
#include <cstddef>
#include <memory>

#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file counter.hpp
 *
 * @brief A potential wrapper.
 */

namespace fly::potential {

  /**
   * @brief Counter potential child class.
   */
  template <typename P>
  class Counter : private P {
  public:
    //
    using P::P;

    using P::energy;

    using P::hessian;
    using P::r_cut;

    /**
     * @brief Compute potential energy gradient.
     *
     * Can assume the neighbour list are ready, force on frozen atoms must be zero.
     *
     *
     * @param out Result is written here.
     * @param in Per-atom data used by potential for computation.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called) configured with a cut-off
     * at least ``r_cut()``.
     * @param threads Number of openMP threads to use.
     */
    auto gradient(system::SoA<PotentialGradient&> out,
                  system::SoA<TypeID const&, Frozen const&> in,
                  neigh::List const& nl,
                  int threads = 1) -> void {}

  private:
    std::shared_ptr<std::atomic<int>> m_count = std::make_shared<std::atomic<int>>(0);
  };

}  // namespace fly::potential