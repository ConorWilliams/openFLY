#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <exception>
#include <memory>

#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"

/**
 * \file base.hpp
 *
 * @brief Virtual interface class definition.
 */

namespace fly::potential {

  /**
   * @brief Thrown from fly::potentials::Base child classes if a pure virtual function is not supported.
   */
  struct unsupported : std::exception {};

  /**
   * @brief Specifies the virtual-interface for potentials in libFLY.
   */
  class Base {
  public:
    /**
     * @brief Get this potentials cut-off radius.
     *
     * This is the maximum distance two atom can interact. The neighbour::List passed to the other functions should be configured with
     * a cut-off equal or greater than this.
     */
    virtual auto r_cut() const noexcept -> double = 0;

    /**
     * @brief Compute the potential energy.
     *
     * Assumes the neighbour list are ready, ignores contributions from the frozen atoms.
     *
     * @param in Per-atom TypeID's and Frozen properties.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     * @return double The potential energy of the system of atoms.
     */
    virtual auto energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int threads = 1) -> double = 0;

    /**
     * @brief Compute potential energy gradient.
     *
     * Assumes the neighbour list are ready, force on frozen atoms will be zero.
     *
     * @param inout Per-atom TypeID's and Frozen properties, potential gradient written to this.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     */
    virtual auto gradient(system::SoA<TypeID const&, Frozen const&, PotentialGradient&> inout, neigh::List const& nl, int threads = 1)
        -> void
        = 0;

    /**
     * @brief Compute hessian matrix of the active atoms.
     *
     * Assumes the neighbour list are ready. The resulting hessian will be m by m and only include contributions from the m active
     * atoms. As hessian matrices are always symmetric this function is only required to compute the lower diagonal portion.
     *
     * @param in Input data.
     * @param out Hessian matrix to write output to.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     */
    virtual auto hessian(system::SoA<TypeID const&, Frozen const&> in, system::Hessian& out, neigh::List const& nl, int threads = 1)
        -> void
        = 0;

    /**
     * @brief Call parent destructor.
     */
    virtual ~Base() {}

  protected:
    /**
     * @brief Protected constructor as this is an interface class.
     */
    constexpr Base() noexcept = default;
  };

}  // namespace fly::potential