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
#include <memory>

#include "libfly/neigh/list.hpp"
#include "libfly/potential/EAM/data.hpp"
#include "libfly/potential/base.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file eam.hpp
 *
 * @brief EAM potential implementation.
 */

namespace fly::potential {

  /**
   * @brief EAM potential child class.
   */
  class EAM : public Base {
  public:
    /**
     * @brief Construct a new EAM object.
     *
     * @param data Tabulated EAM data.
     */
    explicit EAM(std::shared_ptr<DataEAM const> data) : m_data(std::move(data)) {}

    /**
     * @brief Get this potentials cut-off radius.
     *
     * This is the maximum distance two atom can interact. The neighbour::List passed to the other functions should be configured with
     * a cut-off equal or greater than this.
     */
    auto r_cut() const noexcept -> double override;

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
    auto energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, std::size_t threads = 1) -> double override;

    /**
     * @brief Compute potential energy gradient.
     *
     * Assumes the neighbour list are ready, force on frozen atoms will be zero.
     *
     * @param inout Per-atom TypeID's and Frozen properties, potential gradient written to this.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     */
    auto gradient(system::SoA<TypeID const&, Frozen const&, PotentialGradient&> inout, neigh::List const& nl, std::size_t threads = 1)
        -> double override;

  private:
    std::shared_ptr<DataEAM const> m_data;

    struct Fprime : MemTag<double, 1> {};
    struct Rho : MemTag<double, 1> {};
    struct Mu : MemTag<double, spatial_dims> {};
    struct Hidx : MemTag<std::size_t, 1> {};

    AtomArray<Fprime, Rho, Mu, Hidx> m_aux;
  };

}  // namespace fly::potential