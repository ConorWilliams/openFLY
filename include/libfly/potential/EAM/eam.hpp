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
#include "libfly/system/typemap.hpp"
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
     * @param map TypeMap of expected types to find in data.
     * @param data Tabulated EAM data.
     */
    EAM(system::TypeMap<> const& map, std::shared_ptr<DataEAM const> data) : m_data(std::move(data)) {
      if (map.num_types() != m_data->type_map().num_types()) {
        throw error("Different number of types in eam data file: {} != {}", map.num_types(), m_data->type_map().num_types());
      }

      for (TypeID::scalar_t i = 0; i < map.num_types(); i++) {
        std::string_view exp = m_data->type_map().get(i, tp_);
        std::string_view got = map.get(i, tp_);
        if (exp != got) {
          throw error("TypeID mismatch: expecting type {}, got {} at id={}", exp, got, i);
        }
      }
    }

    /**
     * @brief Get this potentials cut-off radius.
     *
     * This is the maximum distance two atom can interact. The neighbour::List passed to the other functions should be configured with
     * a cut-off equal or greater than this.
     */
    auto r_cut() const noexcept -> double override { return m_data->r_cut(); }

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
    auto energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int threads = 1) -> double override;

    /**
     * @brief Compute potential energy gradient.
     *
     * Assumes the neighbour list are ready, force on frozen atoms will be zero.
     *
     * @param inout Per-atom TypeID's and Frozen properties, potential gradient written to this.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     */
    auto gradient(system::SoA<TypeID const&, Frozen const&, PotentialGradient&> inout, neigh::List const& nl, int threads = 1)
        -> void override;

  private:
    std::shared_ptr<DataEAM const> m_data;

    struct Fprime : system::Property<double, 1> {};
    struct Rho : system::Property<double, 1> {};
    struct Mu : system::Property<double, spatial_dims> {};

    struct Hidx : system::Property<std::size_t, 1> {};

    system::SoA<Fprime, Rho, Mu, Hidx> m_aux;
  };

}  // namespace fly::potential