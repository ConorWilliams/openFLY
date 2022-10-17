#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.centroid>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <functional>
#include <memory>

#include "libfly/env/geometry.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file local.hpp
 *
 * @brief Local environments.
 *
 * The local environment (LE) of an atom is the set of atoms within some neighbourhood. It is assumed in OLKMC that the mechanisms
 * accessible to an atom are completely contained-within and solely a-function-of its LE.
 */

namespace fly::env {

  /**
   * @brief An ordered representation of the intra-atomic distances in a Geometry.
   */
  class Fingerprint {
  public:
    /**
     * @brief A fast test to see if two local environments **may** be equivalent.
     *
     * Explicitly, for the LE this fingerprint represents, this function must return ``true`` if the LE is able to be
     * ``permute_onto()`` other with the same ``delta``.
     *
     *
     * @param other The other fingerprint.
     * @param delta The tolerance for the equivalence.
     * @return ``true`` If two local environments **may** be equivalent.
     * @return ``false``  If two local environments are not equivalent.
     */
    auto equiv(Fingerprint const& other, double delta) const -> bool;

    /**
     * @brief Get the smallest intra-atomic separation in this environment.
     */
    auto r_min() const -> double {
      ASSERT(!m_r_0j.empty() && !m_r_ij.empty(), "Not enough atoms for r_min!", 0);
      return std::min(m_r_0j[0], m_r_ij[0]);
    }

  private:
    friend class Local;

    std::vector<double> m_r_0j;
    std::vector<double> m_r_ij;
  };

  /**
   * @brief A local environment is a (localised) geometry augmented with a key and a ``Fingerprint``.
   *
   * Here localised means the centroid is the origin.
   *
   * The key is discrete representation of this topology.
   * The fingerprint is a continuous representation that is cheap(ish)ly comparable.
   */
  class Local : public Geometry<Index> {
  public:
    /**
     * @brief A discrete representation of this environment.
     *
     * This is a hash of the canonical graph + atom colours.
     */
    using Key = std::size_t;

    /**
     * @brief Get a constant reference to the discrete key.
     *
     * This is a hash of the canonical graph + atom colours.
     */
    auto key() const noexcept -> Key const& { return m_key; }

    /**
     * @brief Get a constant reference to the fingerprint.
     */
    auto fingerprint() const noexcept -> Fingerprint const& { return m_fingerprint; }

    /**
     * @brief Rebuild this local environment.
     *
     * This function extracts the atoms in the neighbourhood of atom ``ix`` to build the local environment. It constructs a discrete
     * representation, encoding the colours of the atoms and the structure of the graph built from the geometry. Additionally it
     * constructs two continuous fingerprints: one based on the distance from the centre; one based on the intra-atomic distances.
     *
     * @param ix The index of the atom to centre this LE on.
     * @param atoms The per-atom quantities required.
     * @param nl A neighbour list in the ready state.
     * @param r_env Radius of the local environment.
     * @param r_edge Maximum distance for atoms to be considered connected in neighbour graph.
     * @param num_types The total number of possible types in ``atoms``, equal to ``TypeMap<...>::num_types()``.
     */
    auto rebuild(int ix,
                 system::SoA<TypeID const&, Frozen const&> atoms,
                 neigh::List const& nl,
                 Eigen::Index num_types,
                 double r_env,
                 double r_edge) -> void;

  private:
    Key m_key;
    Fingerprint m_fingerprint;
  };

  /**
   * @brief This class builds and stores a list of local environments from a system.
   */
  class LocalList {
  public:
    /**
     * @brief Used to configure the ``LocalList``.
     */
    struct Options {
      /** @brief Radius of a local environment. */
      double r_env = 5.2;
      /** @brief Maximum distance for atoms to be considered connected in the canonisation neighbour graph. */
      double r_edge = 3.0;
    };

    /**
     * @brief Construct a new Env Cell object.
     */
    LocalList(Options const& opt, system::Box const& box, int num_types) : m_opt{opt}, m_nl(box, opt.r_env), m_num_types(num_types){};

    /**
     * @brief Get the number of Environments in the cell.
     */
    auto size() const noexcept -> std::size_t { return m_envs.size(); }

    /**
     * @brief Get the ``i``th local environment
     */
    auto operator[](std::size_t i) const noexcept -> Local const& { return m_envs[i]; }

    /**
     * @brief Get the ``i``th local environment
     */
    auto operator[](std::size_t i) noexcept -> Local& { return m_envs[i]; }

    /**
     * @brief Construct the local environments around all the atoms.
     */
    auto rebuild(system::SoA<Position const&, TypeID const&, Frozen const&> const& info, int num_threads) -> void;

  private:
    std::vector<Local> m_envs;
    Options m_opt;
    neigh::List m_nl;
    int m_num_types;
  };

}  // namespace fly::env
