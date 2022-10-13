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
#include "libfly/system/VoS.hpp"
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
   * @brief A local environment is a (localised) geometry augmented with a key and a fingerprint.
   *
   * Here localised means the first atom is the central (i.e. fixed) atom and the centroid is the origin.
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
     * @brief An ordered representation of the intra-atomic distances in this Geometry.
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
      bool equiv(Fingerprint const& other, double delta) const;

      /**
       * @brief Get the smallest intra-atomic separation in this environment.
       */
      double r_min() const {
        ASSERT(!m_r_ij.empty(), "Not enough atoms for r_min!", 0);
        return m_r_ij[0];
      }

    private:
      friend class Local;

      std::vector<double> m_r_0j;
      std::vector<double> m_r_ij;
    };

    /**
     * @brief Get a constant reference to the discrete key.
     */
    Key const& key() const noexcept { return m_key; }

    /**
     * @brief Get a constant reference to the fingerprint.
     */
    Fingerprint const& fingerprint() const noexcept { return m_fingerprint; }

    /**
     * @brief Rebuild this local environment.
     *
     * @param i The index of the atom to centre this LE on.
     * @param atoms The per-atom quantities required.
     * @param nl A neighbour list in the ready state.
     * @param r_env Radius of the local environment.
     * @param r_edge Maximum distance for atoms to be considered connected in neighbour graph.
     * @param num_types The total number of types in ``atoms``.
     */
    void rebuild(int i,
                 system::SoA<TypeID const&, Frozen const&> atoms,
                 neigh::List const& nl,
                 Eigen::Index num_types,
                 double r_env,
                 double r_edge);

  private:
    Key m_key;
    Fingerprint m_fingerprint;

    void build_geo(int i,
                   system::SoA<TypeID const&, Frozen const&> atoms,
                   neigh::List const& nl,
                   Eigen::Index num_types,
                   double r_env);
  };

  //   /**
  //    * @brief This class builds and stores a list of local environments from a SimCell.
  //    */
  //   class EnvCell {
  //   public:
  //     struct Options {
  //       /** @brief Radius of environment. */
  //       double r_env;
  //     };

  //     /**
  //      * @brief Construct a new Env Cell object.
  //      */
  //     EnvCell(Options const& opt, OrthoSimBox const& box) : m_nl(box, opt.r_env){};

  //     /**
  //      * @brief Get the number of Environments in the cell.
  //      */
  //     std::size_t size() const noexcept { return m_envs.size(); }

  //     /**
  //      * @brief Get the ith local environment
  //      */
  //     Local const& operator[](std::size_t i) const noexcept { return m_envs[i]; }

  //     /**
  //      * @brief Get the ith local environment
  //      */
  //     Local& operator[](std::size_t i) noexcept { return m_envs[i]; }

  //     /**
  //      * @brief Construct the local environemnts around all the atom.
  //      */
  //     void rebuild(SimCell const& atoms, std::size_t num_threads);

  //   private:
  //     std::vector<Local> m_envs;

  //     neighbour::List m_nl;
  //   };

}  // namespace fly::env
