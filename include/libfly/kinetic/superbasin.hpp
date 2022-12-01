#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/kinetic/basin.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file superbasin.hpp
 *
 * @brief Representation of a collection of basins.
 */

namespace fly::kinetic {

  /**
   * @brief Manages a superbasin: a collection of low-barrier-linked basins. This implements a variation of bac-MRM, the basin
   * auto-constructing mean rate method, to analytically accelerate KMC.
   */
  class SuperBasin {
  public:
    /**
     * @brief Configure the ``SuperBasin`` class.
     */
    struct Options {
      bool debug = false;  ///< Controls printing of debugging data.
    };

    /**
     * @brief Construct a new SuperBasin object containing a single basin.
     */
    SuperBasin(Options const &opt, Basin &&basin) : m_opt(opt) { expand_occupy(std::move(basin)); }

    /**
     * @brief A result type to indicate a randomly selected mechanism in a basin/superbasin.
     */
    class Choice {
    public:
      env::Mechanism const &mech;  ///< A reference to the mechanism stored in the catalogue.
      int atom;                    ///< The index of the atom chosen
      double dt;                   ///< The time elapsed if this mechanisms is carried out.
      std::size_t basin;           ///< The index of the basin that ``mech`` starts from
      bool basin_changed;          ///< True if occupied basin changed during selection.
    };

    /**
     * @brief Choose a mechanism using the modified mean-rate-method.
     *
     * This will pick a non-transient mechanism, e.g. a mechanisms that "escapes" the superbasin (although it may be an internal
     * mechanisms that has not been traversed yet).
     *
     * \rst
     *
     * .. note::
     *
     *   The mechanism may start from a basin that the system is not currently in.
     *
     * \endrst
     *
     * @param pseudo_rng The pseudo random number generator to use.
     */
    auto kmc_choice(Xoshiro &pseudo_rng) const -> Choice;

    /**
     * @brief Get the number of basins in the SuperBasin.
     */
    auto size() const -> std::size_t { return m_super.size(); }

    /**
     * @brief Get the state of the Basin indexed by ``basin``.
     */
    auto state(std::size_t basin) const -> system::SoA<Position> const & {
      ASSERT(basin < m_super.size(), "Basin with index {} is OOB.", basin);
      return m_super[basin].state();
    }

    /**
     * @brief Connect mechanism ``mech`` from basin ``basin`` to the currently occupied basin.
     *
     * This marks ``mech`` as an internal, transient mechanism and fills in the internal transition probability matrix.
     *
     * @param basin The index of the basin which the mechanism starts from.
     * @param atom The index of the atom which the mechanism is centred on.
     * @param m The mechanism that connect ``basin`` to the currently occupied basin.
     */
    auto connect_from(std::size_t basin, int atom, env::Mechanism const &m) -> void;

    /**
     * @brief Expand the SuperBasin by adding ``basin`` to it and setting ``basin`` as the currently occupied basin.
     *
     * @param basin The basin to add to the SuperBasin.
     * @return Returns the previously occupied basin's index.
     */
    auto expand_occupy(Basin &&basin) -> std::size_t {
      m_super.push_back(std::move(basin));
      m_prob.conservativeResizeLike(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(fly::ssize(*this), fly::ssize(*this)));
      return std::exchange(m_occupied, size() - 1);
    }

    //

    /**
     * @brief Find and occupy a basin in the SuperBasin whose state matches ``in``.
     *
     * Searches through the basins in the superbasin, if one is found that matches (all active atoms within L2 tolerance ``tol``)
     * make it the occupied basin.
     *
     * @param hash Of ``in`` as computed by ``fly::kinetic::hash()``.
     * @param in The input state.
     * @param tol The L2 tolerance between ``in`` and a basin state for them to be considered "the same".
     * @return Returns the previously occupied basin's index or ``std::nullopt`` if no match is found.
     */
    auto find_occupy(std::size_t hash, system::SoA<Position const &> in, double tol) -> std::optional<std::size_t>;

  private:
    Options m_opt;
    std::vector<Basin> m_super{};                                    // Collection of basins.
    std::size_t m_occupied{};                                        // Index of current basin in super.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_prob{};  // Transition probability matrix.

    Eigen::Vector<double, Eigen::Dynamic> compute_tau() const;
  };

}  // namespace fly::kinetic