#pragma once

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

#include <fmt/core.h>

#include <cstddef>
#include <limits>
#include <list>
#include <utility>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/kinetic/basin.hpp"
#include "libfly/kinetic/superbasin.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file cache.hpp
 *
 * @brief Representation of a collection of basins.
 */

namespace fly::kinetic {

  /**
   * @brief Manages a collection of superbasins.
   *
   * A SuperCache dynamically manages a collection of superbasin to remove the overhead of exploring/building
   * a superbasin. The is particularly useful when superbasins are frequently re-entered upon exit.
   */
  class SuperCache {
  public:
    /**
     * @brief Configure the ``SuperCache`` class.
     */
    struct Options {
      double state_tol = 0.1;                 ///< L2 norm between atoms for basins to be considered the same.
      double barrier_tol = 0.3;               ///< Tolerance for mechanisms to be considered high-barrier.
      std::size_t cache_size = 128;           ///< Max number of SuperBasins in the cache.
      bool dynamic_tol = true;                ///< If true ``barrier_tol`` is dynamically adjusted.
      std::size_t max_superbasin_size = 256;  ///< Maximum number basins in SB before ``tol_shrink``.
      double tol_grow = 1.5;                  ///< Multiplier by which ``barrier_tol`` increases by.
      double tol_shrink = 0.5;                ///< Multiplier by which ``barrier_tol`` decreases by.
      bool debug = false;                     ///< Controls debug printing.
      Basin::Options opt_basin = {};          ///< Basin::Options options.
      SuperBasin::Options opt_sb = {};        ///< SuperBasin::Options options.
    };

    /**
     * @brief Construct a new SuperCache object initialised with a single superbasin.
     *
     * @param opt The options for this SuperCache.
     * @param sb The initial superbasin.
     */
    SuperCache(Options const &opt, SuperBasin &&sb) : m_opt(opt), m_sb(std::move(sb)) {}

    /**
     * @brief Get the state of the Basin indexed by ``basin`` in the active basin.
     */
    auto const &state(std::size_t basin) const { return m_sb.state(basin); }

    // Pre condition: _sb's occupied basin's state matches 'cell'.
    // Select a mechanism using appropriate KMC algorithm from active superbasin, if mechanism
    // starts from a different basin cell's active atoms are updated accordingly.

    /**
     * @brief Forward the call to the active superbasin's ``fly::kinetic::SuperBasin::kmc_choice()``.
     */
    auto kmc_choice(Xoshiro &psudo_rng) const -> SuperBasin::Choice { return m_sb.kmc_choice(psudo_rng); }

    /**
     * @brief Manages the internal superbasins.
     *
     * Identifies which (if any) basin ``cell`` matches, makes that basins's superbasin the active superbasin
     * and then connects the mechanism ``mech`` with ``fly::kinetic::SuperBasin::connect_from()`` if this was
     * a low barrier mechanism.

     * @param basin The index of the basin which the mechanism starts from in the active superbasin.
     * @param atom The index of the atom which the mechanism is centred on.
     * @param mech The mechanism that was chosen.
     * @param cell The system after reconstructing ``mech`` onto the previous state.
     * @param cat The catalogue in the ready state (i.e. called ``cat.rebuild()`` with argument ``cell``).
     */
    void connect_from(std::size_t basin,
                      int atom,
                      env::Mechanism const &mech,
                      system::SoA<Position const &> cell,
                      env::Catalogue const &cat);

    /**
     * @brief Get the number of superbasins in the SuperCache.
     */
    std::size_t size() const noexcept { return 1 + m_cache.size(); }

    /**
     * @brief Reset the cache to contain only a single superbasin ``sb``.
     */
    void reset(SuperBasin &&sb) {
      m_sb = std::move(sb);
      m_cache.clear();
      m_in_cache_count = 0;
    }

  private:
    Options m_opt;

    SuperBasin m_sb;
    std::list<SuperBasin> m_cache;

    std::size_t m_in_cache_count = 0;

    //////////////////////////////////

    // // Caches SB and removes oldest SB.
    void cache(SuperBasin &&sb) {
      //
      m_cache.push_front(std::move(sb));

      if (m_cache.size() > m_opt.cache_size) {
        m_cache.pop_back();
      }
    }
  };

}  // namespace fly::kinetic