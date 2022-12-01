#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

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
 * @brief d
 */

namespace fly::kinetic {

  /**
   * @brief Manages a collection of superbasins.
   */
  class SuperCache {
  public:
    /**
     * @brief Configure the ``SuperCache`` class.
     */
    struct Options {
      double state_tol = 0.1;                 ///< L2 norm between atoms for basins to be considered the same.
      double barrier_tol = 0.3;               ///< Tolerance for mechanisms to be considered high-barrier.
      std::size_t cache_size = 64;            ///< Max number of SuperBasins in the cache.
      bool dynamic_tol = true;                ///< If true ``barrier_tol`` is dynamically adjusted.
      std::size_t max_superbasin_size = 256;  ///< Maximum number basins in SB before ``tol_shrink``.
      double tol_grow = 1.5;                  ///< Multiplier by which ``barrier_tol`` increases by.
      double tol_shrink = 0.5;                ///< Multiplier by which ``barrier_tol`` decreases by.
      bool debug = false;                     ///< Controls debug printing.
      Basin::Options opt_basin = {};          ///< Basin options.
      SuperBasin::Options opt_sb = {};        ///< Superbasin options.
    };

    /**
     * @brief Construct a new Super Cache object
     *
     * @param opt
     * @param sb
     */
    SuperCache(Options const &opt, SuperBasin &&sb) : m_opt(opt), m_sb(std::move(sb)) {}

    /**
     * @brief
     *
     * @param basin
     * @return auto const&
     */
    auto const &state(std::size_t basin) const { return m_sb.state(basin); }

    // Pre condition: _sb's occupied basin's state matches 'cell'.
    // Select a mechanism using appropriate KMC algorithm from active superbasin, if mechanism
    // starts from a different basin cell's active atoms are updated accordingly.
    /**
     * @brief
     *
     * @param psudo_rng
     * @return SuperBasin::Choice
     */
    auto kmc_choice(Xoshiro &psudo_rng) const -> SuperBasin::Choice { return m_sb.kmc_choice(psudo_rng); }

    // Connect _sb's active basin to the basin 'cell' is currently in via mechanism 'mech', updates
    // active basin to match 'cell'/'env'.
    // Post condition: _sb's occupied basin's state matches 'cell'.

    /**
     * @brief
     *
     * @param basin
     * @param atom
     * @param m
     * @param cell
     * @param cat
     */
    void connect_via(std::size_t basin,
                     int atom,
                     env::Mechanism const &m,
                     system::SoA<Position const &> cell,
                     env::Catalogue const &cat) {
      //

      auto hash = kinetic::hash(cell.size(), cat);

      if (std::optional<std::size_t> found = m_sb.find_occupy(hash, cell, m_opt.state_tol)) {
        // Did internal jump -> connect states
        m_sb.connect_from(basin, atom, m);

        dprint(m_opt.debug, "SuperCache: Existing basin in SB, size={}\n", m_sb.size());

        return;
      }

      if (std::max(m.barrier, m.barrier - m.delta) < m_opt.barrier_tol) {
        // Followed low barrier to get here -> must add basin to SB.

        if (m_opt.dynamic_tol && m_sb.size() >= m_opt.max_superbasin_size) {
          // Overflowing SB -> must lower barrier_tol

          double prev_tol = std::exchange(m_opt.barrier_tol, std::max(0.0, m_opt.barrier_tol * m_opt.tol_shrink));

          dprint(m_opt.debug, "Dynamically adjusting barrier tolerance : {:.3f} -> {:.3f}\n", prev_tol, m_opt.barrier_tol);

          m_sb = SuperBasin{m_opt.opt_sb, {m_opt.opt_basin, cell, cat}};
          m_cache.clear();

        } else {
          // Can expand and occupy
          m_sb.expand_occupy({m_opt.opt_basin, cell, cat});
          m_sb.connect_from(basin, atom, m);

          dprint(m_opt.debug, "SuperCache: New basin in SB, size={}\n", m_sb.size());
        }

        return;
      }

      // Therefore followed high-barrier out of basin

      // Try and retrieve cached SB
      auto cached = [&]() -> std::optional<SuperBasin> {
        for (auto it = m_cache.begin(); it != m_cache.end(); ++it) {
          if (std::optional prev_basin = it->find_occupy(hash, cell, m_opt.state_tol)) {
            SuperBasin tmp = std::move(*it);
            m_cache.erase(it);
            return tmp;
          }
        }
        return {};
      }();

      if (cached) {
        dprint(m_opt.debug, "SuperCache: LOADED CACHED SB with {} basins, cache size = {}\n", cached->size(), 1 + size());

        cache(std::exchange(m_sb, std::move(*cached)));
        m_in_cache_count += 1;
      } else {
        dprint(m_opt.debug, "SuperCache: NEW SB, cache size = {}\n", 1 + size());

        cache(std::exchange(m_sb, SuperBasin{m_opt.opt_sb, {m_opt.opt_basin, cell, cat}}));
        m_in_cache_count = 0;
      }

      if (m_opt.dynamic_tol && m_in_cache_count > m_opt.cache_size) {
        //
        double prev_tol = std::exchange(m_opt.barrier_tol, m_opt.barrier_tol * m_opt.tol_grow);

        dprint(m_opt.debug, "Dynamically adjusting barrier tolerance : {:.3f} -> {:.3f}\n", prev_tol, m_opt.barrier_tol);

        m_sb = SuperBasin{m_opt.opt_sb, {m_opt.opt_basin, cell, cat}};
        m_cache.clear();
      }
    }

    /**
     * @brief
     *
     * @return std::size_t
     */
    std::size_t size() const noexcept { return 1 + m_cache.size(); }

    /**
     * @brief
     *
     * @param sb
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