
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

#include "libfly/kinetic/cache.hpp"

namespace fly::kinetic {

  void SuperCache::connect_from(std::size_t basin,
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

        double prev_tol
            = std::exchange(m_opt.barrier_tol, std::max(0.0, m_opt.barrier_tol * m_opt.tol_shrink));

        dprint(m_opt.debug,
               "Dynamically adjusting barrier tolerance : {:.3f} -> {:.3f}\n",
               prev_tol,
               m_opt.barrier_tol);

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
      dprint(m_opt.debug,
             "SuperCache: LOADED CACHED SB with {} basins, cache size = {}\n",
             cached->size(),
             1 + size());

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

      dprint(m_opt.debug,
             "Dynamically adjusting barrier tolerance : {:.3f} -> {:.3f}\n",
             prev_tol,
             m_opt.barrier_tol);

      m_sb = SuperBasin{m_opt.opt_sb, {m_opt.opt_basin, cell, cat}};
      m_cache.clear();
    }
  }

}  // namespace fly::kinetic