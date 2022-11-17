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

#include "libfly/env/catalogue.hpp"
#include "libfly/saddle//find.hpp"
#include "libfly/system/supercell.hpp"

namespace fly::env {

  /**
   * @brief Return type for ``update_cat()``.
   */
  struct Update {
    bool new_envs = false;  ///< True if new environments encountered.
    bool ref_envs = false;  ///< True if environment refinement occurred.
  };

  /**
   * @brief Updates/rebuilds a catalogue from a ``cell``.
   *
   * This function repeatedly calls ``.rebuild()`` on ``cat`` and coordinates saddle-point searches using ``mast`` to find the
   * mechanisms for the newly encountered environments.
   *
   * @param mast The saddle-point/mechanism finder.
   * @param cat The catalogue to update
   * @param cell The cell to rebuild the catalogue from.
   */
  template <typename Map, typename... T>
  Update update_cat(saddle::Master& mast, env::Catalogue& cat, system::Supercell<Map, T...> const& cell, int num_threads) {
    //

    std::vector<int> ix = cat.rebuild(cell, num_threads);

    if (ix.empty()) {
      return {};
    }

    int refines = 0;

    while (true) {
      //
      std::vector<int> fails;

      fmt::print("New envs @{} with {} refines\n", ix, refines);

      std::vector found = mast.find_mechs(saddle::Master::package({ix}, cat), cell);

      /*
       * If find_mechs has failed, the failed environments must be too symmetric, we must refine them until they are less symmetric.
       */

      for (std::size_t i = 0; i < found.size(); i++) {
        if (found[i]) {
          cat.set_mechs(ix[i], found[i].mechs());
        } else {
          fails.push_back(ix[i]);
        }
      }

      if (fails.empty()) {
        return {true, refines == 0};
      } else {
        ++refines;
      }

      // Refine tolerance's
      for (auto const& f : fails) {
        auto n = cat.calc_self_syms(f).size();
        verify(n > 1, "Env @{} cannot get any less symmetric!", f);
        do {
          double new_tol = cat.refine_tol(f, cat.get_ref(f).delta_max() / 1.5);
          verify(new_tol > 1e-6, "Probably converging to a true symmetry!");
        } while (n == cat.calc_self_syms(f).size());
      }

      std::vector tmp = cat.rebuild(cell, num_threads);

      // The atom whose symmetry tolerance increased still needs to be searched alongside any atoms that now no longer match that
      // environment.
      for (auto const& elem : tmp) {
        verify(std::find(fails.begin(), fails.end(), elem) == fails.end(), "Atom #{} found twice", elem);
      }

      // tmp now stores previous fails + new environments that don't match the refined tolerance.
      tmp.insert(tmp.end(), fails.begin(), fails.end());

      ix = std::move(tmp);
    }
  };

}  // namespace fly::env