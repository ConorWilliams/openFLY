// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/catalogue.hpp"

namespace fly::env {

  std::optional<Catalogue::Pointer> Catalogue::canon_find(Local& mut) {
    // "it" always points to valid bucket (possibly empty)
    auto [it, inserted] = m_cat.try_emplace(mut.key());

    if (!inserted) {
      // Existing key, must search bucket for explicit match
      auto match = std::find_if(it->second.begin(), it->second.end(), [&](Env& ref) {
        //
        double r_min = std::min(mut.fingerprint().r_min(), ref.fingerprint.r_min());

        double delta = std::min(0.4 * ref.delta_mod * r_min, m_opt.delta_max);

        // Test if fuzzy keys match (fast)
        if (!ref.fingerprint.equiv(mut.fingerprint(), delta * m_opt.overfuzz)) {
          return false;
        }

        if (mut.permute_onto(ref, delta)) {
          ref.freq += 1;
          return true;
        } else {
          return false;
        }
      });

      // If found a match, return it
      if (match != it->second.end()) {
        return Pointer(it, match - it->second.begin());
      }
    }

    return std::nullopt;
  }

  Catalogue::Pointer Catalogue::insert(Local& env) {
    //

    ASSERT(!canon_find(env), "Environment already in catalogue.", 0);

    auto [it, inserted] = m_cat.try_emplace(env.key());

    auto& new_ref = it->second.emplace_back(env.fingerprint());

    for (auto&& elem : env) {
      new_ref.emplace_back(elem[r_], elem[col_]);
    }

    m_size++;

    return {it, xise(it->second) - 1};
  }

}  // namespace fly::env