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

#include "libfly/kinetic/basin.hpp"

#include <xxhash.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::kinetic {

  inline constexpr double INV_BOLTZ = 16021766340.0 / 1380649.0;  // eV^{-1}

  auto hash(Eigen::Index num_atoms, fly::env::Catalogue const &cat) -> std::size_t {
    //
    std::vector<int> buff(safe_cast<std::size_t>(num_atoms));

    for (Eigen::Index i = 0; i < num_atoms; ++i) {
      //
      buff[safe_cast<std::size_t>(i)] = cat.get_ref(safe_cast<int>(i)).cat_index();
    }

    return XXH64(buff.data(), buff.size() * sizeof(int), 0);
  }

  Basin::Basin(Options const &opt, system::SoA<Position const &> state, fly::env::Catalogue const &cat)
      : m_opt(opt), m_state(state) {
    //
    for (int i = 0; i < state.size(); ++i) {
      for (env::Mechanism const &mech : cat.get_ref(i).get_mechs()) {
        if (mech.barrier < opt.max_barrier) {
          //
          ASSERT(mech.barrier > 0, "Negative energy barrier for atom {} equal to {}", i, mech.barrier);

          double rate = mech.kinetic_pre * std::exp(mech.barrier / opt.temp * -INV_BOLTZ);

          m_mechs.push_back({i, rate, &mech});

          m_rate_sum += rate;
        }
      }
    }

    if (m_opt.debug) {
      fmt::print("Basin: built a basin with {} mechanisms, rate_sum={}\n", m_mechs.size(), m_rate_sum);
    }

    // verify(m_rate_sum > 1e7, "Check sanity");

    m_state_hash = hash(state.size(), cat);
  }

  auto Basin::kmc_choice(Xoshiro &psudo_rng) const -> Choice {
    //
    std::uniform_real_distribution<double> uniform{0, 1};

    LocalisedMech const &mech = [&]() -> LocalisedMech const & {
      double lim = uniform(psudo_rng) * rate_sum();
      double sum = 0;

      for (auto &&m : m_mechs) {
        sum += m.m_rate;
        if (sum > lim) {
          return m;
        }
      }

      throw error("Hit end of normal choice, lim={}, rate_sum={}", lim, rate_sum());
    }();

    if (m_opt.debug) {
      fmt::print("Basin: KMC choice @atom={} ΔE={:.3f}eV : {:.3f}% of {} choices\n",
                 mech.m_atom_index,
                 mech.m_mech->barrier,
                 mech.m_rate / rate_sum() * 100.,
                 m_mechs.size());
    }

    return {*mech.m_mech, mech.m_atom_index, -std::log(uniform(psudo_rng)) / rate_sum()};
  }

  auto Basin::most_likely(double tol) -> std::vector<Choice> {
    //

    std::vector<Choice> out;

    for (std::size_t i = 0; i < m_mechs.size(); ++i) {
      if (m_mechs[i].m_rate > m_rate_sum * tol) {
        out.push_back({*m_mechs[i].m_mech, m_mechs[i].m_atom_index, -1});
      }
    }

    return out;
  }

}  // namespace fly::kinetic