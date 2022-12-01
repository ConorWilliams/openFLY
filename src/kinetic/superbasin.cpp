// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/kinetic/superbasin.hpp"

#include <fmt/core.h>
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

  // Compute vector of mean residence-times in each basin
  Eigen::Vector<double, Eigen::Dynamic> SuperBasin::compute_tau() const {
    //
    using V = Eigen::Vector<double, Eigen::Dynamic>;

    // Non-allocating Eigen3-objects: theta_{i} = Kroneker_{is}
    auto theta = V::NullaryExpr(fly::ssize(m_super), [o = safe_cast<Eigen::Index>(m_occupied)](auto i) { return i == o; });

    auto identity = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(fly::ssize(m_super), fly::ssize(m_super));

    // Calc theta^{sum}
    V tau = (identity - m_prob).colPivHouseholderQr().solve(theta);

    // Convert to tau_i
    for (std::size_t i = 0; i < size(); ++i) {
      tau[safe_cast<Eigen::Index>(i)] /= m_super[i].rate_sum();
    }

    return tau;
  }

  SuperBasin::Choice SuperBasin::kmc_choice(Xoshiro &psudo_rng) const {
    if (!m_super[m_occupied].connected) {
      Basin::Choice ch = m_super[m_occupied].kmc_choice(psudo_rng);
      return {ch.mech, ch.atom, ch.dt, m_occupied, false};
    }

    //
    std::uniform_real_distribution<double> uniform{0, 1};

    auto tau = compute_tau();

    fmt::print("SuperBasin: density={::.3f}\n", tau.cwiseAbs() / tau.sum());

    int count = 0;

    // Sum over all basin->escape rate times basin modifiers, omit normalising factor of 1/tau.sum()
    double const r_sum = [&]() {
      //
      double sum = 0;

      for (std::size_t i = 0; i < size(); ++i) {
        //
        double basin_exit_sum = 0;

        for (auto const &mech : m_super[i].m_mechs) {
          if (mech.m_exit_mech) {
            basin_exit_sum += mech.m_rate;
            ++count;
          }
        }

        sum += tau[safe_cast<Eigen::Index>(i)] * basin_exit_sum;
      }

      return sum;
    }();

    ASSERT(count > 0, "count = {}", count);
    ASSERT(r_sum > 0, "r_sum = {}", r_sum);

    auto [basin, exit_mech] = [&]() -> std::pair<std::size_t, Basin::LocalisedMech const &> {
      //
      double const lim = uniform(psudo_rng) * r_sum;
      double sum = 0;

      for (std::size_t i = 0; i < size(); ++i) {
        for (auto const &mech : m_super[i].m_mechs) {
          if (mech.m_exit_mech) {
            sum += tau[safe_cast<Eigen::Index>(i)] * mech.m_rate;
            if (sum > lim) {
              return {i, mech};
            }
          }
        }
      }
      throw error("Hit end of super choice");
    }();

    // // Must occupy new basin (this is why method is not const)
    // std::size_t old_basin = std::exchange(m_occupied, basin);

    double const inv_tau = 1 / tau.sum();
    double const eff_rate = tau[safe_cast<Eigen::Index>(basin)] * inv_tau * exit_mech.m_rate;
    double const prob = 100 * eff_rate / (inv_tau * r_sum);

    dprint(m_opt.debug,
           "SuperBasin: SKMC choice @atom={} ΔE={:.3f}eV : {:.3f}% of {}\n",
           exit_mech.m_atom_index,
           exit_mech.m_mech->barrier,
           prob,
           count);

    // Must normalize by inv_tau
    return {
        *exit_mech.m_mech,
        exit_mech.m_atom_index,
        -std::log(uniform(psudo_rng)) / (r_sum * inv_tau),
        basin,
        m_occupied != basin,
    };
  }

  auto SuperBasin::connect_from(std::size_t basin, int atom, env::Mechanism const &m) -> void {
    //
    ASSERT(basin < m_super.size(), "Basin with index {} is OOB.", basin);

    std::vector<Basin::LocalisedMech> &mechs = m_super[basin].m_mechs;

    auto it = std::lower_bound(mechs.begin(), mechs.end(), atom, [](Basin::LocalisedMech const &elem, int val) {
      //
      return elem.m_atom_index < val;
    });

    for (; it != mechs.end() && it->m_atom_index == atom; ++it) {
      if (&m == it->m_mech) {
        m_prob(safe_cast<Eigen::Index>(m_occupied), safe_cast<Eigen::Index>(basin)) = it->m_rate / m_super[basin].rate_sum();
        ASSERT(it->m_exit_mech == true, "Chose an exit mech?", 0);
        it->m_exit_mech = false;
        m_super[basin].connected = true;
        return;
      }
    }

    throw error("Mech/atom not in basin?");
  }

  auto SuperBasin::find_occupy(std::size_t hash, system::SoA<Position const &> in, double tol) -> std::optional<std::size_t> {
    //
    auto com = [](system::SoA<Position const &> x) -> Vec {
      Vec sum = Vec::Zero();
      for (Eigen::Index i = 0; i < x.size(); i++) {
        sum += x(r_, i);
      }
      return sum / x.size();
    };

    Vec com_x = com(in);

    for (std::size_t i = 0; i < size(); ++i) {
      if (hash == m_super[i].m_state_hash) {
        //
        Vec delta = com_x - com(m_super[i].state());

        double sum_sq = 0;

        verify(in.size() == m_super[i].state().size(), "Number of atoms has changed");

        for (Eigen::Index j = 0; j < in.size(); j++) {
          sum_sq += gnorm_sq(in(r_, j) - m_super[i].state()(r_, j) - delta);
        }

        dprint(m_opt.debug, "Super: Hashes match at @{}: dr={}, drift={}\n", i, std::sqrt(sum_sq), gnorm(delta));

        if (sum_sq < tol * tol) {
          return std::exchange(m_occupied, i);
        } else {
          fmt::print("Matching hash but different basin? Check this!\n");
        }
      }
    }
    return std::nullopt;
  }

}  // namespace fly::kinetic