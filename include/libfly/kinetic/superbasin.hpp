#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

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
 * @brief Representation of a collection of transient states in a Markov chain.
 */

namespace fly::kinetic {

  /**
   * @brief Manages a superbasin: a collection of low-barrier-linked basins.
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
     * @brief Get the number of basins in the SuperBasin.
     */
    auto size() const -> std::size_t { return m_super.size(); }

    /**
     * @brief Get the current basin i.e. the basin corresponding to the state of the system.
     */
    auto current_basin() const noexcept -> Basin const & { return m_super[m_occupied]; }

    // Basin const &operator[](std::size_t i) { return _super[i]; }

    /**
     * @brief Get the index of current basin i.e. the basin corresponding to the state of the system.
     */
    auto occupied() const noexcept -> std::size_t { return m_occupied; }

    /**
     * @brief Connect to mechanism 'mech' from basin 'basin' to the currently occupied basin.
     *
     * @param i
     * @param m
     */
    auto connect_from(std::size_t i, env::Mechanism const &m) -> void {
      //
      m_super[i].connected = true;

      for (auto &elem : m_super[i].m_mechs) {
        if (&m == elem.m_mech) {
          m_prob(safe_cast<Eigen::Index>(m_occupied), safe_cast<Eigen::Index>(i)) = elem.m_rate / m_super[i].rate_sum();
          elem.m_exit_mech = false;
          return;
        }
      }
      throw error("Mech not in basin?");
    }

    /**
     * @brief Expand the superbasin by adding 'basin' to it and setting 'basin' as the currently occupied basin. Returns the previously
     * occupied basin's index.
     *
     * @param basin
     * @return std::size_t
     */
    auto expand_occupy(Basin &&basin) -> std::size_t {
      m_super.push_back(std::move(basin));
      m_prob.conservativeResizeLike(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(ssize(*this), ssize(*this)));
      return std::exchange(m_occupied, size() - 1);
    }

    //

    /**
     * @brief Search through the basins in the superbasin, if one is found that matches (all active atoms within L2 tolerance 'tol')
     * make it the occupied basin. Returns the previously occupied basin's index.
     *
     * @param hash
     * @param in
     * @param tol
     * @return std::optional<std::size_t>
     */
    auto find_occupy(std::size_t hash, system::SoA<Position const &> in, double tol) -> std::optional<std::size_t>;

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
     * @brief Choose a mechanism using the modified mean-rate-method, the mechanism may start from a basin that the system is not
     * currently in, if this is the case the occupied basin is updated and choice.basin_changed is set to true.
     *
     * @param psudo_rng
     * @return Choice
     */
    auto kmc_choice(Xoshiro &psudo_rng) -> Choice;

  private:
    Options m_opt;
    std::vector<Basin> m_super{};                                    // Collection of basins.
    std::size_t m_occupied{};                                        // Index of current basin in super.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_prob{};  // Transition probability matrix.

    Eigen::Vector<double, Eigen::Dynamic> compute_tau() const;
  };

}  // namespace fly::kinetic