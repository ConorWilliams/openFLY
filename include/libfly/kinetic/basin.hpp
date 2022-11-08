#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <limits>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file basin.hpp
 *
 * @brief Representation of a set state in a Markov chain
 */

namespace fly::kinetic {

  /**
   * @brief Represents ``Mechanism`` acting on a SPECIFIC atom.
   *
   */
  class LocalisedMech {
  public:
  private:
    friend class Basin;
    friend class SuperBasin;

    int m_atom_index;                   // The index of the atom in the supercell this mechanisms is centred on.
    fly::env::Mechanism const *m_mech;  // The actual mechanism.
    double m_rate;                      // (Hz)
    double m_barrier;                   // Maximum of forward/reverse barriers.
    bool m_exit_mech = true;            // If true it is thought *mech links: inside SB -> outside SB.

    /**
     * @brief Construct a new Localised Mech object.
     *
     * @param i The atom index this mech acts on.
     * @param rate The rate (frequency) of this mechanism.
     * @param mech Pointer to the mechanism.
     */
    LocalisedMech(int i, double rate, fly::env::Mechanism const *mech)
        : m_atom_index(i), m_mech(mech), m_rate(rate), m_barrier(mech->barrier - std::min(0.0, mech->delta)) {
      ASSERT(mech, "Null mechanisms is invalid", 0);
    }
  };

  /**
   * @brief Compute a hash encoding the ID's of the local environments of each atom loaded into the catalogue.
   *
   * @param num_atoms Number of atoms loaded into the catalogue.
   * @param cat A catalogue in the ready state.
   */
  auto hash(int num_atoms, fly::env::Catalogue const &cat) -> std::size_t;

  //   class SuperChoice : private Choice {
  //   public:
  //     /**
  //      * @brief True if the chosen mechanism starts from a basin different than the current one.
  //      */
  //     auto basin_changed() const noexcept -> bool { return m_basin_changed; }

  //   private:
  //     friend class SuperBasin;
  //     friend class SuperCache;

  //     bool m_basin_changed;  // True if mech starts from different basin.
  //     std::size_t m_basin;   // Current basin in superbasin
  //   };

  /**
   * @brief  Represents a basin of the potential energy of the entire system.
   *
   * Stores a reference image of the system (active atoms only) and a list of all mechanisms accessible from the basin. Implements
   * standard kmc algorithm on this list. Stored mechanisms are pointers into catalogue and can be reconstructed onto a supercell.
   *
   */
  class Basin {
  public:
    /**
     * @brief Configure the ``Basin`` class.
     */
    struct Options {
      bool debug = false;                                       ///< Controls printing of debugging data.
      double temp = 300;                                        ///< Of simulation in kelvin.
      double max_barrier = std::numeric_limits<double>::max();  ///< Maximum energy barrier, higher energy mechanisms are ignored.
    };

    /**
     * @brief Construct a new Basin.
     *
     * @param opt The options for this Basin.
     * @param state The state of the system.
     * @param cat A catalogue in the ready state for this system.
     */
    Basin(Options const &opt, system::SoA<Position const &> state, fly::env::Catalogue const &cat);

    /**
     * @brief A result type to indicate a randomly selected mechanism in a basin/superbasin.
     */
    class Choice {
    public:
      env::Mechanism const &mech;  ///< A reference to the mechanism stored in the catalogue.
      int atom;                    ///< The index of the atom chosen
      double dt;                   ///< Fetch the time elapsed if this mechanisms is carried out.
    };

    /**
     * @brief Choose a mechanisms in the basin using standard n-fold-way kmc algorithm.
     */
    auto kmc_choice(Xoshiro &psudo_rng) const -> Choice;

    /**
     * @brief Fetch the sum of the rates of all mechanisms from this basin.
     */
    double rate_sum() const noexcept { return m_rate_sum; }

    /**
     * @brief Fetch every mechanisms that is accessible that with a probability greater than ``tol``.
     */
    auto most_likely(double tol) -> std::vector<Choice>;

  private:
    friend class SuperBasin;
    friend class SuperCache;

    Options m_opt;

    bool connected = false;              // True if any of the mechs have exit_mech = false.
    double m_rate_sum = 0;               // From this basin.
    std::size_t m_state_hash;            // A hash of the state.
    system::SoA<Position> m_state;       // The full xyz state info.
    std::vector<LocalisedMech> m_mechs;  // Pointers invalidated if catalogue refined!
  };

}  // namespace fly::kinetic