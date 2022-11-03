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

#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"

/**
 * \file basin.hpp
 *
 * @brief Representation of a set state in a Markov chain
 */

namespace fly::kinetic {

  /**
   * @brief A result type to indicate a randomly selected mechanism in a basin/superbasin.
   */
  struct Choice {
  public:
    /**
     * @brief True if the chosen mechanism starts from a basin different than the current one.
     */
    bool basin_changed() const noexcept { return m_basin_changed; }

  private:
    friend class Basin;

    bool m_basin_changed;  ///< True if mech starts from different basin.
    double m_delta_t;      ///< Time increment if mech is carried out.

    std::size_t m_mech;   ///< i'th mech in current basin.
    std::size_t m_basin;  ///< Current basin in superbasin
  };

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
      double temp;                                              ///< Of simulation in kelvin.
      double max_barrier = std::numeric_limits<double>::max();  ///< Maximum energy barrier, higher energy mechanisms are ignored.
    };

    // Basin(Options const &, Supercell const &, std::vector<Catalogue::pointer> const &);

    // // Uses standard n-fold-way kmc algorithm to select mechanism inside basin. Choice.basin_changed
    // // is always false,
    // Choice kmc_Choice(pcg64 &, std::size_t basin) const;

    // LocalMech &operator[](std::size_t i) {
    //   CHECK(i < size(), "Invalid mech");
    //   return _mechs[i];
    // }

    // LocalMech const &operator[](std::size_t i) const {
    //   CHECK(i < size(), "Invalid mech");
    //   return _mechs[i];
    // }

    // std::size_t size() const { return _mechs.size(); }

    // VecN<double> const &state() const { return _state; }

    // double rate_sum() const { return _rate_sum; }

  private:
    // Represents generalised-mechanism acting on a SPECIFIC atom.
    // class LocalMech {
    // public:
    //   double rate;            // (Hz)
    //   double barrier;         // Maximum of forward/reverse barriers.
    //   bool exit_mech = true;  // If true it is known *this links: inside SB -> outside SB.

    //   LocalMech(double rate, double barrier, std::size_t atom_idx, std::size_t mech_off, Catalogue::pointer env)
    //       : rate(rate), barrier(barrier), _env(env), _atom_idx(atom_idx), _mech_off(mech_off) {}

    //   // Reconstruct mechanism pointed to by *this onto supercell, returns ref to mech
    //   // reconstructed.
    //   Mechanism const &onto(Supercell &, std::vector<Geometry> &) const;

    //   void refine(std::vector<Geometry> &geo) const;

    // private:
    //   Catalogue::pointer _env;
    //   std::size_t _atom_idx;
    //   std::size_t _mech_off;
    // };

    bool connected = false;  // True if any of the mechs have exit_mech = false
    system::SoA<Position> m_state{10};
    std::vector<int> m_mechs;
    double m_rate_sum = 0;
  };

}  // namespace fly::kinetic