#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see https :  //
// www.gnu.org/licenses/>.
#include <fmt/core.h>
#include <omp.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <random>
#include <vector>

#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file find.hpp
 *
 * @brief Class to coordinate a group of saddle-point searches.
 */

namespace fly::saddle {

  /**
   * @brief Coordinate saddle-point searches on a node.
   *
   * The ``MasterFinder`` coordinates saddle-point finding for a set of openMP threads, typically one ``MasterFinder`` should be
   * instantiated per compute node. Each ``MasterFinder`` stores all the threads reusable objects: minimiser, prng, etc.
   *
   */
  class MasterFinder {
  public:
    /**
     * @brief Configure the saddle-point finding algorithm
     */
    struct Options {
      /**
       * @brief The cut-off radius of the random perturbation (centred on an atom).
       *
       * The initial system randomly perturbed to produce a starting-point for saddle-point finding. The perturbation is Gaussian in
       * each coordinate axis with standard deviation ``stddev``. An envelope function will linearly decrease the size of each atoms
       * perturbation based off its distance to a central atom. Only atoms closer than ``r_pert`` will be perturbed.
       */
      double r_pert = 4;
      /**
       * @brief The standard deviation of the random perturbations, see ``r_pert``.
       */
      double stddev = 0.6;
      /**
       * @brief Tolerance for minimally displaced atom to be frozen (to remove translational error accumulation).
       */
      double freeze_tol = 0.01;
      /**
       * @brief Tolerance for basins to be considered distinct.
       */
      double basin_tol = 0.25;

      /**
       * @brief Tolerance for stationary points to be considered distinct.
       */
      double stationary_tol = basin_tol / 2;
      /**
       * @brief Number of openMP threads to dispatch saddle-point finding to.
       */
      int num_threads = omp_get_max_threads();
      /**
       * @brief Maximum number of searches per environment.
       */
      int max_searches = 75;
      /**
       * @brief Maximum number of consecutive failed searches per environment.
       */
      int max_failed_searches = 50;

      /**
       * @brief Number of simultaneous SP searches per new environment.
       */
      int batch_size = 10;
      /**
       * @brief Fraction of distance between initial position and SP that dimer is displaced along its axis before minimisation.
       */
      double nudge_frac = 0.02;
      /**
       * @brief Print out debugging info.
       */
      bool debug = false;

      /**
       * @brief If provided write debugging data here.
       *
       * It is the users responsibility to ensure the lifetime of ``fout`` is at least as long as the lifetime of the this object
       * and ensure only a single thread writes to ``fout`` at any one time. This really only exist for debugging purposes...
       */
      io::BinaryFile* fout = nullptr;
    };

    /**
     * @brief Construct a new Master Finder object to manage ``opt.num_threads`` openMP threads.
     *
     * This will create ``opt.num_threads`` copies of each of the input parameters, one for each thread.
     *
     * @param opt The configuration options.
     * @param pot The potential to use for minimisations and SP searches.
     * @param min The minimiser to relax either side of a SP.
     * @param dimer Used to find SPs.
     */
    MasterFinder(Options const& opt, potential::Generic const& pot, minimise::LBFGS const& min, saddle::Dimer const& dimer);

    /**
     * @brief Store the SP and final minima of a min->sp->min pathway.
     */
    struct Pathway {
      system::SoA<Position> sp;   ///< The SP.
      system::SoA<Position> fwd;  ///< The forward minima.

      /**
       * @brief Construct a new Pathway object, forwards arguments to ``sp`` and ``fwd``.
       */
      template <typename U, typename V>
      Pathway(U&& u, V&& v) : sp(std::forward<U>(u)), fwd(std::forward<V>(v)) {}
    };

    /**
     * @brief Store a group of pathways all containing mechanisms centred on atom ``index``.
     */
    struct PathGroup {
      int index;                   ///< Index of central atom
      std::vector<Pathway> paths;  ///< All pathways containing mechs centred on atom `index'
    };

    /**
     * @brief Find all the mechanisms centred on the ``unknown`` atoms.
     *
     * @param box Description of the simulation space.
     * @param unknown List of indices of the atoms which the mechanisms must be centred on.
     * @param in Description of system to search in.
     *
     */
    auto find_pathways(system::Box const& box,
                       std::vector<int> const& unknown,
                       system::SoA<Position const&, Frozen const&, TypeID const&> in) -> std::vector<PathGroup>;

  private:
    struct ThreadData {
      potential::Generic pot;
      minimise::LBFGS min;
      saddle::Dimer dimer;
      Xoshiro prng;
    };

    std::vector<ThreadData> m_data;

    Options m_opt;

    /**
     * @brief Attempt to add maybe to found -> true if successful.
     */
    bool push_if_new(std::vector<Pathway>& found, Pathway&& maybe) const noexcept;

    // Find the indices of minimally and maximally separated atoms.
    static std::pair<int, int> min_max(system::SoA<Position const&> a, system::SoA<Position const&> b);

    // Find a single min-SP-min pathway
    // NL must be ready at r_pert.
    std::optional<Pathway> find_one(system::SoA<Position const&, Frozen const&, TypeID const&> in, neigh::List const& nl, int index);

    // Perturb in-place positions around centre and write axis,
    auto perturb(system::SoA<Position&, Axis&, Frozen const&> inout, int centre, neigh::List const& nl) -> void;

    ThreadData& thread() { return m_data[safe_cast<std::size_t>(omp_get_thread_num())]; }

  };  // namespace fly::saddle

}  // namespace fly::saddle