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

#include "libfly/env/geometry.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/hessian.hpp"
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
   * @brief Coordinate saddle-point searches on a (compute) node.
   *
   * The ``Master`` coordinates saddle-point finding for a set of openMP threads, typically one ``Master`` should be
   * instantiated per compute node. Each ``Master`` stores all the threads reusable objects: minimiser, prng, etc.
   */
  class Master {
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
       *
       * If the minimally displaced atom displaces further than this then an error will be raised.
       */
      double freeze_tol = 0.01;
      /**
       * @brief Tolerance for basins to be considered distinct.
       */
      double basin_tol = 0.25;
      /**
       * @brief Tolerance for mechanisms to be considered distinct.
       */
      double mech_tol = 0.25;
      /**
       * @brief Tolerance for stationary points to be considered distinct.
       */
      double stationary_tol = basin_tol / 2;
      /**
       * @brief Maximum error between reconstructed and relaxed saddle-points.
       */
      double sp_relax_tol = 0.5;
      /**
       * @brief Absolute eigen values of the mass-weighted hessian matrix, smaller than this, are considered zero.
       */
      double hessian_eigen_zero_tol = 1e-5;
      /**
       * @brief Number of openMP threads to dispatch saddle-point finding to.
       */
      int num_threads = omp_get_max_threads();
      /**
       * @brief Maximum number of searches per environment.
       */
      int max_searches = 1000;
      /**
       * @brief Maximum number of consecutive failed searches per environment.
       */
      int max_failed_searches = 250;
      /**
       * @brief Number of simultaneous SP searches per new environment.
       */
      int batch_size = std::min(max_failed_searches, num_threads);
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
     * @brief A collection of mechanisms centred on a central atom.
     */
    struct Found {
      Index::scalar_t centre;             ///< The central atom.
      std::vector<env::Mechanism> mechs;  ///< Mechanisms on ``this->centre``.
    };

    /**
     * @brief Construct a new Master Finder object to manage ``opt.num_threads`` openMP threads.
     *
     * This will create ``opt.num_threads`` copies of each of the input parameters, one for each thread.
     *
     * @param opt The configuration options.
     * @param box Description of the simulation space.
     * @param pot The potential to use for minimisations and SP searches.
     * @param min The minimiser to relax either side of a SP.
     * @param dimer Used to find SPs.
     */
    Master(Options const& opt,
           system::Box const& box,
           potential::Generic const& pot,
           minimise::LBFGS const& min,
           saddle::Dimer const& dimer);

    /**
     * @brief Find all the mechanisms centred on the ``unknown`` geometries.
     *
     * This will recursively schedule SP searches on the slave threads.
     *
     * @param geos List of geometries centred on the atoms which the mechanisms must be centred on.
     * @param in Description of system to search in.
     *
     */
    auto find_mechs(std::vector<env::Geometry<Index>> const& geos, system::SoA<Position const&, Frozen const&, TypeID const&> in)
        -> std::vector<Found>;

  private:
    // Per thread variables
    struct ThreadData {
      potential::Generic pot;
      minimise::LBFGS min;
      saddle::Dimer dimer;
      Xoshiro prng;
      neigh::List pot_nl;
      system::Hessian hess;
    };

    // Per input geometry data
    struct Data {
      Index::scalar_t centre;                                        // Central atom
      std::vector<std::pair<Mat, std::vector<Index::scalar_t>>> tr;  // Perm + transformation.
      std::vector<env::Mechanism> mechs;                             // Unique, mechanisms
    };

    // min->sp->min data
    struct Pathway {
      system::SoA<Position> rev;
      system::SoA<Position> sp;
      system::SoA<Position> fwd;
    };

    struct Batch {
      bool collsion = false;
      std::optional<env::Mechanism> mech = {};
      system::SoA<Position, Axis> dimer;

      explicit Batch(Eigen::Index n) : dimer(n){};
    };

    std::vector<ThreadData> m_data;

    Options m_opt;
    system::Box m_box;

    int m_deg_free;           // Degree of freedom in each input
    double m_log_prod_eigen;  // the log(prod e_i) with e_i the the eigen values of the mass weighted hessian matrix.

    ///////////////////////////////////////////////////////////////////////////////////

    // Find all mechs and write to geo_data
    void find_n(env::Geometry<Index> const& geo,
                Data& geo_data,
                system::SoA<Position const&, Frozen const&, TypeID const&> in,
                neigh::List const& nl_pert);

    bool find_batch(int tot,
                    std::vector<Batch>& batch,
                    env::Geometry<Index> const& geo,
                    Data& geo_data,
                    system::SoA<Position const&, Frozen const&, TypeID const&> in,
                    neigh::List const& nl_pert,
                    std::vector<system::SoA<Position>>& sps_cache);

    /**
     * @brief Find a single mechanism
     *
     * @param in Initial basin
     * @param dimer_in_out Perturbed dimer as input and output.
     * @param hist_sp Previous saddle points.
     * @param collision If dimer exits due to rediscovering a previous SP the this is set to true.
     * @param theta_tol forwarded to Dimer::find_sp()
     * @param geo Centred of perturbation.
     */
    std::optional<env::Mechanism> find_one(system::SoA<Position const&, Frozen const&, TypeID const&> in,
                                           system::SoA<Position&, Axis&> dimer_in_out,
                                           bool& collision,
                                           env::Geometry<Index> const& geo,
                                           std::vector<system::SoA<Position>> const& hist_sp,
                                           double theta_tol);

    // Do a saddle point search
    Dimer::Exit find_sp(system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer,
                        system::SoA<Position const&> in,
                        std::vector<system::SoA<Position>> const& hist_sp,
                        double theta_tol);

    /**
     * @brief Given a saddle point produce a min->sp->min pathway
     */
    std::optional<Pathway> do_adj_min(system::SoA<Position const&, Axis const&, Frozen const&, TypeID const&> dimer,
                                      system::SoA<Position const&> in,
                                      Index::scalar_t centre);

    //////////////////////////////////////////////

    // Compute m_deg_free and m_log_prod_eigen.
    void calc_minima_hess(system::SoA<Position const&, Frozen const&, TypeID const&> in);

    std::vector<Data> process_geos(std::vector<env::Geometry<Index>> const& geos) const;

    // Perturb in-place positions around centre and write axis,
    auto perturb(system::SoA<Position&, Axis&> out,
                 system::SoA<Position const&, Frozen const&> in,
                 Index::scalar_t centre,
                 neigh::List const& nl) -> void;

    /**
     * @brief True if mech has been seen before.
     */
    bool is_new_mech(env::Mechanism const& maybe, std::vector<env::Mechanism> const& hist) const;

    /**
     * @brief Compute every symmetric version of new_mech (according to syms), if not in mechs append to mechs.
     */
    std::size_t append_syms(std::vector<std::pair<Mat, std::vector<Index::scalar_t>>> const& syms,
                            env::Mechanism const& new_mech,
                            std::vector<env::Mechanism>& mechs) const;

    // Reconstruct saddle point according to geo and relax system to saddle,
    // check we have not constructed a false SP, if we have mark mechanism as poisoned
    system::SoA<Position> recon_sp_relax(env::Geometry<Index> const& geo,
                                         env::Mechanism& m,
                                         system::SoA<Position const&, TypeID const&, Frozen const&> in);

    ThreadData& thread() { return m_data[safe_cast<std::size_t>(omp_get_thread_num())]; }

  };  // namespace fly::saddle

}  // namespace fly::saddle