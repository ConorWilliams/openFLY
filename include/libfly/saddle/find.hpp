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

#include "libfly/env/catalogue.hpp"
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

// clang-format off

/**
 * \file find.hpp
 *
 * @brief Class to coordinate a group of saddle-point searches.
 *
 *
 * \rst
 *
 * Saddle point (SP) searching is an involved process below is an outline of the stages:
 * 
 * 1. For each central atom ``i`` in the input the following steps are carried out.
 * 2. Saddle point searches are carried out (in batches) until a stopping criterion (one of: reconstruction failure, consecutive failed search limit or hitting a total search limit) is reached. Each SP search requires the following series of steps. 
 * 3. A random perturbation of atomic positions and the dimer axis (centred on ``i``)  is made. 
 * 4. The dimer method is used to move the perturbed state to a SP of the potential. It uses a ``cache`` of previously discovered saddle points (SPs)to avoid previous SPs by cancelling the search if the angle between the current search and the historical SP falls below ``theta``. The angle ``theta`` is a decaying function of the number of SP searches.  
 * 5. It is ensured that the maximally displaced atom is ``i`` and the minimally displaced atom is frozen (if there are translational degrees of freedom).
 * 6. The LBFGS minimiser is used to find the adjacent minima, due to freezing the atom in the previous step we can guarantee no translation in the ``min->sp->min`` pathway.
 * 7. It is ensured than the pathway starts at the initial basin and that all the stationary points are distinct points in space. 
 * 8. A local mechanism is constructed from the displacements of the atoms around the central atom. 
 * 9. It is ensured that the local mechanism can be reconstructed in the orientation and location it was discovered. If not it is marked as a poisoned mechanism.
 * 10. It is ensured that the mechanism has not been encountered before.
 * 11. The set of mechanisms related by the symmetries of the local environment are constructed and it is ensured that when these mechanisms are reconstructed the energy barriers are accurate.
 *
 * \endrst
 */

// clang-format on

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
      double basin_tol = 0.05;
      /**
       * @brief Tolerance for mechanisms to be considered distinct.
       */
      double mech_tol = 0.25;
      /**
       * @brief Tolerance for stationary points to be considered distinct.
       */
      double stationary_tol = basin_tol / 2;
      /**
       * @brief Maximum distance error between reconstructed+relaxed minima/saddle-points and their original discovery.
       *
       * After a min->sp->min pathway is found it is: converted into a mechanism; reconstructed onto the initial minima and then
       * relaxed. These relaxed SP/minima must be within ``capture_r_tol`` of the pathway's sp and final min.
       */
      double capture_r_tol = 0.05;
      /**
       * @brief Maximum energy error between reconstructed+relaxed minima/saddle-points and their original discovery.
       *
       * After a min->sp->min pathway is found it is: converted into a mechanism; reconstructed onto the initial minima and then
       * relaxed. These relaxed SP/minima must have energies within ``capture_E_tol`` of the pathway's sp and final min.
       */
      double capture_E_tol = 0.01;
      /**
       * @brief The absolute energy tolerance for a mechanisms energy barrier & delta when reconstructed onto symmetries.
       */
      double recon_e_tol_abs = 0.01;
      /**
       * @brief The fractional energy tolerance for a mechanisms energy barrier & delta when reconstructed onto symmetries.
       */
      double recon_e_tol_frac = 0.01;
      /**
       * @brief The absolute difference between the expected relaxation and the measured relaxation of a reconstructed mechanism onto
       * its symmetries.
       */
      double recon_norm_abs_tol = 0.25;
      /**
       * @brief The fractional difference between the expected relaxation and the measured relaxation of a reconstructed mechanism onto
       * its symmetries.
       */
      double recon_norm_frac_tol = 0.5;
      /**
       * @brief Absolute eigen values of the mass-weighted hessian matrix, smaller than this, are considered zero.
       */
      double hessian_eigen_zero_tol = 1e-4;
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
    class Found {
    public:
      /**
       * @brief Truthy if search was successful.
       *
       * If during the SPS, a reconstruction was attempted onto a supposedly symmetric state and the dimer-method failed to converge
       * that reconstruction to a SP then the searches in that environment will be cancelled and this will return false.
       *
       * This probably occurred because the symmetry tolerance was not tight enough e.g. this geometry was associated to the wrong
       * reference LE. The best course of action is to tighten the corresponding reference LE's ``delta_max``, reclassify
       * system and re-run the SP searches.
       */
      explicit operator bool() const noexcept { return !m_fail; }

      /**
       * @brief Fetch the mechanisms on ``this->centre``.
       *
       * Throws if SPS failed.
       */
      std::vector<env::Mechanism>& mechs() {
        if (m_fail) {
          throw error("Attempting to fetch mechanisms but the search failed!");
        } else {
          return m_mechs;
        }
      }

    private:
      friend class Master;

      bool m_fail = false;
      std::vector<env::Mechanism> m_mechs{};
    };

    /**
     * @brief Convenience alias.
     */
    using SoA = system::SoA<Position const&, Frozen const&, TypeID const&>;

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
     * @brief A structure which packages the data require to perform saddle point searches.
     */
    struct LocalisedGeo {
      Index::scalar_t centre;                         ///< Central atom.
      env::Geometry<Index> geo;                       ///< Geometry around central atom.
      std::vector<env::Catalogue::SelfSymetry> syms;  ///< Symmetries of ``geo`` onto itself.
    };

    /**
     * @brief Transform a list of atom indices, ``ix``, into a list of ``LocalisedGeo`` for use in ``find_mechs()``.
     *
     * @param ix List of indexes of environments which the mechanisms must be centred on.
     * @param cat Catalogue in ready state.
     * @param num_threads Number of openMP threads to use for this operation.
     */
    static auto package(std::vector<int> const& ix, env::Catalogue const& cat, int num_threads = 1) -> std::vector<LocalisedGeo>;

    /**
     * @brief Find all the mechanisms centred on the ``unknown`` geometries.
     *
     * This will recursively schedule SP searches on the slave threads.
     *
     * @param geos A list of localised geometries encoding the atoms to centre the SP searches on and their local environments.
     * @param in Description of system to search in.
     */
    auto find_mechs(std::vector<LocalisedGeo> const& geos, SoA in) -> std::vector<Found>;

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

    // min->sp->min data
    struct Pathway {
      system::SoA<Position> rev;
      system::SoA<Position> sp;
      system::SoA<Position> fwd;
    };

    struct Batch {
      Dimer::Exit exit = Dimer::Exit::uninit;
      std::optional<env::Mechanism> mech = {};
      system::SoA<Position, Axis> dimer;

      explicit Batch(Eigen::Index n) : dimer(n){};
    };

    std::vector<ThreadData> m_data;

    Options m_opt;
    system::Box m_box;

    int m_num_zero_modes;     // Degree of freedom in each input
    int m_count_frozen;       // The number of frozen atoms in the input.
    double m_log_prod_eigen;  // the log(prod e_i) with e_i the the eigen values of the mass weighted hessian matrix.

    ///////////////////////////////////////////////////////////////////////////////////

    // Find all mechs and write to geo_data
    void find_n(Found& out, LocalisedGeo const& geo_data, SoA in, neigh::List const& nl_pert);

    bool find_batch(int tot,
                    Found& out,
                    std::vector<Batch>& batch,
                    LocalisedGeo const& geo_data,
                    SoA in,
                    neigh::List const& nl_pert,
                    std::vector<system::SoA<Position>>& sps_cache);

    /**
     * @brief Find a single mechanism
     *
     * @param in Initial basin
     * @param dimer_in_out Perturbed dimer as input and output.
     * @param hist_sp Previous saddle points.
     * @param exit Stores the return code of the SP search.
     * @param theta_tol forwarded to Dimer::find_sp()
     * @param geo Centred of perturbation.
     */
    std::optional<env::Mechanism> find_one(SoA in,
                                           system::SoA<Position&, Axis&> dimer_in_out,
                                           Dimer::Exit& exit,
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

    /**
     * @brief Check a reconstructed mechanism is within tolerances, if it is cache the SP if not set out's fail flag and do not cache.
     */
    void check_mech(Found& out,
                    system::SoA<Position>& cache_slot,
                    env::Mechanism const& mech,
                    std::size_t sym_num,
                    LocalisedGeo const& geo_data,
                    SoA in);

    // Compute m_deg_free and m_log_prod_eigen.
    void calc_minima_hess(SoA in);

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
    std::size_t append_syms(std::vector<env::Catalogue::SelfSymetry> const& syms,
                            env::Mechanism const& new_mech,
                            std::vector<env::Mechanism>& mechs) const;

    struct Recon {
      system::SoA<Position> min;                     ///< Reconstructed minima
      system::SoA<Position> sp;                      ///< Reconstructed sp
      std::optional<system::SoA<Position>> rel_min;  ///< Relaxed reconstructed minima
      std::optional<system::SoA<Position>> rel_sp;   ///< Relaxed reconstructed sp
    };

    /**
     * @brief Reconstruct and relax a mechanism's saddle point and minima.
     */
    Recon recon_relax(env::Geometry<Index> const& geo,
                      env::Mechanism const& m,
                      system::SoA<Position const&, TypeID const&, Frozen const&> in);

    // // Reconstruct saddle point according to geo and relax system to saddle,
    // // check we have not constructed a false SP, if we have mark mechanism as poisoned
    // std::optional<system::SoA<Position>> recon_sp_relax(env::Geometry<Index> const& geo,
    //                                                     env::Mechanism& m,
    //                                                     system::SoA<Position const&, TypeID const&, Frozen const&> in);
    // // Check a relaxed minima is close to a true minima.
    // bool recon_min_relax(env::Geometry<Index> const& geo,
    //                      env::Mechanism const& m,
    //                      system::SoA<Position const&, TypeID const&, Frozen const&> in);

    ThreadData& thread() { return m_data[safe_cast<std::size_t>(omp_get_thread_num())]; }

  };  // namespace fly::saddle

}  // namespace fly::saddle