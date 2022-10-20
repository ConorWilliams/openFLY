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

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <random>
#include <vector>

#include "libfly/env/geometry.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/mechanisms.hpp"
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
   * The ``Master`` coordinates saddle-point finding for a set of openMP threads, typically one ``Master`` should be
   * instantiated per compute node. Each ``Master`` stores all the threads reusable objects: minimiser, prng, etc.
   *
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
     * @brief Construct a new Master Finder object to manage ``opt.num_threads`` openMP threads.
     *
     * This will create ``opt.num_threads`` copies of each of the input parameters, one for each thread.
     *
     * @param opt The configuration options.
     * @param box ...
     * @param pot The potential to use for minimisations and SP searches.
     * @param min The minimiser to relax either side of a SP.
     * @param dimer Used to find SPs.
     */
    Master(Options const& opt,
           system::Box const& box,
           potential::Generic const& pot,
           minimise::LBFGS const& min,
           saddle::Dimer const& dimer)
        : m_opt{opt}, m_box{box} {
      //
      std::random_device rd;

      Xoshiro prng({rd(), rd(), rd(), rd()});

      for (int i = 0; i < opt.num_threads; i++) {
        m_data.push_back({pot, min, dimer, prng, {box, pot.r_cut()}});
        prng.long_jump();
      }
    }

    /**
     * @brief Find all the mechanisms centred on the ``unknown`` atoms.
     *
     * @param geos List of indices of the atoms which the mechanisms must be centred on.
     * @param in Description of system to search in.
     *
     */
    auto find_mechs(std::vector<env::Geometry<Index>> const& geos, system::SoA<Position const&, Frozen const&, TypeID const&> in) {
      // Step one get the self symmetries of each geometry.
      std::vector<Data> geo_data = process_geos(geos);

      neigh::List nl_pert{m_box, m_opt.r_pert};

      nl_pert.rebuild(in, m_opt.num_threads);

//
#pragma omp parallel
#pragma omp single nowait
      {
        for (std::size_t j = 0; j < geo_data.size(); ++j) {
#pragma omp task untied default(none) firstprivate(j) shared(in, nl_pert, geo_data, geos)
          {
            auto N = safe_cast<std::size_t>(m_opt.batch_size);

            std::vector<std::optional<system::LocalMech>> mechs(N);           // Batched

            int tot = 0;
            int fail = 0;

            while (tot < m_opt.max_searches && fail < m_opt.max_failed_searches) {
              //  Do batch_size SP searches
              for (std::size_t i = 0; i < N; i++) {
#pragma omp task untied default(none) firstprivate(i, j) shared(in, nl_pert, mechs, geos)
                {
                  system::SoA<Position, Axis> dimer(in.size());
                  perturb(dimer, in, geos[j][0][i_], nl_pert);
                  mechs[i] = find_one(in, dimer, geos[j]);
                }
              }

              // Wait for batch
#pragma omp taskwait
              // Process batch
              for (std::size_t i = 0; i < N; i++) {
                if (mechs[i]) {
                  // SPS succeed.
                  geo_data[j].mechs.push_back(std::move(*mechs[i]));
                  fail = 0;
                } else {
                  fail++;
                }
              }

              tot += m_opt.batch_size;

              if (m_opt.debug) {
                fmt::print("FINDER: Found {} mech(s) at {}: fail={}, tot={}\n", geo_data[j].mechs.size(), geo_data[j].centre, fail, tot);
              }
            }
          }
        }
      }


    //   Output
    if (m_opt.fout) {
      for (std::size_t i = 0; i < geo_data.size(); ++i) {
        for (auto&& mech : geo_data[i].mechs) {
          
          system::SoA<Position> outer(in);
          
          m_opt.fout->commit([&] { m_opt.fout->write(r_, outer); });

          for(int j = 0; j < mech.delta_sp.size(); ++j){
             outer(r_, geos[i][j][i_]) += mech.delta_sp[j][del_];
          }
        
          m_opt.fout->commit([&] { m_opt.fout->write(r_, outer); });

        outer[r_] = in[r_];

          for(int j = 0; j < mech.delta_fwd.size(); ++j){
             outer(r_, geos[i][j][i_]) += mech.delta_fwd[j][del_];
          }
    
          m_opt.fout->commit([&] { m_opt.fout->write(r_, outer); });

          fmt::print("Delta={}eV\n", mech.barrier);
        }
      }
    }

    }

  private:
    struct ThreadData {
      potential::Generic pot;
      minimise::LBFGS min;
      saddle::Dimer dimer;
      Xoshiro prng;
      neigh::List pot_nl;
    };

    std::vector<ThreadData> m_data;

    Options m_opt;
    system::Box m_box;

    ///

    struct Data {
      Index::scalar_t centre;                                        // Central atom
      std::vector<std::pair<Mat, std::vector<Index::scalar_t>>> tr;  // Perm + transformation.
      std::vector<system::LocalMech> mechs;                          // Unique, mechanisms
    };

    std::vector<Data> process_geos(std::vector<env::Geometry<Index>> const& geos) const {
      //
      std::vector<Data> out_data(size_t(geos.size()));

#pragma omp parallel for num_threads(m_opt.num_threads) schedule(static)
      for (std::size_t i = 0; i < geos.size(); i++) {
        //
        out_data[i].centre = geos[i][0][i_];

        double r_min = std::numeric_limits<double>::max();

        for (int j = 0; j < geos[i].size(); j++) {
          for (int k = 0; k < j; k++) {
            if (geos[i][j][col_] == geos[i][k][col_]) {
              r_min = std::min(r_min, gnorm(geos[i][j][r_] - geos[i][k][r_]));
            }
          }
        }

        env::Geometry copy = geos[i];

        double delta = r_min * 0.1;

        env::for_equiv_perms(copy, geos[i], delta, 1, [&](fly::Mat const& O, double) {
          //
          std::vector<Index::scalar_t> perm;

          for (auto const& elem : copy) {
            perm.push_back(elem[i_]);
          }

          out_data[i].tr.emplace_back(O, std::move(perm));

          return false;
        });

        if (m_opt.debug) {
          fmt::print("Env @{} has {} symmetries @tol={}\n", out_data[i].centre, out_data[i].tr.size(), delta);
        }
      }

      return out_data;
    }

    // Find the indices of minimally and maximally separated atoms.
    static  std::pair<int, int> min_max(system::SoA<Position const&> a, system::SoA<Position const&> b) {
    //
    ASSERT(a.size() > 0, "Input has only {} atoms?", a.size());
    ASSERT(a.size() == b.size(), "min_max inputs are different lengths {} != {}", a.size(), b.size());

    int min = 0;
    int max = 0;

    double r_min = std::numeric_limits<double>::max();
    double r_max = 0;

    for (int i = 0; i < a.size(); i++) {
      //
      double r = gnorm_sq(a(r_, i) - b(r_, i));

      if (r < r_min) {
        min = i;
        r_min = r;
      }
      if (r > r_max) {
        max = i;
        r_max = r;
      }
    }

    return {min, max};
  }

    bool find_sp(system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer) {
      //
      ThreadData& thr = thread();

      if (thr.dimer.step(dimer, dimer, thr.pot, 1000, 1)) {
        return false;
      };

      return true;
    }

    struct Pathway {
      system::SoA<Position> rev;
      system::SoA<Position> sp;
      system::SoA<Position> fwd;
    };

    

    /**
     * @brief Given a saddle point produce a min->sp->min pathway
     */
    std::optional<Pathway> find_path(system::SoA<Position const&, Axis const&, Frozen const&, TypeID const&> dimer,
                                     system::SoA<Position const&> in,
                                     Index::scalar_t centre) {
      // Check sp centred on centre and freeze an atom.

      auto [min, max] = min_max(dimer, in);

      if (max != centre) {
        if (m_opt.debug) {
          fmt::print("FINDER: found mech at {} but wanted {}\n", max, centre);
        }
        return {};
      }

      if (double dr = gnorm(in(r_, min) - dimer(r_, min)); dr > m_opt.freeze_tol) {
        throw error("FINDER: Trying to freeze an atom {} that displaced {}", min, dr);
      } else if (m_opt.debug) {
        fmt::print("FINDER: Freezing atom #{}\n", min);
      }

      //   Minimisations

      double disp = gnorm(dimer[r_] - in[r_]);

      system::SoA<Position, PotentialGradient, Frozen, TypeID const&> relax{dimer.size()};
      relax[r_] = dimer[r_] + dimer[ax_] * disp * m_opt.nudge_frac;
      relax[fzn_] = dimer[fzn_];
      relax(fzn_, min) = true;  // Freeze minimally displaced atom.
      relax.rebind(id_, dimer);

      ThreadData& thr = thread();

      if (thr.min.minimise(relax, relax, thr.pot, 1)) {
        if (m_opt.debug) {
          fmt::print("FINDER: minimisation failed\n");
        }
        return {};
      }

      system::SoA<Position> fwd{relax};

      relax[r_] = dimer[r_] - dimer[ax_] * disp * m_opt.nudge_frac;

      if (thr.min.minimise(relax, relax, thr.pot, 1)) {
        if (m_opt.debug) {
          fmt::print("FINDER: minimisation failed\n");
        }
        return {};
      }

      system::SoA<Position> rev{relax};

      // Swap if forward went to final.
      if (gnorm(rev[r_] - in[r_]) > gnorm(fwd[r_] - in[r_])) {
        using std::swap;
        swap(fwd, rev);
      }

      // We now have a min->sp->min pathway with no translation (due to freeze). Need to correct for COM drift during SPS.

      Vec in_com = com(in);
      Vec rev_com = com(rev);

      Vec drift = rev_com - in_com;

      if (m_opt.debug) {
        fmt::print("FINDER: Centre of mass drift={}\n", gnorm(drift));
      }

      system::SoA<Position> sp_{dimer};

      for (int i = 0; i < in.size(); i++) {
        fwd(r_, i) += drift;
        sp_(r_, i) += drift;
        rev(r_, i) += drift;
      }

      return Pathway{std::move(rev), std::move(sp_), std::move(fwd)};
    }

    /**
     * @brief
     *
     * @param in Initial basin
     * @param dimer_in_out Perturbed dimer as input and output.
     * @param geo Centred of perturbation.
     */
    std::optional<system::LocalMech> find_one(
                                              system::SoA<Position const&, Frozen const&, TypeID const&> in,
                                              system::SoA<Position&, Axis&> dimer_in_out,
                                              env::Geometry<Index> const& geo) {
      // Saddle search

      system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer(in.size());

      dimer.rebind(r_, dimer_in_out);
      dimer.rebind(ax_, dimer_in_out);
      dimer.rebind(fzn_, in);
      dimer.rebind(id_, in);

      if(!find_sp(dimer)){
        return {};
      }

      std::optional path = find_path(dimer, in, geo[0][i_]);

      if (!path) {
        return {};
      }

      auto& [rev, sp, fwd] = *path;

      //  CHECK pathway is min->sp->min

      if (double i2r = gnorm(rev[r_] - in[r_]); i2r > m_opt.basin_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Mech starting {}A from initial is unlinked\n", i2r);
        }
        return {};
      }
      if (double r2f = gnorm(fwd[r_] - rev[r_]); r2f < m_opt.basin_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Mech total displacement={} => converged back to initial\n", r2f);
        }
        return {};
      }
      if (double r2w = gnorm(rev[r_] - dimer[r_]); r2w < m_opt.stationary_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Min->Sp displacement={}\n", r2w);
        }
        return {};
      }
      if (double w2f = gnorm(fwd[r_] - dimer[r_]); w2f < m_opt.stationary_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Sp->Min displacement={}\n", w2f);
        }
        return {};
      }

      // We now have a mech via min->sp->min, time to build a global mech

      ThreadData& thr = thread();

      thr.pot_nl.rebuild(rev);
      double E0 = thr.pot.energy(in, thr.pot_nl, 1);

      thr.pot_nl.rebuild(sp);
      double Esp = thr.pot.energy(in, thr.pot_nl, 1);

      thr.pot_nl.rebuild(fwd);
      double Ef = thr.pot.energy(in, thr.pot_nl, 1);

      system::LocalMech mech;

      mech.barrier = Esp - E0;
      mech.delta = Ef - E0;
      mech.kinetic_pre = 5e14;

      for (auto const& atom : geo) {
        //
        auto i = atom[i_];

        mech.delta_sp.emplace_back(sp(r_, i) - rev(r_, i));
        mech.delta_fwd.emplace_back(fwd(r_, i) - rev(r_, i));
      }

      double sum_sq = 0;

      for (auto&& elem : mech.delta_fwd) {
        sum_sq += gnorm_sq(elem[del_]);
      }

      mech.capture_frac = std::sqrt(sum_sq) / gnorm(fwd[r_] - rev[r_]);

      return mech;
    }

    static Vec com(system::SoA<Position const&> p) {
      Vec vsum = Vec::Zero();
      for (int i = 0; i < p.size(); i++) {
        vsum += p(r_, i);
      }
      return vsum / p.size();
    }

    // Perturb in-place positions around centre and write axis,
    auto perturb(system::SoA<Position&, Axis&> out,
                 system::SoA<Position const&, Frozen const&> in,
                 Index::scalar_t centre,
                 neigh::List const& nl) -> void {
      //
      Xoshiro& prng = thread().prng;

      std::normal_distribution normal(0., 1.);

      std::normal_distribution prime(m_opt.stddev, m_opt.stddev / 3.0);

      double stddev = std::abs(prime(prng));

      if (m_opt.debug) {
        fmt::print("FINDER: stddev={}\n", stddev);
      }

      std::normal_distribution<double> gauss(0, stddev);

      out[r_] = in[r_];
      out[ax_] = 0;

      nl.for_neighbours(centre, m_opt.r_pert, [&](auto n, double r, auto const&) {
        if (!in(Frozen{}, n)) {
          out(r_, n) += (1. - r / m_opt.r_pert) * Vec{gauss(prng), gauss(prng), gauss(prng)};
          out(ax_, n) += Vec{normal(prng), normal(prng), normal(prng)};
        }
      });

      ASSERT(!in(Frozen{}, centre), "perturbation centred on a frozen atom {}", centre);

      out(r_, centre) += Vec{gauss(prng), gauss(prng), gauss(prng)};
      out(ax_, centre) += Vec{normal(prng), normal(prng), normal(prng)};

      out[ax_] /= gnorm(out[ax_]);  // normalize
    }

    ThreadData& thread() { return m_data[safe_cast<std::size_t>(omp_get_thread_num())]; }

  };  // namespace fly::saddle

}  // namespace fly::saddle