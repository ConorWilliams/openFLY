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
      int num_threads = 1;
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
    MasterFinder(Options const& opt, potential::Generic const& pot, minimise::LBFGS const& min, saddle::Dimer const& dimer)
        : m_opt{opt} {
      //
      std::random_device rd;

      Xoshiro prng({rd(), rd(), rd(), rd()});

      for (int i = 0; i < opt.num_threads; i++) {
        m_data.push_back({pot, min, dimer, prng});
        prng.long_jump();
      }
    }

    /**
     * @brief Find all the mechanisms centred on the ``unknown`` atoms.
     *
     * @param box Description of the simulation space.
     * @param unknown List of indices of the atoms which the mechanisms must be centred on.
     * @param in Description of system to search in.
     */
    void find_all(system::Box const& box,
                  std::vector<int> const& unknown,
                  system::SoA<Position const&, Frozen const&, TypeID const&> in) {
      //
      neigh::List nl{box, m_opt.r_pert};

      nl.rebuild(in, m_opt.num_threads);

      std::vector<std::vector<Pathway>> mechs(unknown.size());

      //
#pragma omp parallel
#pragma omp single nowait
      {
        for (std::size_t j = 0; j < mechs.size(); ++j) {
            //
            int index = unknown[j];
          
#pragma omp task untied default(none) firstprivate(index, j) shared(in, nl, mechs)
          {
            std::vector<std::optional<Pathway>> batch(safe_cast<std::size_t>(m_opt.batch_size));

            int tot = 0;
            int fail = 0;

            while (tot < m_opt.max_searches && fail < m_opt.max_failed_searches) {
              //  Do batch_size SP searches
              for (std::size_t i = 0; i < batch.size(); i++) {
#pragma omp task untied default(none) firstprivate(index, i) shared(in, nl, batch)
                {
          
                  batch[i] = find_one(in, nl, index);
                
                }
              }

#pragma omp taskwait
              // Process batch
              for (std::optional<Pathway>& elem : batch) {
                if (elem && push_if_new( mechs[j], std::move(*elem))) {
                  fail = 0;
                } else {
                  fail++;
                }
              }

              tot += m_opt.batch_size;

              if(m_opt.debug){
                fmt::print("FINDER: Found {} mech(s) at {}: fail={}, tot={}\n",  mechs[j].size(), index, fail, tot);
              }
            }
          }
        }
      }

    if(m_opt.fout){
        for(auto&& m : mechs){
            for(auto && v : m){
                m_opt.fout->commit([&]{
                    m_opt.fout->write(r_, v.rev);
                });
                m_opt.fout->commit([&]{
                    m_opt.fout->write(r_, v.sp);
                });
                  m_opt.fout->commit([&]{
                    m_opt.fout->write(r_, v.fwd);
                });

                if(m_opt.debug){
                    neigh::List nl2{box, thread().pot.r_cut()};

                    nl2.rebuild(v.rev,  m_opt.num_threads);

                    double E0 = thread().pot.energy(in, nl2, m_opt.num_threads); 

                    nl2.rebuild(v.sp,  m_opt.num_threads);

                    double Esp = thread().pot.energy(in, nl2, m_opt.num_threads); 

                    fmt::print("Delta={}eV\n", Esp - E0);
                }
            }
        }
    }

    }  

  private:
    struct Pathway {
      system::SoA<Position> rev;
      system::SoA<Position> sp;
      system::SoA<Position> fwd;

      template <typename T, typename U, typename V>
      Pathway(T&& t, U&& u, V&& v) : rev(std::forward<T>(t)), sp(std::forward<U>(u)), fwd(std::forward<V>(v)) {}
    };

    struct ThreadData {
      potential::Generic pot;
      minimise::LBFGS min;
      saddle::Dimer dimer;
      Xoshiro prng;
    };

    std::vector<ThreadData> m_data;

    Options m_opt;

    // Returns true if managed to push.
     bool push_if_new(std::vector<Pathway> & found, Pathway&& maybe) const noexcept {
        for(Pathway const & elem : found){

            double dsp = gnorm(elem.sp[r_] - maybe.sp[r_]); // SP displacement vector.
            double dfwd = gnorm(elem.fwd[r_] - maybe.fwd[r_]); // FWD displacement vector.

            if (m_opt.debug)
            {
                fmt::print("FINDER: dsp={} dfwd={}\n", dsp, dfwd);
            }

            // fmt::print("Range: {}\n", (elem.sp[r_] - maybe.sp[r_]).head(25));
            Eigen::Index maxRow;
            auto mm = (elem.sp[r_] - maybe.sp[r_]).abs().maxCoeff(&maxRow);
            fmt::print("I max = {} @ {}, {}\n", mm, maxRow, maxRow /3);

            if(dsp < m_opt.basin_tol && dfwd < m_opt.basin_tol){
                if (m_opt.debug)
                {
                    fmt::print("FINDER: Found same SP as before\n");
                }
                return false;
            }
        }

        found.push_back(std::move(maybe));

        return true;
    }

    // Find the indices of minimally and maximally separated atoms.
    static std::pair<int, int> min_max(system::SoA<Position const&> a, system::SoA<Position const&> b) {
      //
      ASSERT(a.size() > 0, "Input has only {} atoms?", a.size());
      ASSERT(a.size() == b.size(), "min_max inputs are different lengths {} != {}", a.size(), b.size());

      int min = 0;
      int max = 0;

      double r_min = std::numeric_limits<double>::max();
      double r_max = 0;

      for (int i = 1; i < a.size(); i++) {
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



    std::optional<Pathway> find_one(system::SoA<Position const&, Frozen const&, TypeID const&> in, neigh::List const& nl, int index) {
      //
      system::SoA<Position, Axis, Frozen const&, TypeID const&> walker(in.size());

      walker[r_] = in[r_];
      walker.rebind(fzn_, in);
      walker.rebind(id_, in);

      perturb(walker, index, nl);

      ThreadData& thr = thread();

      if (thr.dimer.step(walker, walker, thr.pot, 1000, 1)) {
        return {};
      }

      //   Minimisations

      auto [min, max] = min_max(walker, in);

      if (max != index) {
        if (m_opt.debug) {
          fmt::print("FINDER: found mech at {} but wanted {}\n", max, index);
        }
        return {};
      }

      if (double dr = gnorm(in(r_, min) - walker(r_, min)); dr > m_opt.freeze_tol) {
        throw error("FINDER: Trying to freeze an atom {} that displaced {}", min, dr);
      } else if (m_opt.debug) {
        fmt::print("FINDER: Freezing atom #{}\n", min);
      }

      double disp = gnorm(walker[r_] - in[r_]);

      system::SoA<Position, PotentialGradient, Frozen, TypeID const&> fwd{walker.size()};
      fwd[r_] = walker[r_] + walker[ax_] * disp * m_opt.nudge_frac;
      fwd[fzn_] = walker[fzn_];
      fwd(fzn_, min) = true;  // Freeze minimally displaced atom.
      fwd.rebind(id_, in);

      if (thr.min.minimise(fwd, fwd, thr.pot, 1)) {
        return {};
      }

      system::SoA<Position, PotentialGradient, Frozen, TypeID const&> rev = fwd;
      rev[r_] = walker[r_] - walker[ax_] * disp * m_opt.nudge_frac;

      if (thr.min.minimise(rev, rev, thr.pot, 1)) {
        return {};
      }

      // ///////////////

      double i2f = gnorm(fwd[r_] - in[r_]);
      double i2r = gnorm(rev[r_] - in[r_]);

      // Swap if forward went to final.
      if (i2r > i2f) {
        using std::swap;
        swap(fwd, rev);
        swap(i2f, i2r);
      }

      if (i2r > m_opt.basin_tol) {
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

      if (double r2w = gnorm(rev[r_] - walker[r_]); r2w < m_opt.stationary_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Min->Sp displacement={}\n", r2w);
        }
        return {};
      }

      if (double w2f = gnorm(fwd[r_] - walker[r_]); w2f < m_opt.stationary_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Sp->Min displacement={}\n", w2f);
        }
        return {};
      }

     

      walker[r_] = in[r_] + walker[r_] - rev[r_];
      fwd[r_] = in[r_] + fwd[r_] - rev[r_];

      //

      return Pathway(in, walker, fwd);
    }

    auto perturb(system::SoA<Position&, Axis&, Frozen const&> inout, int centre, neigh::List const& nl) -> void {
      //
      Xoshiro& prng = thread().prng;

      std::normal_distribution<double> normal(0, 1);

      std::normal_distribution<double> prime(m_opt.stddev, m_opt.stddev / 3.0);

      double stddev = std::abs(prime(prng));

      if (m_opt.debug) {
        fmt::print("FINDER: stddev={}\n", stddev);
      }

      std::normal_distribution<double> gauss(0, stddev);

      inout[ax_] = 0;

      nl.for_neighbours(centre, m_opt.r_pert, [&](auto n, double r, auto const&) {
        if (!inout(Frozen{}, n)) {
          inout(r_, n) += (1. - r / m_opt.r_pert) * Vec{gauss(prng), gauss(prng), gauss(prng)};
          inout(ax_, n) += Vec{normal(prng), normal(prng), normal(prng)};
        }
      });

      ASSERT(!inout(Frozen{}, centre), "perturbation centred on a frozen atom {}", centre);

      inout(r_, centre) += Vec{gauss(prng), gauss(prng), gauss(prng)};
      inout(ax_, centre) += Vec{normal(prng), normal(prng), normal(prng)};

      inout[ax_] /= gnorm(inout[ax_]);  // normalize
    }

    MasterFinder(MasterFinder const&) { fmt::print("copy\n"); }

    ThreadData& thread() {
      //   fmt::print("thread {} accesses its data\n", omp_get_thread_num());
      return m_data[safe_cast<std::size_t>(omp_get_thread_num())];
    }

  };  // namespace fly::saddle

}  // namespace fly::saddle