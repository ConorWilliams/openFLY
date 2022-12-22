#pragma once

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

#include <fstream>
#include <utility>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/kinetic/cache.hpp"
#include "libfly/kinetic/superbasin.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/saddle/find.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file skmc.hpp
 *
 * @brief
 */

namespace fly::kinetic {

  /**
   * @brief Updates/rebuilds a catalogue from a ``cell``.
   *
   * This function repeatedly calls ``.rebuild()`` on ``cat`` and coordinates saddle-point searches using
   * ``mast`` to find the mechanisms for the newly encountered environments.
   *
   * @param mast The saddle-point/mechanism finder.
   * @param cat The catalogue to update
   * @param cell The cell to rebuild the catalogue from.
   * @param num_threads The number of openMP threads to use.
   * @return true If new environments encountered.
   * @return false If no new environments encountered.
   */
  template <typename Map, typename... T>
  bool update_cat(saddle::Master& mast,
                  env::Catalogue& cat,
                  system::Supercell<Map, T...> const& cell,
                  int num_threads) {
    //
    std::vector<int> ix = cat.rebuild(cell, num_threads);

    if (ix.empty()) {
      return false;
    }

    int refines = 0;

    while (true) {
      //
      std::vector<int> fails;

      fmt::print("Update: New envs @{} with {} refines\n", ix, refines);

      std::vector found = mast.find_mechs(saddle::Master::package({ix}, cat), cell);

      /*
       * If find_mechs has failed, the failed environments must be too symmetric, we must refine them until
       * they are less symmetric.
       */

      for (std::size_t i = 0; i < found.size(); i++) {
        if (found[i]) {
          cat.set_mechs(ix[i], found[i].mechs());
        } else {
          fails.push_back(ix[i]);
        }
      }

      if (fails.empty()) {
        return true;
      } else {
        ++refines;
      }

      // Refine tolerance's, Refines only occur at new environments
      for (auto const& f : fails) {
        auto n = cat.calc_self_syms(f).size();
        verify(n > 1, "Env @{} cannot get any less symmetric!", f);
        do {
          double new_tol = cat.refine_tol(f, cat.get_ref(f).delta_max() / 1.5);
          verify(new_tol > 1e-6, "Probably converging to a true symmetry!");
        } while (n == cat.calc_self_syms(f).size());
      }

      std::vector tmp = cat.rebuild(cell, num_threads);

      // The atom whose symmetry tolerance increased still needs to be searched alongside any atoms that now
      // no longer match that environment.
      for (auto const& elem : tmp) {
        verify(std::find(fails.begin(), fails.end(), elem) == fails.end(), "Atom #{} found twice", elem);
      }

      // tmp now stores previous fails + new environments that don't match the refined tolerance.
      tmp.insert(tmp.end(), fails.begin(), fails.end());

      ix = std::move(tmp);
    }
  }

  /**
   * @brief Coordinate the running of a OLKMC simulation.
   */
  class SKMC {
  public:
    /**
     * @brief Configure the ``SKMC`` class.
     */
    struct Options {
      /**
       * @brief Configure debug printing.
       */
      bool debug = false;
      /**
       * @brief The name of the file to read the catalogue from.
       */
      std::string fread = "cat.bin";
      /**
       * @brief The name of the file to write the catalogue to.
       */
      std::string fwrite = fread;
      /**
       * @brief Options for the superbasin.
       */
      SuperCache::Options opt_cache = {};
      /**
       * @brief Options for the saddle-point finder.
       */
      saddle::Master::Options opt_master = {};
    };

    /**
     * @brief Construct a new SKMC object.
     */
    SKMC(Options const& opt,
         system::Box const& box,
         minimise::LBFGS const& min,
         potential::Generic const& pot,
         saddle::Dimer const& dimer)
        : m_opt(opt),
          m_mast{
              opt.opt_master,
              box,
              pot,
              min,
              dimer,
          },
          m_minimiser(min),
          m_pot(pot),
          m_cat([&] {
            if (std::ifstream fcat(m_opt.fread); fcat.good()) {
              dprint(m_opt.debug, "SKMC: Opening existing catalogue: \"{}\"\n", m_opt.fread);
              return env::Catalogue{{}, fcat};
            } else {
              dprint(m_opt.debug, "SKMC: Could not open catalogue: \"{}\"\n", m_opt.fread);
              return env::Catalogue{{}};
            }
          }()) {
      m_cat.optimize();
    }

    /**
     * @brief
     *
     * @param cell
     * @param num_threads
     * @param f
     */
    template <typename Map, typename... T, typename F>
    auto skmc(system::Supercell<Map, T...> const& cell, int num_threads, F const& f) -> void;

  private:
    Options m_opt;
    saddle::Master m_mast;
    minimise::LBFGS m_minimiser;
    potential::Generic m_pot;
    env::Catalogue m_cat;

    void dump_cat() const {
      std::ofstream file(m_opt.fwrite);
      m_cat.dump(file);
    };
  };

  template <typename Map, typename... T, typename F>
  auto SKMC::skmc(system::Supercell<Map, T...> const& in_cell, int num_threads, F const& f) -> void {
    //
    system::Supercell<Map, T...> cell = in_cell;

    neigh::List neigh_list(cell.box(), m_pot.r_cut());

    {  // Initial minimisation
      system::SoA<Position, PotentialGradient> out(cell.size());
      bool err = timeit("Minimise", [&] { return m_minimiser.minimise(out, cell, m_pot, num_threads); });
      cell[r_] = out[r_];
      verify(!err, "Minimiser failed");
    }

    auto energy = [&](system::SoA<Position const&> x) {
      neigh_list.rebuild(x, num_threads);
      return m_pot.energy(cell, neigh_list, num_threads);
    };

    if (kinetic::update_cat(m_mast, m_cat, cell, num_threads)) {
      dump_cat();
    }

    std::random_device rd;

    Xoshiro rng(rd);

    double time = 0;

    // For system that is reconstructed but not relaxed.
    system::SoA<Position, Frozen const&, TypeID const&> raw_recon(cell.size());
    raw_recon.rebind(fzn_, cell);
    raw_recon.rebind(id_, cell);
    // For relaxed raw_recon
    system::SoA<Position, PotentialGradient> rel_recon(cell.size());

    SuperCache super(m_opt.opt_cache, {m_opt.opt_cache.opt_sb, {m_opt.opt_cache.opt_basin, cell, m_cat}});

    for (int i = 0;; ++i) {
      //
      bool stop = false;

      timeit("SKMC-loop\n", [&] {
        ///////////// Select mechanism /////////////
        auto const& [m, atom, dt, basin, changed] = super.kmc_choice(rng);

        cell[r_] = super.state(basin)[r_];

        ///////////// Reconstruct mech /////////////

        verify(!m.poison, "KMC chose a poisoned mechanisms with dE={}", m.barrier);

        double E0 = energy(cell);  // Energy before mechanism

        m_cat.reconstruct(raw_recon, m, atom, cell, !changed, num_threads);

        auto err = m_minimiser.minimise(rel_recon, raw_recon, m_pot, num_threads);

        centroid_align(rel_recon, raw_recon);

        double Ef = energy(rel_recon);  // Energy after relax

        double dR_err = std::abs(gnorm(rel_recon[r_] - raw_recon[r_]) - m.err_fwd);
        double dE_err = std::abs(Ef - E0 - m.delta);

        double dR_err_frac = dR_err / m.err_fwd;
        double dE_err_frac = dE_err / std::abs(m.delta);

        dprint(m_opt.debug,
               "SKMC: dE_err={:.3f}[{:.3f}], dR_err={:.3f}[{:.3f}]\n",
               dE_err,
               dE_err_frac,
               dR_err,
               dR_err_frac);

        {
          bool fail = false;

          auto const& opt = m_mast.get_options();

          if (err) {
            fail = true;
          } else if (dE_err > opt.recon_e_tol_abs && dE_err_frac > opt.recon_e_tol_frac) {
            fail = true;
          } else if (dR_err > opt.recon_norm_frac_tol && dR_err_frac > opt.recon_norm_abs_tol) {
            fail = true;
          }

          if (fail) {
            dprint(m_opt.debug, "SKMC: Reconstruction failed @{}\n", atom);

            auto initial_assign = m_cat.get_ref(atom).cat_index();
            bool new_envs = false;

            do {
              double new_tol = m_cat.refine_tol(atom);
              dprint(m_opt.debug, "SKMC: Refined tolerance to {}\n", new_tol);
              new_envs = new_envs || kinetic::update_cat(m_mast, m_cat, cell, num_threads);
            } while (initial_assign == m_cat.get_ref(atom).cat_index());

            if (new_envs) {
              dump_cat();
            }

            super.reset({m_opt.opt_cache.opt_sb, {m_opt.opt_cache.opt_basin, cell, m_cat}});

            return;
          }
        }

        ///////////// Time /////////////

        time += dt;

        stop = std::invoke(f, time, std::as_const(cell), atom, m, system::SoA<Position const&>{rel_recon});

        cell[r_] = rel_recon[r_];

        ///////////// Update catalogue /////////////

        if (kinetic::update_cat(m_mast, m_cat, cell, num_threads)) {
          dump_cat();
        }

        ///////////// Update SuperBasin /////////////

        super.connect_from(basin, atom, m, cell, m_cat);

        dprint(m_opt.debug, "SKMC: Iteration #{} time = {:.3e}\n", i, time);
      });

      if (stop) {
        break;
      }
    }
  }

}  // namespace fly::kinetic