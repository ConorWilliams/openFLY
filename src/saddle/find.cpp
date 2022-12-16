// Copyright © 2020-2022 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see https :
// // www.gnu.org/licenses/>.

#include "libfly/saddle/find.hpp"

#include <fmt/core.h>
#include <fmt/format.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/geometry.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

/**
 * \file find.hpp
 *
 * @brief Class to coordinate a group of saddle-point searches.
 */

namespace fly::saddle {

  /////////// static functions ///////////

  // Reconstruct sp/minima according to reference geometry
  static system::SoA<Position> reconstruct(env::Geometry<Index> const& ref,
                                           system::VoS<Delta> const& m,
                                           system::SoA<Position const&> in) {
    system::SoA<Position> sp(in);

    for (int j = 0; j < m.size(); ++j) {
      sp(r_, ref[j][i_]) += m[j][del_];
    }

    return sp;
  }

  /**
   * @brief Reconstruct and relax a mechanism's saddle point and minima.
   */
  Master::Recon Master::recon_relax(env::Geometry<Index> const& geo,
                                    env::Mechanism const& m,
                                    system::SoA<Position const&, TypeID const&, Frozen const&> in) {
    //
    struct Recon res {
      reconstruct(geo, m.delta_fwd, in), reconstruct(geo, m.delta_sp, in), {}, {},
    };

    auto& thr = thread();

    {  // Minima
      system::SoA<Position, PotentialGradient, TypeID const&, Frozen const&> minima(in.size());
      minima.rebind(id_, in);
      minima.rebind(fzn_, in);
      minima[r_] = res.min[r_];

      if (!thr.min.minimise(minima, minima, thr.pot, 1)) {
        res.rel_min = std::move(minima);
      }
    }

    {  // SP
      system::SoA<Position, Axis, TypeID const&, Frozen const&> dim(in.size());
      dim.rebind(id_, in);
      dim.rebind(fzn_, in);
      dim[r_] = res.sp[r_];

      std::normal_distribution<double> dist;
      // Axis random and prop_to displacement.
      for (int j = 0; j < in.size(); ++j) {
        if (!in(fzn_, j)) {
          dim(ax_, j)
              = Vec::NullaryExpr([&] { return dist(thr.prng); }) * gnorm_sq(res.sp(r_, j) - in(r_, j));
        }
      }
      // ... and randomise.
      dim[ax_] /= gnorm(dim[ax_]);

      if (!thr.dimer.find_sp(dim, dim, in, thr.pot, {}, 1)) {
        res.rel_sp = std::move(dim);
      }
    }

    return res;
  }

  // Compute angle between VoSs (in degrees)
  static double theta_mech(system::VoS<Delta> const& a, system::VoS<Delta> const& b) {
    //
    verify(a.size() == b.size(), "");

    auto n = a.size();

    double ct = 0;
    double sa = 0;
    double sb = 0;

    for (int i = 0; i < n; i++) {
      ct += gdot(a[i][del_], b[i][del_]);

      sa += gnorm_sq(a[i][del_]);
      sb += gnorm_sq(b[i][del_]);
    }

    return std::acos(ct / std::sqrt(sa * sb)) / 2 / M_PI * 360.;
  }

  // Find the indices of minimally and maximally separated atoms.
  static std::pair<int, int> min_max(system::SoA<Position const&> a, system::SoA<Position const&> b) {
    //
    ASSERT(a.size() > 0, "Input has only {} atoms?", a.size());
    ASSERT(a.size() == b.size(), "min_max()'s inputs are different lengths {} != {}", a.size(), b.size());

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

  ////////////////////////////////////////

  Master::Master(Master::Options const& opt,
                 system::Box const& box,
                 potential::Generic const& pot,
                 minimise::LBFGS const& min,
                 saddle::Dimer const& dimer)
      : m_opt{opt}, m_box{box} {
    //
    std::random_device rd;

    Xoshiro prng(rd);

    for (int i = 0; i < opt.num_threads; i++) {
      m_data.push_back({pot, min, dimer, prng, {box, pot.r_cut()}, {}});
      prng.long_jump();
    }
  }

  std::vector<Master::LocalisedGeo> Master::package(std::vector<int> const& ix,
                                                    env::Catalogue const& cat,
                                                    int num_threads) {
    //
    std::vector<LocalisedGeo> out_data(ix.size());

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t i = 0; i < ix.size(); i++) {
      //
      out_data[i].geo = cat.get_geo(ix[i]);
      out_data[i].centre = out_data[i].geo[0][i_];
      out_data[i].syms = cat.calc_self_syms(ix[i]);
    }

    return out_data;
  }

  std::vector<Master::Found> Master::find_mechs(std::vector<LocalisedGeo> const& geo_data, SoA in) {
    // Early exit!
    if (geo_data.empty()) {
      return {};
    }

    neigh::List nl_pert{m_box, m_opt.r_pert};

    nl_pert.rebuild(in, m_opt.num_threads);

    calc_minima_hess(in);

    {
      m_sep_list.resize(safe_cast<std::size_t>(in.size()));

      auto mi = m_box.slow_min_image_computer();

#pragma omp parallel for num_threads(m_opt.num_threads)
      for (Eigen::Index i = 0; i < in.size(); i++) {
        //

        Eigen::Index k = 0;
        double max = 0;

        for (Eigen::Index j = 0; j < in.size(); j++) {
          if (double dr = mi(in(r_, i), in(r_, j)); dr > max) {
            k = j;
            max = dr;
          }
        }

        m_sep_list[safe_cast<std::size_t>(i)] = k;
      }
    }

    std::vector<Found> out(geo_data.size());

//
#pragma omp parallel
#pragma omp single nowait
    {
      for (std::size_t j = 0; j < geo_data.size(); ++j) {
#pragma omp task untied default(none) firstprivate(j) shared(out, in, nl_pert, geo_data)
        {
          find_n(out[j], geo_data[j], in, nl_pert);
          fmt::print("FINDER: Done @{:<4} found {} mechs\n",
                     geo_data[j].centre,
                     out[j] ? int(out[j].mechs().size()) : -1);
        }
      }
    }

    // Compute pre factors
    for (auto& f : out) {
      if (f) {
        for (auto& m : f.m_mechs) {
          m.kinetic_pre
              = std::sqrt(std::exp(m_log_prod_eigen - m.kinetic_pre) / (2 * M_PI * 1.6605390666050e-27));
        }
      }
    }

    if (m_opt.fout || m_opt.debug) {
      auto c = m_opt.fout->n_frames();

      for (std::size_t i = 0; i < out.size(); ++i) {
        if (out[i]) {
          for (auto& mech : out[i].m_mechs) {
            if (m_opt.fout) {
              m_opt.fout->commit(
                  [&] { m_opt.fout->write(r_, reconstruct(geo_data[i].geo, mech.delta_sp, in)); });
              m_opt.fout->commit(
                  [&] { m_opt.fout->write(r_, reconstruct(geo_data[i].geo, mech.delta_fwd, in)); });
            }
            fmt::print("FINDER: @{} frame={}, Delta={}eV, A={:e}Hz\n",
                       geo_data[i].centre,
                       c,
                       mech.barrier,
                       mech.kinetic_pre);

            c += 2;
          }
        }
      }
    }

    return out;
  }

  ///////////////////////////////////////////////////////

  void Master::find_n(Found& out, LocalisedGeo const& geo_data, SoA in, neigh::List const& nl_pert) {
    //

    auto N = safe_cast<std::size_t>(m_opt.batch_size);

    std::vector<Batch> batch(N, Batch{in.size()});

    std::vector<system::SoA<Position>> cache;  // Cache saddle-points

    int tot = 0;
    int fail = 0;

    int mod = in(id_, geo_data.centre) == 1 ? 2 : 1;

    while (tot < m_opt.max_searches * mod && fail < m_opt.max_failed_searches * mod) {
      //
      if (find_batch(tot, out, batch, geo_data, in, nl_pert, cache)) {
        fail = 0;
      } else {
        fail += m_opt.batch_size;
      }

      if (out.m_fail) {
        return;
      }

      tot += m_opt.batch_size;

      if (m_opt.debug) {
        //
        double t_min = 360;

        for (auto const& a : out.m_mechs) {
          for (auto const& b : out.m_mechs) {
            if (&a != &b) {
              t_min = std::min(t_min, theta_mech(a.delta_sp, b.delta_sp));
            }
          }
        }

        fmt::print("FINDER: {} mech(s) @{}: fail={}, tot={}, t_min={}deg, n-cache={}\n",
                   out.m_mechs.size(),
                   geo_data.centre,
                   fail,
                   tot,
                   t_min,
                   cache.size());
      }
    }
  }

  void Master::dump_recon(system::SoA<Position const&> in,
                          Index::scalar_t centre,
                          Recon const& recon,
                          system::SoA<Position const&> dimer) const {
    //
    constexpr char const* fname = "crash.check.gsd";

    fly::io::BinaryFile file(fname, fly::io::create);

    fmt::print(stderr,
               "ERROR: First symmetry (identity) should be guaranteed to reconstruct at atom {}. Failure was "
               "written to \"{}\" in the current working directory whose frames are: initial configuration, "
               "dimer final configuration, reconstructed saddle-point, reconstructed final-minima.\n",
               centre,
               fname);

    file.commit([&] {
      file.write("particles/N", fly::safe_cast<std::uint32_t>(in.size()));
      file.write(r_, in);
    });

    file.commit([&] { file.write(r_, dimer); });

    file.commit([&] { file.write(r_, recon.sp); });

    file.commit([&] { file.write(r_, recon.min); });

    if (recon.rel_sp) {
      fmt::print(stderr, "ERROR: and relaxed saddle, \n");
      file.commit([&] { file.write(r_, *recon.rel_sp); });
    }

    if (recon.rel_min) {
      fmt::print(stderr, "ERROR: and relaxed min, \n");
      file.commit([&] { file.write(r_, *recon.rel_min); });
    }
  }

  void Master::check_mech(Found& out,
                          system::SoA<Position>& cache_slot,
                          system::SoA<Position const&> dimer,
                          env::Mechanism const& mech,
                          std::size_t sym_indx,
                          LocalisedGeo const& geo_data,
                          SoA in) {
    bool fail;

#pragma omp atomic read
    fail = out.m_fail;

    if (fail) {
      return;
    }

    //
    auto recon = recon_relax(geo_data.geo, mech, in);

    auto set_fail_flag = [&, sym_indx] {
      if (sym_indx == 0) {
#pragma omp critical
        dump_recon(in, geo_data.centre, recon, dimer);
        throw error("First symmetry (identity) should be guaranteed to reconstruct");
      }

#pragma omp atomic write
      out.m_fail = true;
    };

    if (!recon.rel_min || !recon.rel_sp) {
      dprint(m_opt.debug,
             "FINDER: @ symmetry #{} recon_relax()->[{},{}]\n",
             sym_indx,
             bool(recon.rel_min),
             bool(recon.rel_sp));

      set_fail_flag();
      return;
    }

    centroid_align(*recon.rel_min, recon.min);
    centroid_align(*recon.rel_sp, recon.sp);

    ThreadData& thr = thread();

    thr.pot_nl.rebuild(in);
    double E0 = thr.pot.energy(in, thr.pot_nl, 1);

    thr.pot_nl.rebuild(*recon.rel_min);
    double Ef = thr.pot.energy(in, thr.pot_nl, 1);

    thr.pot_nl.rebuild(*recon.rel_sp);
    double Esp = thr.pot.energy(in, thr.pot_nl, 1);

    double re_barrier = Esp - E0;
    double re_delta = Ef - E0;

    double d_barrier = std::abs(re_barrier - mech.barrier);
    double d_delta = std::abs(re_delta - mech.delta);

    double frac_barrier = d_barrier / std::abs(mech.barrier);
    double frac_delta = d_delta / std::abs(mech.delta);

    double err_fwd = gnorm((*recon.rel_min)[r_] - recon.min[r_]);
    double err_sp = gnorm((*recon.rel_sp)[r_] - recon.sp[r_]);

    double d_fwd = std::abs(mech.err_fwd - err_fwd);
    double d_sp = std::abs(mech.err_sp - err_sp);

    double frac_fwd = d_fwd / mech.err_fwd;
    double frac_sp = d_sp / mech.err_sp + err_sp;

    dprint(m_opt.debug,
           "FINDER: Recon @{:>4} Δ(ΔE#)={:.3f} [{:.4f}] Δ(ΔE)={:.3f} [{:.4f}], ΔR_min={:.3f} [{:.4f}], "
           "ΔR_sp={:.3f} [{:.4f}]\n",
           geo_data.centre,
           d_barrier,
           frac_barrier,
           d_delta,
           frac_delta,
           d_fwd,
           frac_fwd,
           d_sp,
           frac_sp);

    if (d_barrier > m_opt.recon_e_tol_abs && frac_barrier > m_opt.recon_e_tol_frac) {
      set_fail_flag();
    } else if (d_delta > m_opt.recon_e_tol_abs && frac_delta > m_opt.recon_e_tol_frac) {
      set_fail_flag();
    } else if (frac_fwd > m_opt.recon_norm_frac_tol && d_fwd > m_opt.recon_norm_abs_tol) {
      set_fail_flag();
    } else if (frac_sp > m_opt.recon_norm_frac_tol && d_sp > m_opt.recon_norm_abs_tol) {
      set_fail_flag();
    } else {
      centroid_align(*recon.rel_sp, in);
      cache_slot = std::move(*recon.rel_sp);
    }
  }

  bool Master::find_batch(int tot,
                          Found& out,
                          std::vector<Batch>& batch,
                          LocalisedGeo const& geo_data,
                          SoA in,
                          neigh::List const& nl_pert,
                          std::vector<system::SoA<Position>>& cache) {
    // Abort SPS if the cosine of the angle between the dimer and a known SP is greater than this.
    double theta_tol = ((30 - 5) * std::exp(-0.02 * tot) + 5) / 360. * 2. * M_PI;

    dprint(m_opt.debug, "FINDER: Theta tolerance = {}\n", theta_tol / 2. / M_PI * 360.);

    //  Do batch_size SP searches
    for (Batch& elem : batch) {
#pragma omp task untied default(none) firstprivate(theta_tol) \
    shared(elem, in, nl_pert, batch, cache, geo_data)
      {
        perturb(elem.dimer, in, geo_data.centre, nl_pert);
        elem.mech = find_one(in, elem.dimer, elem.exit, geo_data.geo, cache, theta_tol);
      }
    }
#pragma omp taskwait

    bool found_one_or_more = false;

    // Process batch
    for (auto&& elem : batch) {
      if (!elem.mech) {
        ASSERT(elem.exit != Dimer::Exit::uninit, "Uninitialised exit/return code exit={}!", elem.exit);

        if (elem.exit == Dimer::Exit::success) {
          dprint(m_opt.debug, "FINDER: Caching failure\n");
          cache.emplace_back(elem.dimer);
        } else {
          dprint(m_opt.debug, "FINDER: Not caching SPS collision/fail\n");
        }
        continue;
      }

      // Found a mechanism.

      if (!is_new_mech(*elem.mech, out.m_mechs)) {
        dprint(m_opt.debug, "FINDER: Duplicate mech, caching\n");
        cache.emplace_back(elem.dimer);
        continue;
      }

      // Found a new mechanism.

      if (elem.mech->poison) {
        dprint(m_opt.debug, "FINDER: caching poisoned\n");
        out.m_mechs.push_back(std::move(*elem.mech));
        cache.emplace_back(elem.dimer);
        continue;
      }

      // Found a new non-poisoned mechanism.

      std::size_t num_new;

      try {
        num_new = append_syms(geo_data.syms, *elem.mech, out.m_mechs);
      } catch (...) {
#pragma omp critical
        {
          fly::io::BinaryFile file("crash.gsd", fly::io::create);

          fmt::print(
              stderr,
              "Append symmetries threw @{}, this implies a symmetrical saddle-point but sym-breaking minima, "
              "perhaps the minimiser reached the wrong minima i.e. failed to follow the minimum mode or the "
              "sp was higher order. Failure was written to \"crash.gsd\" in the current working directory\n",
              geo_data.centre);

          file.commit([&] {
            file.write("particles/N", fly::safe_cast<std::uint32_t>(in.size()));
            file.write(r_, in);
          });
          file.commit([&] { file.write(r_, reconstruct(geo_data.geo, elem.mech->delta_sp, in)); });
          file.commit([&] { file.write(r_, elem.dimer); });
          file.commit([&] { file.write(r_, reconstruct(geo_data.geo, elem.mech->delta_fwd, in)); });
        }
        throw;
      }

      verify(num_new > 0, "Only found {} but should have had at least 1", num_new);

      std::size_t n = cache.size();

      cache.resize(n + num_new);

      for (std::size_t k = 0; k < num_new; k++) {
#pragma omp task untied default(none) firstprivate(n, k, num_new) shared(out, in, cache, geo_data, elem)
        check_mech(
            out, cache[n + k], elem.dimer, out.m_mechs[out.m_mechs.size() - num_new + k], k, geo_data, in);
      }
#pragma omp taskwait

      found_one_or_more = true;
    }

    return found_one_or_more;
  }

  std::optional<env::Mechanism> Master::find_one(
      system::SoA<Position const&, Frozen const&, TypeID const&> in,
      system::SoA<Position&, Axis&> dimer_in_out,
      Dimer::Exit& exit,
      env::Geometry<Index> const& geo,
      std::vector<system::SoA<Position>> const& hist_sp,
      double theta_tol) {
    // Saddle search

    system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer(in.size());

    dimer.rebind(r_, dimer_in_out);
    dimer.rebind(ax_, dimer_in_out);
    dimer.rebind(fzn_, in);
    dimer.rebind(id_, in);

    if ((exit = find_sp(dimer, in, hist_sp, theta_tol))) {
      return {};
    }

    std::optional path = do_adj_min(dimer, in, geo[0][i_]);

    if (!path) {
      return {};
    }

    auto& [rev, sp, fwd] = *path;

    //  CHECK pathway is min->sp->min

    double i2r = gnorm(rev[r_] - in[r_]);
    double r2f = gnorm(fwd[r_] - rev[r_]);
    double r2w = gnorm(rev[r_] - sp[r_]);
    double w2f = gnorm(fwd[r_] - sp[r_]);

    dprint(m_opt.debug, "FINDER: i2r={:.4f}, r2f={:.4f} r2sp={:.4f} sp2f={:.4f}\n", i2r, r2f, r2w, w2f);

    if (i2r > m_opt.basin_tol) {
      dprint(m_opt.debug, "FINDER: Mech starting {}A from initial is unlinked\n", i2r);
      return {};
    }
    if (r2f < m_opt.basin_tol) {
      dprint(m_opt.debug, "FINDER: Mech total displacement={} => converged back to initial\n", r2f);
      return {};
    }
    if (r2w < m_opt.stationary_tol) {
      dprint(m_opt.debug, "FINDER: Min->Sp displacement={}\n", r2w);
      return {};
    }
    if (w2f < m_opt.stationary_tol) {
      dprint(m_opt.debug, "FINDER: Sp->Min displacement={}\n", w2f);
      return {};
    }

    // We now have a mech via min->sp->min, time to build a mechanism.

    ThreadData& thr = thread();

    thr.pot_nl.rebuild(rev);
    double E0 = thr.pot.energy(in, thr.pot_nl, 1);

    thr.pot_nl.rebuild(sp);
    double Esp = thr.pot.energy(in, thr.pot_nl, 1);

    thr.pot_nl.rebuild(fwd);
    double Ef = thr.pot.energy(in, thr.pot_nl, 1);

    env::Mechanism mech;

    mech.barrier = Esp - E0;
    mech.delta = Ef - E0;

    for (auto const& atom : geo) {
      //
      auto i = atom[i_];

      mech.delta_sp.emplace_back(sp(r_, i) - rev(r_, i));
      mech.delta_fwd.emplace_back(fwd(r_, i) - rev(r_, i));
    }

    // Test reconstruction

    system::SoA<Position const&, TypeID const&, Frozen const&> ref(rev.size());
    ref.rebind(r_, rev);
    ref.rebind(fzn_, in);
    ref.rebind(id_, in);

    auto [re_min, re_sp, rel_min, rel_sp] = recon_relax(geo, mech, ref);

    if (!rel_min || !rel_sp) {
      dprint(m_opt.debug, "FINDER: Mech poisoned!\n");
      mech.poison = true;
    } else {
      mech.err_fwd = gnorm((*rel_min)[r_] - re_min[r_]);
      mech.err_sp = gnorm((*rel_sp)[r_] - re_sp[r_]);

      dprint(m_opt.debug,
             "FINDER: In-place reconstruct, unaligned: err_fwd={}, err_sp={}\n",
             mech.err_fwd,
             mech.err_sp);

      centroid_align(*rel_min, fwd);
      centroid_align(*rel_sp, sp);

      thr.pot_nl.rebuild(*rel_min);
      double re_Ef = thr.pot.energy(in, thr.pot_nl, 1);

      thr.pot_nl.rebuild(*rel_sp);
      double re_Esp = thr.pot.energy(in, thr.pot_nl, 1);

      double err_fwd = gnorm((*rel_min)[r_] - fwd[r_]);
      double err_sp = gnorm((*rel_sp)[r_] - sp[r_]);

      double r_tol = m_opt.capture_r_tol;
      double e_tol = m_opt.capture_E_tol;

      double del_Esp = std::abs(re_Esp - Esp);
      double del_Efwd = std::abs(re_Ef - Ef);

      dprint(m_opt.debug,
             "FINDER: In-place reconstruct: err_fwd={} err_sp={} dEsp={} dEfwd={}\n",
             err_fwd,
             err_sp,
             del_Esp,
             del_Efwd);

      if (err_fwd > r_tol || err_sp > r_tol || del_Efwd > e_tol || del_Esp > e_tol) {
        dprint(m_opt.debug, "FINDER: Mech poisoned!\n");
        mech.poison = true;
      } else {
        mech.poison = false;
      }
    }

    if (mech.poison) {
      //   verify(mech.barrier > 2., "Mechanism @{} with energy barrier = {}eV is poisoned!", geo[0][i_],
      //   mech.barrier);
    }

    //////////////// Partial hessian compute. ////////////////

    thr.pot_nl.rebuild(sp, 1);

    thr.pot.hessian(thr.hess, in, thr.pot_nl, 1);

    thr.pot.mw_hessian(thr.hess, in, 1);

    system::Hessian::Vector const& freq = thr.hess.eigenvalues();

    // Must have at least one neg
    verify(freq[0] < -m_opt.hessian_eigen_zero_tol, "Saddle-point with minimum mode={}", freq[0]);

    int count_zeros = 0;
    mech.kinetic_pre = 0;

    for (int i = 1; i < freq.size(); i++) {
      //
      if (freq[i] < -m_opt.hessian_eigen_zero_tol) {
        dprint(m_opt.debug, "FINDER: Second order SP or higher, modes {} = {}\n", i, freq.head(i));
        return {};
      } else if (freq[i] > m_opt.hessian_eigen_zero_tol) {
        mech.kinetic_pre += std::log(freq[i]);
      } else {
        ++count_zeros;
      }
    }

    if (count_zeros != m_num_zero_modes) {
      for (int j = 0; j < std::max(m_num_zero_modes, count_zeros) + 3; j++) {
        fmt::print(stderr, "FINDER: [ERR] min mode {} = {}\n", j, freq[j]);
      }
      throw error("Expecting {} zero modes but hessian has {}", m_num_zero_modes, count_zeros);
    }

    if (m_opt.debug) {
      fmt::print(
          "FINDER: Mech ΔEsp={:.3e}, ΔEfwd={:.3e}, N_zeros={}\n", mech.barrier, mech.delta, count_zeros);

      for (int i = 0; i < m_num_zero_modes + 3; i++) {
        fmt::print("FINDER: sp mode {} = {}\n", i, freq[i]);
      }
    }

    return mech;
  }

  Dimer::Exit Master::find_sp(system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer,
                              system::SoA<Position const&> in,
                              std::vector<system::SoA<Position>> const& hist_sp,
                              double theta_tol) {
    //
    ThreadData& thr = thread();

    auto err = thr.dimer.find_sp(dimer, dimer, in, thr.pot, hist_sp, theta_tol, 1);

    if (err && m_opt.debug) {
      switch (err) {
        case Dimer::collision:
          fmt::print("FINDER: SPS fail - collision\n");
          break;
        case Dimer::iter_max:
          fmt::print("FINDER: SPS fail - iter_max\n");
          break;
        case Dimer::convex:
          fmt::print("FINDER: SPS fail - convex\n");
          break;
        default:
          ASSERT(false, "Unknown error = {}", err);
      }
    }

    return err;
  }

  /**
   * @brief Given a saddle point produce a min->sp->min pathway
   */
  std::optional<Master::Pathway> Master::do_adj_min(
      system::SoA<Position const&, Axis const&, Frozen const&, TypeID const&> dimer,
      system::SoA<Position const&> in,
      Index::scalar_t centre) {
    // Check sp centred on centre and freeze an atom.

    auto [_, max] = min_max(dimer, in);

    if (max != centre) {
      dprint(m_opt.debug, "FINDER: found mech at {} but wanted {}\n", max, centre);
      return {};
    }

    //   Minimisations

    double disp = gnorm(dimer[r_] - in[r_]);

    system::SoA<Position, PotentialGradient, Frozen, TypeID const&> relax{dimer.size()};
    relax[r_] = dimer[r_] + dimer[ax_] * disp * m_opt.nudge_frac;
    relax[fzn_] = dimer[fzn_];
    relax.rebind(id_, dimer);

    if (m_count_frozen == 0) {
      // Freeze furthest atom to remove translational degrees of freedom.

      auto min = m_sep_list[safe_cast<size_t>(max)];

      dprint(m_opt.debug, "FINDER: Freezing atom #{}\n", min);

      relax(fzn_, min) = true;
    }

    ThreadData& thr = thread();

    if (thr.min.minimise(relax, relax, thr.pot, 1)) {
      dprint(m_opt.debug, "FINDER: minimisation failed\n");
      return {};
    }

    system::SoA<Position> fwd{relax};

    relax[r_] = dimer[r_] - dimer[ax_] * disp * m_opt.nudge_frac;

    if (thr.min.minimise(relax, relax, thr.pot, 1)) {
      dprint(m_opt.debug, "FINDER: minimisation failed\n");
      return {};
    }

    system::SoA<Position> rev{relax};

    // Swap if forward went to final.
    if (gnorm(rev[r_] - in[r_]) > gnorm(fwd[r_] - in[r_])) {
      using std::swap;
      swap(fwd, rev);
    }

    // We now have a min->sp->min pathway with no translation (due to freeze). Need to correct for COM drift
    // during SPS.

    Vec in_com = centroid(in);
    Vec rev_com = centroid(rev);

    Vec drift = rev_com - in_com;

    dprint(m_opt.debug, "FINDER: Centre of mass drift={}\n", gnorm(drift));

    system::SoA<Position> sp_{dimer};

    for (int i = 0; i < in.size(); i++) {
      fwd(r_, i) -= drift;
      sp_(r_, i) -= drift;
      rev(r_, i) -= drift;
    }

    ASSERT(gnorm(centroid(rev) - centroid(in)) < 1e-8, "Drift correction failed", 0);

    return Pathway{std::move(rev), std::move(sp_), std::move(fwd)};
  }

  ///////////////////////////////////////////////////////////

  void Master::calc_minima_hess(system::SoA<Position const&, Frozen const&, TypeID const&> in) {
    //

    m_count_frozen = 0;

    for (int i = 0; i < in.size(); i++) {
      if (in(fzn_, i)) {
        m_count_frozen++;
      }
    }
    // If all unfrozen then zero modes due to translational degrees of freedom.
    int exp_zero_modes = m_count_frozen > 0 ? spatial_dims * m_count_frozen : spatial_dims;

    m_data[0].pot_nl.rebuild(in);

    m_data[0].pot.hessian(m_data[0].hess, in, m_data[0].pot_nl, m_opt.num_threads);

    m_data[0].pot.mw_hessian(m_data[0].hess, in, m_opt.num_threads);

    system::Hessian::Vector const& freq = m_data[0].hess.eigenvalues();

    m_log_prod_eigen = 0;
    m_num_zero_modes = 0;

    for (int i = 0; i < freq.size(); i++) {
      // Negative modes = not a minima.
      verify(freq[i] > -m_opt.hessian_eigen_zero_tol, "Mode {} has eigenvalue = {}", i, freq[i]);

      // Sum non-zero modes.
      if (freq[i] > m_opt.hessian_eigen_zero_tol) {
        m_log_prod_eigen += std::log(freq[i]);
      } else {
        ++m_num_zero_modes;
      }
    }

    if (m_num_zero_modes != exp_zero_modes) {
      for (int j = 0; j < std::max(m_num_zero_modes, exp_zero_modes) + 5; j++) {
        fmt::print(stderr, "FINDER: [ERR ]min mode {} = {}\n", j, freq[j]);
      }
      throw error("Expecting {} zero modes but fount {}", exp_zero_modes, m_num_zero_modes);
    }

    if (m_opt.debug) {
      fmt::print(
          "FINDER: Null space of Hessian has dimension  = {}, exp = {}\n", m_num_zero_modes, exp_zero_modes);

      for (int i = 0; i < m_num_zero_modes + 5; i++) {
        fmt::print("FINDER: Initial min mode {} = {}\n", i, freq[i]);
      }
    }
  }

  // Perturb in-place positions around centre and write axis,
  auto Master::perturb(system::SoA<Position&, Axis&> out,
                       system::SoA<Position const&, Frozen const&> in,
                       Index::scalar_t centre,
                       neigh::List const& nl) -> double {
    //
    Xoshiro& prng = thread().prng;

    std::normal_distribution normal(0., 1.);

    std::normal_distribution prime(m_opt.stddev, m_opt.stddev / 3.0);

    double stddev = -1;

    while (stddev < 0) {
      stddev = prime(prng);
    }

    dprint(m_opt.debug, "FINDER: standard deviation={}\n", stddev);

    std::normal_distribution<double> gauss(0, stddev);

    out[r_] = in[r_];
    out[ax_] = 0;

    nl.for_neighbours(centre, m_opt.r_pert, [&](auto n, double r, auto const&) {
      if (!in(Frozen{}, n)) {
        out(r_, n) += Vec::NullaryExpr([&] { return gauss(prng); }) * (1. - r / m_opt.r_pert);
        out(ax_, n) += Vec::NullaryExpr([&] { return normal(prng); });
      }
    });

    /* Experimental quasi-random centre

    constexpr double g = 1.32471795724474602596;  // Golden ratio
    constexpr double a1 = 1.0 / g;
    constexpr double a2 = 1.0 / (g * g);

    static int n = 0;

    double u = 0.5 + a1 * n - std::floor(0.5 + a1 * n);
    double v = 0.5 + a2 * n - std::floor(0.5 + a2 * n);

    n++;

    ASSERT(u >= 0 && u < 1, "u={}", u);
    ASSERT(u >= 0 && u < 1, "v={}", v);

    //    Lambert projection

    double phi = 2 * M_PI * v;
    double sin_lam = 2 * u - 1;
    double cos_lam = std::cos(std::asin(sin_lam));

    Vec p = {
        cos_lam * std::cos(phi),
        cos_lam * std::sin(phi),
        sin_lam,
    };

    */

    ASSERT(!in(Frozen{}, centre), "perturbation centred on a frozen atom {}", centre);

    out(r_, centre)
        += Vec::NullaryExpr([&] { return gauss(prng); });  // p * std::abs(gauss(prng)) * std::sqrt(3);
    out(ax_, centre) += Vec::NullaryExpr([&] { return normal(prng); });

    out[ax_] /= gnorm(out[ax_]);  // normalize

    return stddev;
  }

  bool Master::is_new_mech(env::Mechanism const& maybe, std::vector<env::Mechanism> const& hist) const {
    // All searches have been around same atom hence orientation is fixed.

    if (hist.empty()) {
      dprint(m_opt.debug, "FINDER: First mech must be new\n");
      return true;
    }

    double min_d_sp = std::numeric_limits<double>::max();
    double min_d_fwd = std::numeric_limits<double>::max();

    for (env::Mechanism const& m : hist) {
      //

      double d_sp = env::rmsd<Delta>(maybe.delta_sp, m.delta_sp);
      double d_fwd = env::rmsd<Delta>(maybe.delta_fwd, m.delta_fwd);

      if (d_sp < m_opt.mech_tol && d_fwd < m_opt.mech_tol) {
        dprint(m_opt.debug, "FINDER: Mech is NOT new, distance: old-sp={:.5f} old-min={:.5f}\n", d_sp, d_fwd);
        return false;
      }

      min_d_sp = std::min(min_d_sp, d_sp);
      min_d_fwd = std::min(min_d_fwd, d_fwd);
    }

    dprint(m_opt.debug,
           "FINDER: Mech is new, min distance: old-sp={:.5f} old-min={:.5f}\n",
           min_d_sp,
           min_d_fwd);
    return true;
  }

  std::size_t Master::append_syms(std::vector<env::Catalogue::SelfSymetry> const& syms,
                                  env::Mechanism const& new_mech,
                                  std::vector<env::Mechanism>& mechs) const {
    //
    env::Mechanism mech = new_mech;

    std::size_t count = 0;
    //
    for (auto const& [O, perm] : syms) {
      for (std::size_t i = 0; i < perm.size(); ++i) {
        mech.delta_sp[Eigen::Index(i)][del_].noalias() = O * new_mech.delta_sp[perm[i]][del_];
        mech.delta_fwd[Eigen::Index(i)][del_].noalias() = O * new_mech.delta_fwd[perm[i]][del_];
      }

      if (is_new_mech(mech, mechs)) {
        mechs.push_back(mech);
        ++count;
      }
    }

    dprint(m_opt.debug, "FINDER: Added {} mechs related by symmetry\n", count);

    return count;
  }

}  // namespace fly::saddle