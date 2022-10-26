// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see https :  //
// www.gnu.org/licenses/>.

#include "libfly/saddle/find.hpp"

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

  /////////// static functions ///////////

  // Reconstruct saddle point according to reference geometry
  static system::SoA<Position> reconstruct_sp(env::Geometry<Index> const& ref,
                                              env::Mechanism const& m,
                                              system::SoA<Position const&> in) {
    system::SoA<Position> sp(in);

    for (int j = 0; j < m.delta_sp.size(); ++j) {
      sp(r_, ref[j][i_]) += m.delta_sp[j][del_];
    }

    return sp;
  }

  // Compute angle between mechanisms saddle points
  static double theta_mech(env::Mechanism const& a, env::Mechanism const& b) {
    //
    verify(a.delta_sp.size() == b.delta_sp.size(), "");

    auto n = a.delta_sp.size();

    double ct = 0;
    double sa = 0;
    double sb = 0;

    for (int i = 0; i < n; i++) {
      ct += gdot(a.delta_sp[i][del_], b.delta_sp[i][del_]);

      sa += gnorm_sq(a.delta_sp[i][del_]);
      sb += gnorm_sq(b.delta_sp[i][del_]);
    }

    return std::acos(ct / std::sqrt(sa * sb)) / 2 / M_PI * 360.;
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

  // Centre of mass of a SoA
  static Vec com(system::SoA<Position const&> p) {
    Vec vsum = Vec::Zero();
    for (int i = 0; i < p.size(); i++) {
      vsum += p(r_, i);
    }
    return vsum / p.size();
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

    Xoshiro prng({rd(), rd(), rd(), rd()});

    for (int i = 0; i < opt.num_threads; i++) {
      m_data.push_back({pot, min, dimer, prng, {box, pot.r_cut()}, {}});
      prng.long_jump();
    }
  }

  auto Master::find_mechs(std::vector<env::Geometry<Index>> const& geos, system::SoA<Position const&, Frozen const&, TypeID const&> in)
      -> std::vector<Master::Found> {
    //
    neigh::List nl_pert{m_box, m_opt.r_pert};

    nl_pert.rebuild(in, m_opt.num_threads);

    calc_minima_hess(in);

    // Step one get the self symmetries of each geometry.
    std::vector<Data> geo_data = process_geos(geos);

//
#pragma omp parallel
#pragma omp single nowait
    {
      for (std::size_t j = 0; j < geo_data.size(); ++j) {
#pragma omp task untied default(none) firstprivate(j) shared(in, nl_pert, geo_data, geos)
        find_n(geos[j], geo_data[j], in, nl_pert);
      }
    }

    std::vector<Found> out;

    for (auto& elem : geo_data) {
      if (!elem.mechs.empty()) {
        out.push_back({elem.centre, std::move(elem.mechs)});
      }
    }

    // Compute pre factors
    for (auto& f : out) {
      for (auto& m : f.mechs) {
        m.kinetic_pre = std::sqrt(std::exp(m_log_prod_eigen - m.kinetic_pre) / (2 * M_PI * 1.6605390666050e-27));
      }
    }

    if (m_opt.debug) {
      //
      int c = 1;

      for (auto& f : out) {
        //
        env::Geometry<Index> const* geo = nullptr;

        for (auto const& elem : geos) {
          if (elem[0][i_] == f.centre) {
            geo = &elem;
            break;
          }
        }

        ASSERT(geo, "Cannot locate geo @{}!", f.centre);

        for (auto& mech : f.mechs) {
          if (m_opt.fout) {
            m_opt.fout->commit([&] { m_opt.fout->write(r_, reconstruct_sp(*geo, mech, in)); });
          }
          fmt::print("FINDER: @{} frame={}, Delta={}eV, A={:e}Hz\n", f.centre, c++, mech.barrier, mech.kinetic_pre);
        }
      }
    }
    return out;
  }

  ///////////////////////////////////////////////////////

  void Master::find_n(env::Geometry<Index> const& geo,
                      Data& geo_data,
                      system::SoA<Position const&, Frozen const&, TypeID const&> in,
                      neigh::List const& nl_pert) {
    //

    auto N = safe_cast<std::size_t>(m_opt.batch_size);

    std::vector<Batch> batch(N, Batch{in.size()});

    std::vector<system::SoA<Position>> cache;  // Cache saddle-points

    int tot = 0;
    int fail = 0;

    while (tot < m_opt.max_searches && fail < m_opt.max_failed_searches) {
      //
      if (find_batch(tot, batch, geo, geo_data, in, nl_pert, cache)) {
        fail = 0;
      } else {
        fail += m_opt.batch_size;
      }

      tot += m_opt.batch_size;

      if (m_opt.debug) {
        //
        double t_min = 180;

        for (auto const& a : geo_data.mechs) {
          for (auto const& b : geo_data.mechs) {
            if (&a != &b) {
              t_min = std::min(t_min, theta_mech(a, b));
            }
          }
        }

        fmt::print("FINDER: {} mech(s) @{}: fail={}, tot={}, t_min={}, n-cache={}\n",
                   geo_data.mechs.size(),
                   geo_data.centre,
                   fail,
                   tot,
                   t_min,
                   cache.size());
      }
    }
  }

  bool Master::find_batch(int tot,
                          std::vector<Batch>& batch,
                          env::Geometry<Index> const& geo,
                          Data& geo_data,
                          system::SoA<Position const&, Frozen const&, TypeID const&> in,
                          neigh::List const& nl_pert,
                          std::vector<system::SoA<Position>>& cache) {
    // Abort SPS if the cosine of the angle between the dimer and a known SP is greater than this.
    double theta_tol = ((30 - 5) * std::exp(-0.02 * tot) + 5) / 360. * 2. * M_PI;

    if (m_opt.debug) {
      fmt::print("FINDER: Theta tolerance = {}\n", theta_tol / 2. / M_PI * 360.);
    }

    //  Do batch_size SP searches
    for (Batch& elem : batch) {
#pragma omp task untied default(none) firstprivate(theta_tol) shared(elem, in, nl_pert, batch, geo, cache, geo_data)
      {
        perturb(elem.dimer, in, geo_data.centre, nl_pert);
        elem.mech = find_one(in, elem.dimer, elem.collsion, geo, cache, theta_tol);
      }
    }
#pragma omp taskwait

    bool found_one_or_more = false;

    // Process batch
    for (auto&& elem : batch) {
      if (elem.mech && is_new_mech(*elem.mech, geo_data.mechs)) {
        // Found a new mechanism, append all symmetries and cache them all

        std::size_t num_new = append_syms(geo_data.tr, *elem.mech, geo_data.mechs);

        ASSERT(num_new > 0, "Only found {} but should have had at least 1", num_new);

        std::size_t n = cache.size();

        cache.resize(n + num_new);

        for (std::size_t k = 0; k < num_new; k++) {
#pragma omp task untied default(none) firstprivate(n, k, num_new) shared(in, geo, cache, geo_data)
          { cache[n + k] = (recon_sp_relax(geo, geo_data.mechs[geo_data.mechs.size() - num_new + k], in)); }
        }
#pragma omp taskwait

        found_one_or_more = true;

      } else if (elem.mech && m_opt.debug) {
        // Found a previously known mechanism

        fmt::print("FINDER: Duplicate mech\n");

      } else if (!elem.mech && elem.collsion == false) {
        // Failure by collision with known mechanisms

        if (m_opt.debug) {
          fmt::print("FINDER: Caching failure\n");
        }

        cache.emplace_back(elem.dimer);  // Add failures

      } else if (m_opt.debug) {
        //
        fmt::print("FINDER: Not caching SPS collision/fail\n");
      }
    }

    return found_one_or_more;
  }

  std::optional<env::Mechanism> Master::find_one(system::SoA<Position const&, Frozen const&, TypeID const&> in,
                                                 system::SoA<Position&, Axis&> dimer_in_out,
                                                 bool& collision,
                                                 env::Geometry<Index> const& geo,
                                                 std::vector<system::SoA<Position>> const& hist_sp,
                                                 double theta_tol) {
    // Saddle search

    system::SoA<Position&, Axis&, Frozen const&, TypeID const&> dimer(in.size());

    dimer.rebind(r_, dimer_in_out);
    dimer.rebind(ax_, dimer_in_out);
    dimer.rebind(fzn_, in);
    dimer.rebind(id_, in);

    if (find_sp(dimer, in, hist_sp, theta_tol)) {
      collision = true;
      return {};
    } else {
      collision = false;
    }

    std::optional path = do_adj_min(dimer, in, geo[0][i_]);

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
    if (double r2w = gnorm(rev[r_] - sp[r_]); r2w < m_opt.stationary_tol) {
      if (m_opt.debug) {
        fmt::print("FINDER: Min->Sp displacement={}\n", r2w);
      }
      return {};
    }
    if (double w2f = gnorm(fwd[r_] - sp[r_]); w2f < m_opt.stationary_tol) {
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

    env::Mechanism mech;

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

    double cap = std::sqrt(sum_sq);
    double tot = gnorm(fwd[r_] - rev[r_]);

    mech.capture_frac = cap / tot;

    //////////////// Partial hessian compute. ////////////////

    thr.pot_nl.rebuild(sp, 1);

    thr.pot.hessian(thr.hess, in, thr.pot_nl);

    system::Hessian::Vector const& freq = thr.hess.eigenvalues();

    ASSERT(freq.size() > 1 + m_deg_free, "Only {} eigenvalues", freq.size());

    // Must have at least one neg
    verify(freq[0] < -m_opt.hessian_eigen_zero_tol, "Saddle-point with minimum mode={}", freq[0]);

    for (int i = 1; i < 1 + m_deg_free; i++) {
      // More than one is a higher order.
      if (freq[i] < -m_opt.hessian_eigen_zero_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Second order SP or higher, modes = {}\n", freq.head(10));
        }
        return {};
      }
      verify(std::abs(freq[i]) < m_opt.hessian_eigen_zero_tol, "Expecting zero eigen-value, got modes = {}", freq.head(10));
    }
    verify(freq[1 + m_deg_free] > m_opt.hessian_eigen_zero_tol, "Seem to have too many zero modes = {}", freq.head(10));

    mech.kinetic_pre = 0;

    for (int i = 1 + m_deg_free; i < freq.size(); i++) {
      mech.kinetic_pre += std::log(freq[i]);
    }

    if (m_opt.debug) {
      fmt::print("FINDER: Mech ΔE={:.2e} f={:.4f}, uncap={:.2e}, modes = {::.2e}, \n",
                 mech.barrier,
                 mech.capture_frac,
                 tot - cap,
                 freq.head(5));
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
  std::optional<Master::Pathway> Master::do_adj_min(system::SoA<Position const&, Axis const&, Frozen const&, TypeID const&> dimer,
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

  ///////////////////////////////////////////////////////////

  void Master::calc_minima_hess(system::SoA<Position const&, Frozen const&, TypeID const&> in) {
    //
    m_deg_free = 0;

    for (int i = 0; i < spatial_dims; i++) {
      if (m_box.periodic(i)) {
        m_deg_free++;
      }
    }

    for (auto&& elem : in[fzn_]) {
      if (elem) {
        m_deg_free = 0;
        break;
      }
    }

    if (m_opt.debug) {
      fmt::print("FINDER: Degrees of freedom = {}\n", m_deg_free);
    }

    m_data[0].pot_nl.rebuild(in);

    m_data[0].pot.hessian(m_data[0].hess, in, m_data[0].pot_nl);

    system::Hessian::Vector const& freq = m_data[0].hess.eigenvalues();

    ASSERT(freq.size() > m_deg_free, "Only {} normal modes", freq.size());

    for (int i = 0; i < m_deg_free; i++) {
      verify(std::abs(freq[i]) < m_opt.hessian_eigen_zero_tol, "Master's input is not minima, modes (1)= {}", freq.head(10));
    }

    verify(freq[m_deg_free] > m_opt.hessian_eigen_zero_tol, "Master input modes (2)= {}", freq.head(10));

    double sum = 0;

    for (int i = m_deg_free; i < freq.size(); i++) {
      sum += std::log(freq[i]);
    }

    m_log_prod_eigen = sum;
  }

  std::vector<Master::Data> Master::process_geos(std::vector<env::Geometry<Index>> const& geos) const {
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

      for (int j = 0; j < copy.size(); j++) {
        copy[j][i_] = j;
      }

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
        fmt::print("FINDER: Env @{} has {} symmetries @tol={}\n", out_data[i].centre, out_data[i].tr.size(), delta);
      }
    }

    return out_data;
  }

  // Perturb in-place positions around centre and write axis,
  auto Master::perturb(system::SoA<Position&, Axis&> out,
                       system::SoA<Position const&, Frozen const&> in,
                       Index::scalar_t centre,
                       neigh::List const& nl) -> void {
    //
    Xoshiro& prng = thread().prng;

    std::normal_distribution normal(0., 1.);

    std::normal_distribution prime(m_opt.stddev, m_opt.stddev / 3.0);

    double stddev = -1;

    while (stddev < 0) {
      stddev = prime(prng);
    }

    if (m_opt.debug) {
      fmt::print("FINDER: standard deviation={}\n", stddev);
    }

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

    out(r_, centre) += Vec::NullaryExpr([&] { return gauss(prng); });  // p * std::abs(gauss(prng)) * std::sqrt(3);
    out(ax_, centre) += Vec::NullaryExpr([&] { return normal(prng); });

    out[ax_] /= gnorm(out[ax_]);  // normalize
  }

  bool Master::is_new_mech(env::Mechanism const& maybe, std::vector<env::Mechanism> const& hist) const {
    // All searches have been around same atom hence orientation is fixed.
    for (env::Mechanism const& m : hist) {
      //
      double d_sp = env::rmsd<Delta>(maybe.delta_sp, m.delta_sp);
      double d_fwd = env::rmsd<Delta>(maybe.delta_fwd, m.delta_fwd);

      if (d_sp < m_opt.mech_tol) {
        verify(d_fwd < m_opt.mech_tol, "Same sp, different final, d_fwd={}, d_sp={}", d_fwd, d_sp);
        return false;
      }
    }
    return true;
  }

  std::size_t Master::append_syms(std::vector<std::pair<Mat, std::vector<Index::scalar_t>>> const& syms,
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

    if (m_opt.debug) {
      fmt::print("FINDER: Added {} mechs related by symmetry\n", count);
    }

    return count;
  }

  system::SoA<Position> Master::recon_sp_relax(env::Geometry<Index> const& geo,
                                               env::Mechanism& m,
                                               system::SoA<Position const&, TypeID const&, Frozen const&> in) {
    //
    system::SoA<Position> recon = reconstruct_sp(geo, m, in);

    system::SoA<Position, Axis, TypeID const&, Frozen const&> dim(in.size());
    dim.rebind(id_, in);
    dim.rebind(fzn_, in);
    dim[r_] = recon[r_];

    std::normal_distribution<double> dist;

    auto& thr = thread();

    // Axis random and prop_to displacement.
    for (int j = 0; j < in.size(); ++j) {
      dim(ax_, j) = Vec::NullaryExpr([&] { return dist(thr.prng); }) * gnorm_sq(recon(r_, j) - in(r_, j));
    }
    // ... and randomise.
    dim[ax_] /= gnorm(dim[ax_]);

    auto err = thr.dimer.find_sp(dim, dim, in, thr.pot, {}, 1);
    auto lax = gnorm(recon[r_] - dim[r_]);

    if (err || lax > m_opt.sp_relax_tol) {
#pragma omp critical
      {
        if (m_opt.fout) {
          m_opt.fout->commit([&] { m_opt.fout->write(r_, recon); });
          m_opt.fout->commit([&] { m_opt.fout->write(r_, dim); });
        }
      }
      if (m_opt.debug) {
        fmt::print("Reconstructed saddle relaxing failed with: err={}, delta={}, cap_frac={}, disp0={}\n",
                   err,
                   lax,
                   m.capture_frac,
                   gnorm(m.delta_sp[0][del_]));
      }
      if (err) {
        throw error("Relaxing mech's SP failed with err={}", err);
      }
      m.poison = true;
    }

    return system::SoA<Position>(std::move(dim));
  }

}  // namespace fly::saddle