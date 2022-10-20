// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/saddle/find.hpp"

#include <vector>

#include "libfly/system/SoA.hpp"
#include "libfly/utility/core.hpp"

namespace fly::saddle {

  MasterFinder::MasterFinder(Options const& opt, potential::Generic const& pot, minimise::LBFGS const& min, saddle::Dimer const& dimer)
      : m_opt{opt} {
    //
    std::random_device rd;

    Xoshiro prng({rd(), rd(), rd(), rd()});

    for (int i = 0; i < opt.num_threads; i++) {
      m_data.push_back({pot, min, dimer, prng});
      prng.long_jump();
    }
  }

  std::vector<MasterFinder::PathGroup> MasterFinder::find_pathways(system::Box const& box,
                                                                   std::vector<int> const& unknown,
                                                                   system::SoA<Position const&, Frozen const&, TypeID const&> in) {
    //
    neigh::List nl{box, m_opt.r_pert};

    nl.rebuild(in, m_opt.num_threads);

    std::vector<std::vector<Pathway>> paths(unknown.size());

    //
#pragma omp parallel
#pragma omp single nowait
    {
      for (std::size_t j = 0; j < paths.size(); ++j) {
        //
        int index = unknown[j];

#pragma omp task untied default(none) firstprivate(index, j) shared(in, nl, paths)
        {
          std::vector<std::optional<Pathway>> batch(safe_cast<std::size_t>(m_opt.batch_size));

          int tot = 0;
          int fail = 0;

          while (tot < m_opt.max_searches && fail < m_opt.max_failed_searches) {
            //  Do batch_size SP searches
            for (std::size_t i = 0; i < batch.size(); i++) {
#pragma omp task untied default(none) firstprivate(index, i) shared(in, nl, batch)
              { batch[i] = find_one(in, nl, index); }
            }

            // Wait for batch
#pragma omp taskwait
            // Process batch
            for (std::optional<Pathway>& elem : batch) {
              if (elem && push_if_new(paths[j], std::move(*elem))) {
                fail = 0;
              } else {
                fail++;
              }
            }

            tot += m_opt.batch_size;

            if (m_opt.debug) {
              fmt::print("FINDER: Found {} mech(s) at {}: fail={}, tot={}\n", paths[j].size(), index, fail, tot);
            }
          }
        }
      }
    }

    if (m_opt.fout) {
      for (auto&& m : paths) {
        for (auto&& v : m) {
          m_opt.fout->commit([&] { m_opt.fout->write(r_, in); });
          m_opt.fout->commit([&] { m_opt.fout->write(r_, v.sp); });
          m_opt.fout->commit([&] { m_opt.fout->write(r_, v.fwd); });

          if (m_opt.debug) {
            neigh::List nl2{box, thread().pot.r_cut()};

            nl2.rebuild(in, m_opt.num_threads);

            double E0 = thread().pot.energy(in, nl2, m_opt.num_threads);

            nl2.rebuild(v.sp, m_opt.num_threads);

            double Esp = thread().pot.energy(in, nl2, m_opt.num_threads);

            fmt::print("Delta={}eV\n", Esp - E0);
          }
        }
      }
    }

    std::vector<PathGroup> ret{};

    for (std::size_t i = 0; i < paths.size(); i++) {
      if (!paths[i].empty()) {
        ret.push_back({safe_cast<int>(i), std::move(paths[i])});
      }
    }

    return ret;
  }

  bool MasterFinder::push_if_new(std::vector<Pathway>& found, Pathway&& maybe) const noexcept {
    for (Pathway const& elem : found) {
      double dsp = gnorm(elem.sp[r_] - maybe.sp[r_]);     // SP displacement vector.
      double dfwd = gnorm(elem.fwd[r_] - maybe.fwd[r_]);  // FWD displacement vector.

      if (m_opt.debug) {
        fmt::print("FINDER: dsp={} dfwd={}\n", dsp, dfwd);
      }

      if (dsp < m_opt.basin_tol && dfwd < m_opt.basin_tol) {
        if (m_opt.debug) {
          fmt::print("FINDER: Found same SP as before\n");
        }
        return false;
      }
    }

    found.push_back(std::move(maybe));

    return true;
  }

  std::pair<int, int> MasterFinder::min_max(system::SoA<Position const&> a, system::SoA<Position const&> b) {
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

  static Vec com(system::SoA<Position const&> p) {
    Vec vsum = Vec::Zero();
    for (int i = 0; i < p.size(); i++) {
      vsum += p(r_, i);
    }
    return vsum / p.size();
  }

  std::optional<MasterFinder::Pathway> MasterFinder::find_one(system::SoA<Position const&, Frozen const&, TypeID const&> in,
                                                              neigh::List const& nl,
                                                              int index) {
    // Saddle search

    system::SoA<Position, Axis, Frozen const&, TypeID const&> dimer(in.size());

    dimer[r_] = in[r_];
    dimer.rebind(fzn_, in);
    dimer.rebind(id_, in);

    perturb(dimer, index, nl);

    ThreadData& thr = thread();

    if (thr.dimer.step(dimer, dimer, thr.pot, 1000, 1)) {
      return {};
    }

    auto [min, max] = min_max(dimer, in);

    if (max != index) {
      if (m_opt.debug) {
        fmt::print("FINDER: found mech at {} but wanted {}\n", max, index);
      }
      return {};
    }

    //   Minimisations

    if (double dr = gnorm(in(r_, min) - dimer(r_, min)); dr > m_opt.freeze_tol) {
      throw error("FINDER: Trying to freeze an atom {} that displaced {}", min, dr);
    } else if (m_opt.debug) {
      fmt::print("FINDER: Freezing atom #{}\n", min);
    }

    nl.for_neighbours(index, m_opt.r_pert, [min = min, index](auto n, double r, auto const&) {
      verify(n != min, "Frozen atom #{} is only {} from centre #{}", n, r, index);
    });

    double disp = gnorm(dimer[r_] - in[r_]);

    system::SoA<Position, PotentialGradient, Frozen, TypeID const&> relax{dimer.size()};
    relax[r_] = dimer[r_] + dimer[ax_] * disp * m_opt.nudge_frac;
    relax[fzn_] = dimer[fzn_];
    relax(fzn_, min) = true;  // Freeze minimally displaced atom.
    relax.rebind(id_, in);

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
    }

    // We now have a min->sp->min pathway with no translation (due to freeze). Need to correct for COM drift during SPS.

    Vec in_com = com(in);
    Vec rev_com = com(rev);

    Vec drift = rev_com - in_com;

    if (m_opt.debug) {
      fmt::print("FINDER: Centre of mass drift={}\n", gnorm(drift));
    }

    for (int i = 0; i < in.size(); i++) {
      if (!in(fzn_, i)) {
        fwd(r_, i) += drift;
        dimer(r_, i) += drift;
        rev(r_, i) += drift;
      }
    }

    // /////////////// CHECK pathway is min(in)->sp->min

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

    return Pathway(dimer, fwd);
  }

  auto MasterFinder::perturb(system::SoA<Position&, Axis&, Frozen const&> inout, int centre, neigh::List const& nl) -> void {
    //
    Xoshiro& prng = thread().prng;

    std::normal_distribution normal(0., 1.);

    std::normal_distribution prime(m_opt.stddev, m_opt.stddev / 3.0);

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

}  // namespace fly::saddle