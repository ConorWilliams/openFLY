
#include <fmt/core.h>
#include <omp.h>

//

#include <algorithm>
#include <array>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/kinetic/basin.hpp"
#include "libfly/kinetic/superbasin.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle//find.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/saddle/perturb.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

//

#include "update.hpp"

using namespace fly;

template <typename... T>
system::Supercell<system::TypeMap<>, Position, Frozen, T...> bcc_iron_motif() {
  //
  system::TypeMap<> FeH(2);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");

  Mat basis{
      {2.855300, 0.000000, 0.000000},
      {0.000000, 2.855300, 0.000000},
      {0.000000, 0.000000, 2.855300},
  };

  system::Supercell motif = system::make_supercell<Position, Frozen, T...>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

int main() {
  //

  system::Supercell cell = remove_atoms(motif_to_lattice(bcc_iron_motif<Hash>(), {6, 6, 6}), {1});

  cell = add_atoms(
      cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14}, false, 0)});

  //   cell = add_atoms(
  //       cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 * 2 + 3.14}, false,
  //       0)});

  //   cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {0.992116, 6.01736, 4.56979}, false, 0)});

  //   cell(r_, 431) = Vec{4.5896, 4.5892, 5.91909};
  //   cell(r_, 0) = Vec{4.45211, 4.45172, 4.2526};

  fly::io::BinaryFile file("build/gsd/sim.gsd", fly::io::create);

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());
    file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));
    file.write(id_, cell);

    file.write(r_, cell);
  });

  //   exit(0);

  minimise::LBFGS minimiser({}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  neigh::List neigh_list(cell.box(), pot.r_cut());

  {
    system::SoA<Position, PotentialGradient> out(cell.size());
    bool err = timeit("Minimise", [&] { return minimiser.minimise(out, cell, pot, omp_get_max_threads()); });
    cell[r_] = out[r_];
    verify(!err, "Minimiser failed");
  }

  file.commit([&] { file.write(r_, cell); });

  static constexpr char fname[] = "build/gsd/cat.bin";

  env::Catalogue cat = [&] {
    if (std::ifstream fcat(fname); fcat.good()) {
      fmt::print("Opening existing catalogue: \"{}\"\n", fname);
      return env::Catalogue{{}, fcat};
    } else {
      fmt::print("Could not open catalogue: \"{}\"\n", fname);
      return env::Catalogue{{}};
    }
  }();

  auto dump_cat = [&] {
    // exit(0);
    std::ofstream fcat(fname);
    cat.dump(fcat);
  };

  saddle::Master::Options opt{
      .num_threads = omp_get_max_threads(),
      .max_searches = 100,
      .max_failed_searches = 50,
      .debug = false,
  };

  saddle::Master mast{
      opt,
      cell.box(),
      pot,
      minimiser,
      saddle::Dimer{{}, {}, cell.box()},
  };
  //

  if (update_cat(mast, cat, cell, omp_get_max_threads())) {
    dump_cat();
  }

  std::random_device rd;

  Xoshiro rng(rd);

  double time = 0;

  system::SoA<Position, Frozen const&, TypeID const&> raw_recon(cell.size());
  raw_recon.rebind(fzn_, cell);
  raw_recon.rebind(id_, cell);
  system::SoA<Position, PotentialGradient> rel_recon(cell.size());

  auto energy = [&](system::SoA<Position const&> x) {
    neigh_list.rebuild(x, omp_get_max_threads());
    return pot.energy(cell, neigh_list, omp_get_max_threads());
  };

  double temp = 300;

  /**
   * @brief Build a new superbasin from the cell/catalogue.
   */
  auto new_sb = [&]() -> kinetic::SuperBasin {
    fmt::print("New superbasin!\n");
    return {{.debug = true}, {{.debug = true, .temp = temp}, cell, cat}};
  };

  kinetic::SuperBasin sb = new_sb();

  //   Trap detection!

  struct Traps {
  public:
    double try_push(double E, Vec h_atom) {
      //
      if (reset && E < E_trap) {
        fmt::print("Resetting traps\n");
        m_traps.clear();
        reset = false;
      }

      double min = std::numeric_limits<double>::max();

      for (auto const& elem : m_traps) {
        min = std::min(min, gnorm(h_atom - elem));
      }

      if (E < E_trap) {
        if (min < lim) {
          fmt::print("Known trap\n");
          return min;
        } else {
          fmt::print("New trap\n");
          m_traps.push_back(h_atom);
          return 0;
        }
      } else {
        fmt::print("H is {} from trap\n", min);
        // Not a trapping state return min.
        return min;
      }
    }

    // Once this function has been called at the next try_push that is a trapping state the traps will be cleared
    void reset_at_next_trap() { reset = true; }

    std::size_t size() const noexcept { return m_traps.size(); }

  private:
    bool reset = false;
    double lim = 0.1;
    double E_trap = -1731.56;  // v2 = -1726.05;  // v3 = -1720.68;  // V1 = -1731.56;
    std::vector<Vec> m_traps;
  };

  Traps traps;

  fmt::print("E0 = {}\n", energy(cell));

  double dtrap = traps.try_push(energy(cell), cell(r_, cell.size() - 1));

  for (int i = 0; dtrap < 6; i++) {
    timeit("TOTAL", [&] {
      ///////////// Select mechanism /////////////

      auto const& [m, atom, dt, basin, changed] = timeit("Super::choice", [&] { return sb.kmc_choice(rng); });

      if (changed) {
        fmt::print("Basin changed to {}\n", basin);
        cell[r_] = sb.current_basin().state()[r_];
        file.commit([&] { file.write(r_, cell); });
      } else {
        fmt::print("Basin did NOT changed, still at {}\n", basin);
      }

      ///////////// Reconstruct mech /////////////

      verify(!m.poison, "KMC chose a poisoned mechanisms with dE={}", m.barrier);

      double E0 = energy(cell);  // Energy before mechanism

      cat.reconstruct(raw_recon, m, atom, cell, !changed, omp_get_max_threads());

      auto err = timeit("Minimise", [&] { return minimiser.minimise(rel_recon, raw_recon, pot, omp_get_max_threads()); });

      double Ef = energy(rel_recon);  // Energy after relax

      double dR_err = std::abs(gnorm(rel_recon[r_] - raw_recon[r_]) - m.err_fwd);
      double dE_err = std::abs(Ef - E0 - m.delta);

      double dR_err_frac = dR_err / m.err_fwd;
      double dE_err_frac = dE_err / std::abs(m.delta);

      fmt::print("dE_err={:.3f}[{:.3f}], dR_err={:.3f}[{:.3f}]\n", dE_err, dE_err_frac, dR_err, dR_err_frac);

      bool fail = false;

      if (err) {
        fail = true;
      } else if (dE_err > opt.recon_e_tol_abs && dE_err_frac > opt.recon_e_tol_frac) {
        fail = true;
      } else if (dR_err > opt.recon_norm_frac_tol && dR_err_frac > opt.recon_norm_abs_tol) {
        fail = true;
      } else {
        cell[r_] = rel_recon[r_];
        file.commit([&] { file.write(r_, cell); });
      }

      if (fail) {
        fmt::print("Reconstruction failed @{}\n", atom);
        fmt::print("Failed reconstruction written to frame {}\n", file.n_frames());

        file.commit([&] { file.write(r_, rel_recon); });

        auto initial_assign = cat.get_ref(atom).cat_index();
        bool new_envs = false;

        do {
          double new_tol = cat.refine_tol(atom);
          fmt::print("Refined tolerance to {}\n", new_tol);
          new_envs = new_envs || update_cat(mast, cat, cell, omp_get_max_threads());
        } while (initial_assign == cat.get_ref(atom).cat_index());

        if (new_envs) {
          dump_cat();
        }

        sb = new_sb();

        // throw error("superbasin not ready for this|");

        return;  // Effectively continue.
      }

      ///////////// Time /////////////

      time += dt;

      ///////////// Custom /////////////

      fmt::print("E_eff={}\n", Ef);

      if (atom != cell.size() - 1) {
        traps.reset_at_next_trap();
      }

      dtrap = traps.try_push(Ef, cell(r_, cell.size() - 1));

      fmt::print("The superbasin has {} traps\n", traps.size());

      ///////////// Update catalogue /////////////

      if (timeit("Update catalogue", [&] { return update_cat(mast, cat, cell, omp_get_max_threads()); })) {
        dump_cat();
      }

      ///////////// Update SuperBasin /////////////

      if (sb.find_occupy(kinetic::hash(cell.size(), cat), cell, 0.1)) {
        fmt::print("Old basin in SuperBasin, size={}\n", sb.size());
      } else {
        sb.expand_occupy({{.debug = true, .temp = temp}, cell, cat});
        fmt::print("New basin in SuperBasin, size={}\n", sb.size());
      }

      sb.connect_from(basin, atom, m);
    });

    fmt::print("Iteration #{} time = {:.3e}, frames={}\n\n", i, time, file.n_frames());
  }

  return 0;
}
