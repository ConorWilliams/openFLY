
#include <fmt/core.h>
#include <omp.h>

//

#include <algorithm>
#include <array>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <string_view>
#include <utility>
#include <variant>

#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/kinetic/basin.hpp"
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

struct Update {
  bool new_envs = false;  ///< True if new environments encountered.
  bool ref_envs = false;  ///< True if environment refinement occurred.
};

// Returns true it catalogue has been refined.
template <typename Map, typename... T>
Update update_cat(saddle::Master& mast, env::Catalogue& cat, system::Supercell<Map, T...> const& cell) {
  //

  std::vector<int> ix = timeit("cat.rebuild()", [&] { return cat.rebuild(cell, omp_get_max_threads()); });

  if (ix.empty()) {
    return {};
  }

  int refines = 0;

  while (true) {
    //
    std::vector<int> fails;

    fmt::print("New envs @{} with {} refines\n", ix, refines);

    std::vector found = mast.find_mechs(mast.package(ix, cat), cell);

    /*
     * If find_mechs has failed, the failed environments must be too symmetric, we must refine them until they are less symmetric.
     */

    for (std::size_t i = 0; i < found.size(); i++) {
      if (found[i]) {
        cat.set_mechs(ix[i], found[i].mechs());
      } else {
        fails.push_back(ix[i]);
      }
    }

    if (fails.empty()) {
      return {true, refines == 0};
    } else {
      ++refines;
    }

    // Refine tolerance's
    for (auto const& f : fails) {
      auto n = cat.calc_self_syms(f).size();
      do {
        cat.refine_tol(f, cat.get_ref(f).delta_max() / 1.5);
      } while (n == cat.calc_self_syms(f).size());
    }

    std::vector tmp = cat.rebuild(cell, omp_get_max_threads());

    // The atom whose symmetry tolerance increased still needs to be searched alongside any atoms that now no longer match that
    // environment.
    for (auto const& elem : tmp) {
      verify(std::find(fails.begin(), fails.end(), elem) == fails.end(), "Atom #{} found twice", elem);
    }

    // tmp now stores previous fails + new environments that don't match the refined tolerance.
    tmp.insert(tmp.end(), fails.begin(), fails.end());

    ix = std::move(tmp);
  }
};

struct Options {};

int main() {
  //

  system::Supercell cell = remove_atoms(motif_to_lattice(bcc_iron_motif<Hash>(), {6, 6, 6}), {1});

  //   cell
  //       = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 0.5 + 3.14}, false,
  //       0)});

  fly::io::BinaryFile file("build/gsd/sim.gsd", fly::io::create);

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());
    file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));
    file.write(id_, cell);

    file.write(r_, cell);
  });

  minimise::LBFGS minimiser({}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  neigh::List neigh_list(cell.box(), pot.r_cut());

  system::SoA<Position, PotentialGradient> out(cell.size());

  auto do_min = [&] {
    bool err = timeit("Minimise", [&] { return minimiser.minimise(out, cell, pot, omp_get_max_threads()); });
    verify(!err, "Minimiser failed");
    cell[r_] = out[r_];
  };

  do_min();

  file.commit([&] { file.write(r_, cell); });

  env::Catalogue cat = [&] {
    if (std::ifstream fcat("build/gsd/cat.h.bin"); fcat.good()) {
      return env::Catalogue{{}, fcat};
    } else {
      return env::Catalogue{{}};
    }
  }();

  auto dump_cat = [&] {
    std::ofstream fcat("build/gsd/cat.h.bin");
    cat.dump(fcat);
  };

  saddle::Master mast{
      {.num_threads = omp_get_max_threads(), .max_searches = 500, .max_failed_searches = 125},
      cell.box(),
      pot,
      minimiser,
      saddle::Dimer{{}, {}, cell.box()},
  };

  auto [new_env, _] = update_cat(mast, cat, cell);

  if (new_env) {
    dump_cat();
  }

  std::random_device rd;

  Xoshiro rng(rd);

  double time = 0;

  for (int i = 0; i < 1000; i++) {
    timeit("TOTAL\n", [&] {
      ///////////// Select mechanism /////////////

      kinetic::Basin basin = timeit("kinetic::Basin", [&] { return kinetic::Basin({.debug = true}, cell, cat); });

      auto const& [m, atom, dt] = basin.kmc_choice(rng);

      time += dt;

      ///////////// Reconstruct mech /////////////

      // Energy before mechanism
      neigh_list.rebuild(cell, omp_get_max_threads());
      double E0 = pot.energy(cell, neigh_list, omp_get_max_threads());

      // recon
      cat.reconstruct(cell, m, atom, cell, true, omp_get_max_threads());

      // Energy after reconstruction
      neigh_list.rebuild(cell, omp_get_max_threads());
      double Etmp = pot.energy(cell, neigh_list, omp_get_max_threads());

      do_min();

      // Energy after relax
      neigh_list.rebuild(cell, omp_get_max_threads());
      double Ef = pot.energy(cell, neigh_list, omp_get_max_threads());

      ///////////// Output data/logs /////////////

      fmt::print("E0={}, Etmp={}, Ef={}\n", E0, Etmp, Ef);

      file.commit([&] { file.write(r_, cell); });

      ///////////// Update catalogue /////////////

      auto [new_envs, refined] = timeit("Update catalogue", [&] { return update_cat(mast, cat, cell); });

      if (new_envs) {
        dump_cat();
      }
    });
  }

  return 0;
}
