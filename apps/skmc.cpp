
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
#include "libfly/kinetic/skmc.hpp"
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

  system::Supercell motif
      = system::make_supercell<Position, Frozen, T...>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

int main() {
  //

  system::Supercell cell = remove_atoms(motif_to_lattice(bcc_iron_motif<Hash>(), {6, 6, 6}), {1, 3});

  cell = add_atoms(cell,
                   {system::Atom<TypeID, Position, Frozen, Hash>(
                       1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14}, false, 0)});

  //   cell = add_atoms(
  //       cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 /
  //       4 * 2 + 3.14}, false, 0)});

  //   cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {0.992116, 6.01736, 4.56979},
  //   false, 0)});

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

  minimise::LBFGS minimiser({}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

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

  saddle::Master mast{
      {
          .num_threads = omp_get_max_threads(),
          .max_searches = 100,
          .max_failed_searches = 50,
          .debug = false,
      },
      cell.box(),
      pot,
      minimiser,
      saddle::Dimer{{}, {}, cell.box()},
  };

  kinetic::skmc(cell,
                minimiser,
                pot,
                mast,
                cat,
                {.debug = true, .opt_basin = {.debug = true, .temp = 300}, .opt_sb = {.debug = true}},
                omp_get_max_threads(),
                [&](system::SoA<Position const &> pre, system::SoA<Position const &> post) {
                  file.commit([&] { file.write(r_, pre); });
                  file.commit([&] { file.write(r_, post); });
                });

  return 0;
}
