
#include <fmt/core.h>
#include <omp.h>

#include <array>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/env/heuristics.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle//find.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/saddle/find.hpp"
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

system::Supercell<system::TypeMap<>, PotentialGradient, Position, Frozen, Hash> bcc_motif() {
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
      = system::make_supercell<PotentialGradient, Position, Frozen, Hash>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;
  motif[hash_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

int main() {
  system::Supercell cell = remove_atoms(motif_to_lattice(bcc_motif(), {6, 6, 6}), {1});

  //   system::Supercell old = cell;

  //   cell.destructive_resize(old.size() + 2);

  //   for (int i = 0; i < old.size(); i++) {
  //     cell(r_, i) = old(r_, i);
  //     cell(id_, i) = old(id_, i);
  //   }

  //   cell(r_, old.size()) = Vec{0.4, 1.4, 1.4};
  //   cell(id_, old.size()) = 1;

  //   cell(r_, old.size() + 1) = Vec{2.4, 1.4, 1.4};
  //   cell(id_, old.size() + 1) = 1;

  //   cell[fzn_] = false;
  //   cell[r_] += 0.2;

  //   Minimise.

  minimise::LBFGS minimiser({.debug = false, .fout = nullptr}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  fmt::print("FoundMin?={}\n", !timeit("Minimise", [&] { return minimiser.minimise(cell, cell, pot, omp_get_max_threads()); }));

  env::Catalogue cat({.delta_max = 0.1, .debug = false});

  std::vector ix = timeit("cat.rebuild()", [&] { return cat.rebuild(cell, omp_get_max_threads()); });

  for (auto &&elem : ix) {
    fmt::print("New env @{}\n", elem);
  }

  for (int i = 0; i < cell.size(); ++i) {
    cell(hash_, i) = cat.get_ref(i).cat_index();
  }

  fly::io::BinaryFile fout("geo.gsd", fly::io::create);

  fout.commit([&] {
    fout.write(cell.box());                                                 //< Write the box to frame 0.
    fout.write(cell.map());                                                 //< Write the map to frame 0.
    fout.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));  //< Write the number of atoms to frame 0.
    fout.write(fly::id_, cell);                                             //< Write the TypeID's of the atoms to frame 0.
    fout.write(fly::r_, cell);                                              //< Write the TypeID's of the atoms to frame 0.
    fout.write(fly::hash_, cell);                                           //< Write the TypeID's of the atoms to frame 0.
  });

  // ////////////////////// find 2 ////////////////////////

  fly::io::BinaryFile dout("finder.gsd", fly::io::create);

  dout.write(cell.box());                                                 //< Write the box to frame 0.
  dout.write(cell.map());                                                 //< Write the map to frame 0.
  dout.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));  //< Write the number of atoms to frame 0.
  dout.write(fly::id_, cell);
  dout.write(fly::r_, cell);
  dout.commit();

  std::vector<env::Geometry<Index>> tmp;

  for (auto const &elem : ix) {
    tmp.push_back(cat.get_geo(elem));
  }

  saddle::Dimer dimer{
      {},  //{.debug = true, .fout = &dout},
      {},  //{.debug = true},
      cell.box(),
  };

  saddle::Master mast{
      {.num_threads = omp_get_max_threads(), .debug = true, .fout = &dout},
      cell.box(),
      pot,
      minimiser,
      dimer,
  };

  mast.find_mechs({cat.get_geo(2)}, cell);

  /////////////////////////// IO ///////////////////////////

  return 0;
}
