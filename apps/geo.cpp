
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
#include "libfly/utility/random.hpp"

using namespace fly;

template <typename... Ts>
system::Supercell<system::TypeMap<>, Ts...> supercell_from(std::string_view fname, std::uint64_t frame) {
  //
  io::BinaryFile file(fname, io::read_only);

  system::TypeMap out_map = file.read_map(frame);

  system::Box out_box = file.read_box(frame);

  system::Supercell out_cell = system::make_supercell<Ts...>(out_box, out_map, file.read<std::uint32_t>(frame, "particles/N"));

  auto read = [&](auto property) {
    try {
      file.read_to(frame, property, out_cell);
    } catch (fly::RuntimeError const &err) {
      fmt::print("Ignoring error, what(): {}\n", err.what());
    }
  };

  (read(Ts{}), ...);

  return out_cell;
}

#include "libfly/env/heuristics.hpp"

int main() {
  system::Supercell cell = supercell_from<Position, Frozen, PotentialGradient, Axis, Hash>("data/xyz/V1-unrelaxed.gsd", 0);

  system::Supercell old = cell;

  cell.destructive_resize(old.size() + 2);

  for (int i = 0; i < old.size(); i++) {
    cell(r_, i) = old(r_, i);
    cell(id_, i) = old(id_, i);
  }

  cell(r_, old.size()) = Vec{0.4, 1.4, 1.4};
  cell(id_, old.size()) = 1;

  cell(r_, old.size() + 1) = Vec{2.4, 1.4, 1.4};
  cell(id_, old.size() + 1) = 1;

  cell[fzn_] = false;
  cell[r_] += 0.2;

  //   Minimise.

  minimise::LBFGS minimiser({.debug = false, .fout = nullptr}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  fmt::print("FoundMin?={}\n", !timeit("Minimise", [&] { return minimiser.minimise(cell, cell, pot, omp_get_max_threads()); }));

  env::Catalogue cat({.delta_max = 500, .debug = false});

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

  mast.find_mechs({cat.get_geo(113)}, cell);

  /////////////////////////// IO ///////////////////////////

  return 0;
}
