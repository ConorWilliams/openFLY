
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

#include "libfly/io/gsd.hpp"
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
    } catch (fly::RuntimeError const& err) {
      fmt::print("Ignoring error, what(): {}\n", err.what());
    }
  };

  (read(Ts{}), ...);

  return out_cell;
}

#include "libfly/env/local.hpp"

int main() {
  system::Supercell cell = supercell_from<Position, Frozen, PotentialGradient, Axis, Hash>("data/xyz/V1-unrelaxed.gsd", 0);

  //   Minimise.

  minimise::LBFGS minimiser({.debug = false, .fout = nullptr}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  bool done = timeit("Minimise", [&] { return minimiser.minimise(cell, cell, pot, omp_get_max_threads()); });

  fmt::print("FoundMin?={}\n", !done);

  double r_env = 5.2;

  neigh::List nl(cell.box(), r_env);

  nl.rebuild(cell, omp_get_max_threads());

  //   env::Local le1;

  Vector<env::Local> les(cell.size());

  std::map<std::size_t, std::size_t> mp;

  int count = 0;

  for (int i = 0; i < cell.size(); i++) {
    les[i].rebuild(i, cell, nl, cell.map().num_types(), r_env, 3);

    auto h = les[i].key();

    if (mp.count(h) == 0) {
      mp[h] = count++;
    }

    cell(hash_, i) = mp[h];
  }

  // IO

  fly::io::BinaryFile fout("geo.gsd", fly::io::create);

  fout.commit([&] {
    fout.write(cell.box());                                                 //< Write the box to frame 0.
    fout.write(cell.map());                                                 //< Write the map to frame 0.
    fout.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));  //< Write the number of atoms to frame 0.
    fout.write(fly::id_, cell);                                             //< Write the TypeID's of the atoms to frame 0.
    fout.write(fly::r_, cell);                                              //< Write the TypeID's of the atoms to frame 0.
    fout.write(fly::hash_, cell);                                           //< Write the TypeID's of the atoms to frame 0.
  });

  return 0;
}
