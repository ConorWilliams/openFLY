
#include <fmt/core.h>
#include <omp.h>

#include <array>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <string_view>
#include <utility>
#include <variant>

#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/ROT/dimer.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle//find.hpp"
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

int main() {
  //

  system::Supercell cell = supercell_from<Position, Frozen, PotentialGradient, Axis>("data/xyz/V1-unrelaxed.gsd", 0);

  //   auto cell = make_super();

  // IO

  fly::io::BinaryFile fout("lbfgs.gsd", fly::io::create);

  fout.write(cell.box());                                                 //< Write the box to frame 0.
  fout.write(cell.map());                                                 //< Write the map to frame 0.
  fout.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));  //< Write the number of atoms to frame 0.
  fout.write(fly::id_, cell);                                             //< Write the TypeID's of the atoms to frame 0.

  //   WORK

  minimise::LBFGS::Options opt;

  opt.debug = true;
  opt.fout = &fout;

  minimise::LBFGS minimiser(opt, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  bool done = timeit("Minimise", [&] { return minimiser.minimise(cell, cell, pot, omp_get_max_threads()); });

  fmt::print("FoundMin?={}\n", !done);

  //   /////////////////////////////////////////////////

  system::Supercell<system::TypeMap<>, Position, Axis, Frozen, PotentialGradient> dcell(cell.box(), cell.map(), cell.size());

  dcell[r_] = cell[r_];
  dcell[fzn_] = cell[fzn_];
  dcell[id_] = cell[id_];
  dcell[ax_] = 1;
  dcell[ax_] /= gnorm(dcell[ax_]);

  std::random_device dev;

  Xoshiro rng({0, 0, 1, 3});

  saddle::perturb(dcell, rng, dcell.box(), dcell(r_, 113), dcell, 4, 0.6);

  fout.commit([&] {
    fout.write(fly::r_, dcell);  //< Write the positions of perturbed the atoms.
    fout.write(fly::ax_, dcell);
  });

  potential::Dimer::Options opt2;

  opt2.debug = true;
  //   opt2.relax_in_convex = false;

  potential::Generic dimer{
      potential::Dimer{
          opt2,
          pot,
      },
  };

  fly::system::SoA<Position, Axis> init{dcell};

  bool sp = timeit("warmup", [&] { return minimiser.minimise(dcell, dcell, dimer, omp_get_max_threads()); });

  for (size_t i = 0; i < 0; i++) {
    dcell[r_] = init[r_];
    dcell[ax_] = init[ax_];
    sp = timeit("FindSP", [&] { return minimiser.minimise(dcell, dcell, dimer, omp_get_max_threads()); });
  }

  fout.commit([&] {
    fout.write(fly::r_, dcell);  //< Write the position of the SP.
    fout.write(fly::ax_, dcell);
  });

  fmt::print("FoundSP?={}\n", !sp);

  return 0;
}
