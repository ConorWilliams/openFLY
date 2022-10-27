
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

system::Supercell<system::TypeMap<>, Position, Frozen> bcc_iron_motif() {
  //
  fly::system::TypeMap<> FeH(2);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");

  Mat basis{
      {2.855300, 0.000000, 0.000000},
      {0.000000, 2.855300, 0.000000},
      {0.000000, 0.000000, 2.855300},
  };

  system::Supercell motif = system::make_supercell<Position, Frozen>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

system::Supercell<system::TypeMap<>, Position, Frozen> fcc_iron_motif() {
  //
  fly::system::TypeMap<> FeH(2);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");

  Mat basis{
      {3.650000, 0.000000, 0.000000},
      {0.000000, 3.650000, 0.000000},
      {0.000000, 0.000000, 3.650000},
  };

  system::Supercell motif = system::make_supercell<Position, Frozen>({basis, Arr<bool>::Constant(true)}, FeH, 4);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec{0.0, 0.5, 0.5};
  motif(r_, 2) = Vec{0.5, 0.0, 0.5};
  motif(r_, 3) = Vec{0.5, 0.5, 0.0};

  return motif;
}

int main() {
  //

  system::Supercell cell0 = motif_to_lattice(fcc_iron_motif(), {5, 5, 5});

  system::Supercell cell = remove_atoms(cell0, {172});

  fly::io::BinaryFile tmp("test.gsd", fly::io::create);

  tmp.commit([&] {
    tmp.write(cell.box());
    tmp.write(cell.map());
    tmp.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));
    tmp.write(id_, cell);

    tmp.write(r_, cell);
  });

  fly::io::BinaryFile fout("lbfgs.gsd", fly::io::create);

  //   WORK

  minimise::LBFGS minimiser({.debug = true, .fout = &tmp}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  system::SoA<Position, PotentialGradient> out(cell.size());

  bool done = timeit("Minimise", [&] { return minimiser.minimise(out, cell, pot, omp_get_max_threads()); });

  fmt::print("FoundMin?={}\n", !done);

  tmp.commit([&] { tmp.write(r_, out); });

  return 0;
}
