
#include <fmt/core.h>
#include <omp.h>

#include <array>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string_view>
#include <utility>
#include <variant>

#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

using namespace fly;

template <typename... Ts>
system::Supercell<system::TypeMap<>, Ts...> supercell_from(std::string_view fname, std::uint64_t frame) {
  //
  io::BinaryFile file(fname, io::read_only);

  system::TypeMap out_map = file.read_map(frame);

  system::Box out_box = file.read_box(frame);

  system::Supercell out_cell = system::make_supercell<Ts...>(out_box, out_map, file.read<std::uint32_t>(frame, "particles/N"));

  (file.read_to(frame, Ts{}, out_cell), ...);

  return out_cell;
}

auto make_super(bool erase = true) {
  // Generate a BBC lattice, from old code base.

  struct MotifPt {
    TypeID::scalar_t num;
    Vec off;
  };

  // Fractional motif
  std::array BCC = {
      MotifPt{0, Vec::Zero()},
      MotifPt{0, Vec::Constant(0.5)},
  };

  double T = 300.0;

  double V = 11.64012 + T * (9.37798e-5 + T * (3.643134e-7 + T * (1.851593e-10 + T * 5.669148e-14)));

  double a = std::pow(2 * V, 1.0 / 3);  // lat param

  std::vector<MotifPt> lat;

  Arr<int> shape = Arr<int>::Constant(15);  // In unit cells

  template_for<int>(Arr<int>::Zero(), shape, [&](auto... i) {
    for (auto const& mot : BCC) {
      //

      Vec lp = Arr<int>{i...}.cast<double>().matrix() + mot.off;

      lat.push_back({mot.num, lp * a});
    }
  });

  if (erase) {
    lat.erase(lat.begin() + 1);
  }

  system::TypeMap<> map{2};

  map.set(0, "Fe");
  map.set(1, "H");

  Mat basis = Mat::Zero();

  for (int i = 0; i < spatial_dims; i++) {
    basis(i, i) = a * shape[i];
  }

  auto box = fly::system::Box(basis, Arr<bool>::Constant(true));

  auto cell = fly::system::make_supercell<Position, Frozen, PotentialGradient>(box, map, xise(lat));

  for (int i = 0; i < xise(lat); i++) {
    cell(r_, i) = lat[safe_cast<std::size_t>(i)].off;
    cell(id_, i) = lat[safe_cast<std::size_t>(i)].num;
    cell(fzn_, i) = false;
  }

  return cell;
}

int main() {
  //

  system::Supercell cell = supercell_from<Position, Frozen>("data/xyz/V1-unrelaxed.gsd", 0);

  cell(fzn_, 0) = true;
  cell(fzn_, 113) = true;

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
  //   opt.fout = &fout;

  opt.skin_frac = 1.05;

  minimise::LBFGS minimiser(opt, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  bool done = timeit("Minimise", [&] {
    //
    return minimiser.minimise(cell, cell, pot, omp_get_max_threads());
  });

  fmt::print("FoundMin?={}\n", done);

  return 0;
}
