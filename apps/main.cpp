
#include <fmt/core.h>

#include <array>
#include <ctime>
#include <fstream>
#include <iostream>
#include <utility>
#include <variant>

#include "libfly/io/gsd.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/utility/core.hpp"

using namespace fly;

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

  Arr<int> shape = Arr<int>::Constant(7);  // In unit cells

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

  std::ifstream eam_tab{"data/wen.eam.fs"};

  system::Supercell cell = make_super();

  potential::Generic pot(potential::EAM{cell.map(), std::make_shared<potential::DataEAM>(std::move(eam_tab))});

  fmt::print("r_cut={}\n", pot.r_cut());

  neigh::List nl(cell.box(), pot.r_cut());

  timeit("rebuild", [&] { nl.rebuild(cell, omp_get_max_threads()); });

  pot.energy(cell, nl, 4);

  timeit("grad", [&] { pot.gradient(cell, cell, nl, omp_get_max_threads()); });

  //   potential::EAM p2{cell.map(), std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"})};

  //   p2.gradient(cell.soa(), nl, 8);

  fmt::print("pot={}\n", cell[g_].head(10));

  //   std::unique_ptr<potential::Base> pot = std::make_unique<potential::EAM>(cell.map(), data);

  fmt::print("Done\n");

  return 0;
}
