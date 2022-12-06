// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This eam_tab is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

#include "libfly/potential/EAM/eam.hpp"

#include <Eigen/src/Core/Diagonal.h>
#include <fmt/core.h>

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <memory>

#include "libfly/potential/EAM/data.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
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

  auto cell = fly::system::make_supercell<Position, Frozen, PotentialGradient>(box, map, fly::ssize(lat));

  for (int i = 0; i < fly::ssize(lat); i++) {
    cell(r_, i) = lat[safe_cast<std::size_t>(i)].off;
    cell(id_, i) = lat[safe_cast<std::size_t>(i)].num;
    cell(fzn_, i) = false;
  }

  return cell;
}

// Testing away from a minima
TEST_CASE("EAM compute", "[potential]") {
  //
  std::ifstream eam_tab{"../../../data/wen.eam.fs"};

  if (eam_tab.good()) {
    //
    auto data = std::make_shared<potential::DataEAM>(std::move(eam_tab));

    system::Supercell cell = make_super();

    auto pot = std::make_unique<potential::EAM>(cell.map(), data);

    neigh::List nl(cell.box(), pot->r_cut());

    timeit("rebuild", [&] { nl.rebuild(cell, omp_get_max_threads()); });

    double E0;

    timeit("energy", [&] { E0 = pot->energy(cell, nl, omp_get_max_threads()); });

    REQUIRE(near(E0, -2748.5339306065985));  // Known good data for 3D

    timeit("potential", [&] { pot->gradient(cell, cell, nl, omp_get_max_threads()); });

    Vec head = cell(g_, 0);

    REQUIRE(gnorm(head - Vec::Constant(-0.19182439085483005)) < 1e-10);

    system::Hessian H;

    timeit("hessian", [&] { pot->hessian(H, cell, nl, omp_get_max_threads()); });

    Mat ngd = Mat{
        {+9.0654715551510830, -1.0875103661087880, -1.0875103661087877},
        {-1.0875103661087880, +9.0654715551510830, -1.0875103661087873},
        {-1.0875103661087877, -1.0875103661087873, +9.0654715551510830},
    };

    // ngd /= 55.847;  // Mass weighted.

    std::cout << H(0, 0) << std::endl;

    REQUIRE(gnorm(ngd - H(0, 0)) < 1e-10);

  } else {
    fmt::print("WARN - tests could not open eam file\n");
  }
}

// Testing at a minima
TEST_CASE("EAM hess", "[potential]") {
  //
  std::ifstream eam_tab{"../../../data/wen.eam.fs"};

  if (eam_tab.good()) {
    //
    auto data = std::make_shared<potential::DataEAM>(std::move(eam_tab));

    system::Supercell cell = make_super(false);

    auto pot = std::make_unique<potential::EAM>(cell.map(), data);

    neigh::List nl(cell.box(), pot->r_cut());

    timeit("rebuild", [&] { nl.rebuild(cell, omp_get_max_threads()); });

    timeit("potential", [&] { pot->gradient(cell, cell, nl, omp_get_max_threads()); });

    REQUIRE(gnorm(cell[g_]) < 1e-7);

    system::Hessian H;

    timeit("hessian", [&] { pot->hessian(H, cell, nl, omp_get_max_threads()); });

    auto const& ev = H.eigenvalues();

    fmt::print("eigen_values = {}", ev.head(10));

    for (int i = 0; i < spatial_dims; i++) {
      REQUIRE(near(ev[i], 0.0));
    }

    for (int i = spatial_dims; i < ev.size(); i++) {
      REQUIRE(!near(ev[i], 0.0));
    }

  } else {
    fmt::print("WARN - tests could not open eam file\n");
  }
}