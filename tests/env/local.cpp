// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/local.hpp"

#include <catch2/catch_test_macros.hpp>

#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"

using namespace fly;

static system::VoS<Position> hypercube() {
  //
  system::VoS<Position> x;
  // Put atom at each corner of a hypercube.
  template_for<int>(Arr<int>::Zero(), Arr<int>::Constant(2), [&x](auto... is) {
    //
    Vec point = Vec{static_cast<double>(is)...};

    x.emplace_back(point);
  });

  return x;
}

TEST_CASE("Local", "[env]") {
  //
  system::TypeMap<> map(1);

  map.set(0, tp_, "Null");

  system::Box box(10 * Mat::Identity(), Arr<bool>::Constant(false));

  system::VoS cube = hypercube();

  auto cell = system::make_supercell<Position, Frozen>(box, map, cube.size());

  cell[fzn_] = false;
  cell[id_] = 0;

  for (int i = 0; i < cell.size(); i++) {
    cell(r_, i) = cube[i][r_];
  }

  neigh::List nl(cell.box(), 5.2);

  nl.rebuild(cell, omp_get_max_threads());

  env::Local le1;

  le1.rebuild(0, cell, nl, cell.map().num_types(), 5.2, 1.1);

  REQUIRE(le1.size() == cell.size());

  env::Local le2;

  le2.rebuild(3, cell, nl, cell.map().num_types(), 5.2, 1.1);

  double delta = 0.1;

  REQUIRE(le1.fingerprint().equiv(le2.fingerprint(), delta));

  std::optional res = le1.permute_onto(le2, delta);

  REQUIRE(res);
  REQUIRE(res->rmsd < 1e-10);
  REQUIRE(le2.key() == le1.key());
}