// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/local.hpp"

#include <fmt/core.h>

#include <catch2/catch_test_macros.hpp>
#include <random>

#include "libfly/env/geometry.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

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

/**
 * @brief get a supercell
 *
 * All atoms are the same type and not frozen.
 * There are two types (A and B).
 *
 * Atoms are on the corners of a hypercube with side length 1.
 *
 * The box is cubic and non periodic 10 x 10 x ...
 *
 * @return auto
 */
static auto build_hc_cell() {
  system::TypeMap<> map(2);

  map.set(0, tp_, "A");
  map.set(1, tp_, "B");

  system::Box box(10 * Mat::Identity(), Arr<bool>::Constant(false));

  system::VoS cube = hypercube();

  auto cell = system::make_supercell<Position, Frozen>(box, map, cube.size());

  cell[fzn_] = false;
  cell[id_] = 0;

  for (int i = 0; i < cell.size(); i++) {
    cell(r_, i) = cube[i][r_];
  }

  return cell;
}

// static auto box() {
//   system::TypeMap<> map(2);

//   map.set(0, tp_, "A");
//   map.set(1, tp_, "B");

//   system::Box box(10 * Mat::Identity(), Arr<bool>::Constant(false));

//   auto cell = system::make_supercell<Position, Frozen>(box, map, 4);

//   cell[fzn_] = false;
//   cell[id_] = 0;

//   cell(r_, 0) = Vec{0, 0, 0};
//   cell(r_, 1) = Vec{1, 0, 0};
//   cell(r_, 2) = Vec{0, 1, 0};
//   cell(r_, 3) = Vec{1, 1, 0};

//   return cell;
// }

TEST_CASE("Local", "[env]") {
  //

  system::Supercell cell = build_hc_cell();

  neigh::List nl(cell.box(), 5.2);

  nl.rebuild(cell, omp_get_max_threads());

  env::Local le1;

  le1.rebuild(0, cell, nl, cell.map().num_types(), 5.2, 1.1);

  fmt::print("LE1:\n");
  for (auto&& elem : le1) {
    fmt::print("{}\n", elem[r_]);
  }

  //   REQUIRE(le1.size() == cell.size());

  env::Local le2;

  le2.rebuild(3, cell, nl, cell.map().num_types(), 5.2, 1.1);

  fmt::print("LE2:\n");
  for (auto&& elem : le2) {
    fmt::print("{}\n", elem[r_]);
  }

  double delta = 0.1;

  REQUIRE(le1.fingerprint().equiv(le2.fingerprint(), delta));
  REQUIRE(le2.key() == le1.key());

  Mat O = env::ortho_onto(le1, le2);

  double rmsd = env::grmsd(O, le1, le2);

  REQUIRE(rmsd < 1e-10);
}
