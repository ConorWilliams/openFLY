// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/supercell.hpp"

#include <catch2/catch_test_macros.hpp>

#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("TypeMap + Supercell", "[system]") {
  using namespace fly;

  system::TypeMap<Mass, Position> map{4};

  map.set(tp_, 0, "Fe");

  CHECK(map.get(tp_, 0) == "Fe");

  map.set(r_, 0, {0, 1, 0});

  //   Slicing

  system::TypeMap<Position> pmap(map);

  CHECK(pmap.get(tp_, 0) == "Fe");
  CHECK(pmap.get(r_, 0) == Vec{0, 1, 0});

  // Supercell

  system::Box box(Mat::Identity(), {true, false, true});

  auto cell = system::make_supercell<Position>(box, map, 10);

  cell[r_] = 9;

  CHECK(cell.map().get(tp_, 0) == "Fe");
}
