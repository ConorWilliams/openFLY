// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/boxes/hypergrid.hpp"

#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <random>

#include "libfly/system/atom.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("HyperGrid::cell_idx", "[system]") {
  //
  using namespace fly;

  system::HyperGrid grid{Arr<Position::scalar_t>::Constant(10), 3};

  int A = grid.cell_idx(Vec<Position::scalar_t>::Constant(0));
  int B = grid.cell_idx(Vec<Position::scalar_t>::Constant(5));
  int C = grid.cell_idx(Vec<Position::scalar_t>::Constant(9.5));

  if constexpr (spatial_dims == 3) {
    REQUIRE(A == 1 + 1 * 5 + 1 * 5 * 5);
    REQUIRE(B == 2 + 2 * 5 + 2 * 5 * 5);
    REQUIRE(C == 3 + 3 * 5 + 3 * 5 * 5);
  }

  if constexpr (spatial_dims == 2) {
    REQUIRE(A == 1 + 1 * 5);
    REQUIRE(B == 2 + 2 * 5);
    REQUIRE(C == 3 + 3 * 5);
  }
}