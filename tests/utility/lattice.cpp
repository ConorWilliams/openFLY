// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY.
// If not, see <https://www.gnu.org/licenses/>.

#include "libfly/utility/lattice.hpp"

#include <catch2/catch.hpp>

#include "libfly/utility/core.hpp"

using namespace fly;

constexpr auto norm = [](auto const& a, auto const& b) { return gnorm(a - b); };

TEST_CASE("kruskal-1D", "[utility]") {
  //
  std::vector<Vec> points = {};

  REQUIRE(near(-1., kruskal_max(points, norm)));

  if constexpr (spatial_dims == 3) {
    //
    points.emplace_back(0, 0, 0);

    REQUIRE(near(-1., kruskal_max(points, norm)));

    points.emplace_back(1, 0, 0);

    REQUIRE(near(1., kruskal_max(points, norm)));

    points.emplace_back(3, 0, 0);

    REQUIRE(near(2., kruskal_max(points, norm)));

    points.emplace_back(0, 9, 0);

    REQUIRE(near(9., kruskal_max(points, norm)));
  }
}

TEST_CASE("kruskal-V3", "[utility]") {
  if constexpr (spatial_dims == 3) {
    //
    std::vector<Vec> points = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
    };

    REQUIRE(near(1., kruskal_max(points, norm)));
  }
}

TEST_CASE("kruskal-V4", "[utility]") {
  if constexpr (spatial_dims == 3) {
    //
    std::vector<Vec> points = {
        {0, 0, 0},
        {0, 1, 0},
        {0, 0, 9},
        {0, 1, 9},
    };

    REQUIRE(near(9., kruskal_max(points, norm)));
  }
}
