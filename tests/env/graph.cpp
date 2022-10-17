// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "../../src/env/graph.hpp"

#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "libfly/env/geometry.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"

using namespace fly;

TEST_CASE("offsets", "[env]") {
  env::Geometry<> geo;

  geo.emplace_back(Vec::Zero(), 2);
  geo.emplace_back(Vec::Constant(1), 0);
  geo.emplace_back(Vec::Constant(2), 2);
  geo.emplace_back(Vec::Constant(3), 0);
  geo.emplace_back(Vec::Constant(4), 2);
  geo.emplace_back(Vec::Constant(5), 2);

  std::vector<std::size_t> off = env::offsets(geo, 4);

  REQUIRE(off == std::vector<std::size_t>{1, 3, 3, 6, 2});
}