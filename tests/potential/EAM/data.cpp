// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

#include "libfly/potential/EAM/data.hpp"

#include <fmt/core.h>

#include <catch2/catch.hpp>

#include "libfly/utility/core.hpp"

using namespace fly;

TEST_CASE("EAM parsing", "[potential]") {
  //
  std::ifstream file{"../../../data/wen.eam.fs"};

  if (file.good()) {
    potential::DataEAM data{{}, std::move(file)};

    REQUIRE(data.type_map().num_types() == 2);

    auto const& map = data.type_map();

    CHECK(map.get(0, tp_) == "Fe");
    CHECK(map.get(1, tp_) == "H");

    CHECK(near(data.f(0).f(0), -0.1414216239132407E+00));
    CHECK(near(data.f(0).f(1.9999999999999999E-02), -0.2000010747904311E+00));
  } else {
    fmt::print("WARN - tests could not open eam file\n");
  }
}