// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/boxes/orthorhombic.hpp"

#include <catch2/catch.hpp>
#include <optional>
#include <random>

#include "libfly/system/atom.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("Orthorhombic::canon_image", "[system]") {
  //
  using namespace fly;

  std::mt19937 gen(33);
  std::uniform_real_distribution<double> dis(0, 1);

  auto vrand = [&] { return Arr<double>::NullaryExpr([&]() { return dis(gen); }); };

  for (std::size_t i = 0; i < 100'000; i++) {
    //
    Arr<double> extents = vrand() + 1;
    Arr<bool> periodic = vrand() < .5;

    fly::system::Orthorhombic box{extents, periodic};

    // Random points inside simbox
    Vec a = vrand() * extents;
    Vec b = vrand() * extents;

    // Displace by integral number of random extents in each periodic periodic direction
    Vec b_prime = periodic.select(b.array() + (10 * vrand()).floor() * extents, b);

    // Check in same position
    REQUIRE(std::abs(gnorm(box.canon_image(b_prime) - a) - fly::gnorm(a - b)) < 0.001);
  }
}

TEST_CASE("OrthoGrid::gen_image", "[system]") {
  using namespace fly;

  system::Orthorhombic box{Arr<double>::Constant(10), Arr<bool>::Constant(true)};

  auto grid = box.make_grid(3);

  {
    std::optional im = grid.gen_image<Sign::plus>(Vec::Constant(5), 0);

    REQUIRE(!im);
  }

  {
    Vec origin = Vec::Constant(0);

    origin[0] += 0.2;

    REQUIRE(!grid.gen_image<Sign::minus>(origin, 0));

    auto im = grid.gen_image<Sign::plus>(origin, 0);

    REQUIRE(im);

    Vec correct = origin;

    correct[0] += 10;

    REQUIRE(gnorm(correct - *im) < 0.001);
  }

  {
    Vec origin = Vec::Constant(0);

    origin[0] = 10 - 0.2;

    REQUIRE(!grid.gen_image<Sign::plus>(origin, 0));

    auto im = grid.gen_image<Sign::minus>(origin, 0);

    REQUIRE(im);

    Vec correct = origin;

    correct[0] -= 10;

    REQUIRE(gnorm(correct - *im) < 0.001);
  }
}