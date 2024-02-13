// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/utility/spline.hpp"

#include <fmt/core.h>

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <random>
#include <vector>

#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

TEST_CASE("Spline basics", "[utility]") {
  //
  std::vector<double> f = {1, 1, 1, 1};

  fly::Spline sp(f, 1);

  REQUIRE(fly::near(sp.f(0), 1.0));
  REQUIRE(fly::near(sp.fp(0), 0.0));
  REQUIRE(fly::near(sp.fpp(0), 0.0));
}

TEST_CASE("Spline fuzzed", "[utility]") {
  //

  fly::Xoshiro gen({1, 3, 3, 4});

  std::uniform_real_distribution<double> dis(-1, 1);

  std::vector<double> tab(100);

  for (int i = 0; i < 100; i++) {
    //
    double x0 = dis(gen);
    double x1 = dis(gen);
    double x2 = dis(gen);
    double x3 = dis(gen);

    auto f = [=](double x) { return x0 * x * x * x + x1 * x * x + x2 * x + x3; };
    auto fp = [=](double x) { return 3 * x0 * x * x + 2 * x1 * x + x2; };
    auto fpp = [=](double x) { return 3 * 2 * x0 * x + 2 * x1; };

    double dx = std::abs(1.0 / static_cast<double>(tab.size()));

    for (std::size_t j = 0; j < tab.size(); j++) {
      tab[j] = f(static_cast<double>(j) * dx);
    }

    // Miss end points where derivatives deviate.
    std::uniform_real_distribution<double> dis2(0.1, 0.9);

    fly::Spline sp(tab, dx);

    for (int j = 0; j < 1000; j++) {
      //

      double x = dis2(gen);

      fmt::print("x={} {} {}\n", x, sp.f(x), f(x));
      fmt::print("x={} {} {}\n", x, sp.fp(x), fp(x));
      fmt::print("x={} {} {}\n", x, sp.fpp(x), fpp(x));

      REQUIRE(fly::near(sp.f(x), f(x), 1e-5, 0.01));
      REQUIRE(fly::near(sp.fp(x), fp(x), 1e-5, 0.01));
      REQUIRE(fly::near(sp.fpp(x), fpp(x), 1e-5, 0.01));
    }
  }
}
