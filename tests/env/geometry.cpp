// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/geometry.hpp"

#include <fmt/core.h>

#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <random>

#include "libfly/system/VoS.hpp"
#include "libfly/system/property.hpp"
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

TEST_CASE("centroid", "[env]") {
  //
  system::VoS<Position> x = hypercube();

  Vec com = env::centroid(x);

  REQUIRE(gnorm(com - Vec::Constant(0.5)) < 1e-8);
}

TEST_CASE("rmsd", "[env]") {
  //
  system::VoS<Position> x = hypercube();

  REQUIRE(near(env::rmsd(x, x), 0.));

  auto y = x;

  Vec disp = Vec::Constant(1.);

  for (auto& elem : y) {
    elem[r_] += disp;
  }

  REQUIRE(near(env::rmsd(x, y), std::sqrt(static_cast<double>(y.size()) * gnorm_sq(disp))));
}

// See https://stackoverflow.com/a/38430739
static auto random_ortho(Xoshiro& rng) -> Mat {
  //
  std::normal_distribution gauss(0., 1.);
  //
  Mat O = Mat::NullaryExpr([&rng, &gauss]() { return gauss(rng); });

  Mat Q = Eigen::HouseholderQR<Mat>(O).householderQ();

  Mat Ip = Mat::Zero();

  for (int i = 0; i < spatial_dims; i++) {
    Ip(i, i) = gauss(rng) < 0 ? 1.0 : -1.0;
  }

  return Q * Ip;
}

TEST_CASE("ortho_onto", "[env]") {
  //
  system::VoS<Position> x = hypercube();

  std::random_device dev;

  Xoshiro rng({dev(), dev(), dev(), dev()});

  for (int i = 0; i < 100; i++) {
    //
    auto y = x;

    Mat O = random_ortho(rng);

    REQUIRE(gnorm(O.transpose() * O - Mat::Identity()) < 1e-8);

    for (auto& elem : y) {
      elem[r_] = O * elem[r_];
    }

    Mat Op = env::ortho_onto(x, y);

    REQUIRE(gnorm(Op - O) < 1e-8);

    //
  }
}

void make_centroid_origin(system::VoS<Position, Colour>& x) {
  //
  Vec com = env::centroid(x);

  for (auto& elem : x) {
    elem[r_] -= com;
  }

  REQUIRE(near(gnorm(env::centroid(x)), 0.0));
}

// Generate a random set of n atoms in a ball.
system::VoS<Position, Colour> rand_geo(Xoshiro& rng, int num_atoms) {
  //
  double rad = 5;

  std::uniform_real_distribution<> uni(-rad, rad);

  system::VoS<Position, Colour> out;

  out.emplace_back(Vec::Zero(), 0);

  while (out.size() < num_atoms) {
    //
    Vec p = Vec::NullaryExpr([&uni, &rng] { return uni(rng); });

    if (gnorm(p) < rad) {  // Rejection sampling.
      out.emplace_back(p, 0);
    }
  }

  make_centroid_origin(out);

  return out;
}

TEST_CASE("for_equiv_perms one", "[env]") {
  //

  std::random_device dev;

  Xoshiro rng({dev(), dev(), dev(), dev()});

  system::VoS<Position, Colour> ref = rand_geo(rng, 10);

  for (int i = 0; i < 1000; i++) {
    //
    Mat O = random_ortho(rng);

    system::VoS<Position, Colour> mut = ref;

    for (auto& elem : mut) {
      elem[r_] = O * elem[r_];
    }

    std::optional<Mat> info = {};

    env::for_equiv_perms(mut, ref, 1e-10, 1, [&](Mat const& Op, double) {
      info = Op;
      return false;
    });
    //
    REQUIRE(info);
    REQUIRE(gnorm(*info - O.transpose()) < 1e-10);
  }
}
