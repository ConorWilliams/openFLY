// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/heuristics.hpp"

#include <fmt/core.h>

#include <catch2/catch.hpp>
#include <random>

#include "libfly/env/geometry.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

using namespace fly;

static env::Geometry<Index> hypercube() {
  //
  env::Geometry<Index> x;

  int i = 0;

  // Put atom at each corner of a hypercube.
  template_for<int>(Arr<int>::Zero(), Arr<int>::Constant(2), [&x, &i](auto... is) {
    //
    Vec point = Vec{static_cast<double>(is)...};

    x.emplace_back(point, 0, i++);
  });

  x.centre();

  return x;
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

// Shuffle + transform geometry
static void shuffle_rotate_geo(env::Geometry<Index>& geo, Xoshiro& rng) {
  //
  Mat O = random_ortho(rng);

  for (auto& elem : geo) {
    elem[r_] = O * elem[r_];
  }

  std::shuffle(geo.begin() + 1, geo.end(), rng);
}

TEST_CASE("Fingerprint", "[env]") {
  //
  std::random_device dev;

  Xoshiro rng(dev);

  env::Geometry cell = hypercube();
  env::Geometry copy = cell;

  env::Fingerprint f1;
  env::Fingerprint f2;

  f1.rebuild(cell);

  for (size_t i = 0; i < 1000; i++) {
    copy = cell;
    shuffle_rotate_geo(copy, rng);
    f2.rebuild(copy);
    REQUIRE(f1.equiv(f2, 1e-10));
  }
}

TEST_CASE("nauty", "[env]") {
  //
  std::random_device dev;

  Xoshiro rng(dev);

  env::Geometry ref = hypercube();

  auto h1 = env::canon_hash(ref, 1.1, 2);

  env::Geometry copy = ref;

  for (size_t i = 0; i < 1000; i++) {
    copy = ref;
    shuffle_rotate_geo(copy, rng);

    auto h2 = env::canon_hash(copy, 1.1, 2);

    REQUIRE(h1 == h2);

    Mat O = env::ortho_onto(ref, copy);
    double rmsd = env::grmsd(O, ref, copy);
    REQUIRE(rmsd < 1e-10);
  }
}