// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/boxes/triclinic.hpp"

#include <catch2/catch_test_macros.hpp>
#include <ios>
#include <iostream>
#include <random>

#include "libfly/system/atom.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("Triclinic::canon_image", "[system]") {
  //
  using namespace fly;

  std::mt19937 gen(33);
  std::uniform_real_distribution<double> dis(0, 1);

  // Random numbers on interval 0 to 1
  auto vrand = [&] { return Vec::NullaryExpr([&]() { return dis(gen); }); };

  using M = Mat;

  for (std::size_t i = 0; i < 1'000; i++) {
    //
    M basis = [&]() {
      M tmp = M::NullaryExpr([&]() { return dis(gen); });
      tmp.array() += 1;
      tmp = tmp.triangularView<Eigen::Upper>();
      return tmp;
    }();

    Arr<bool> periodic = vrand().array() > .5;

    fly::system::Triclinic box{basis, periodic};

    // // Random points inside simbox
    Vec a = basis * vrand();
    Vec b = basis * vrand();

    Vec integer_randoms = (10 * vrand()).array().floor();

    Eigen::DiagonalMatrix<double, spatial_dims, spatial_dims> diag(integer_randoms);

    M offsets = basis * diag;

    Vec b_prime = b;

    for (int j = 0; j < spatial_dims; j++) {
      if (periodic[j]) {
        b_prime += offsets.col(j);
      }
    }

    Vec b_undo = box.canon_image(b_prime);

    // // Check in same position
    REQUIRE(std::abs(gnorm(b_undo - a) - fly::gnorm(a - b)) < 0.001);
  }
}

TEST_CASE("TriGrid::gen_image", "[system]") {
  using namespace fly;

  Mat basis = Mat::Constant(10).triangularView<Eigen::Upper>();

  system::Triclinic box{basis, Arr<bool>::Constant(true)};

  auto grid = box.make_grid(3);

  {
    std::optional im = grid.gen_image<Sign::plus>(basis * Vec::Constant(0.5), 0);

    REQUIRE(!im);
  }

  {
    Vec origin = basis * Vec::Constant(0.1);

    for (int i = 0; i < spatial_dims; i++) {
      REQUIRE(!grid.gen_image<Sign::minus>(origin, i));

      auto im = grid.gen_image<Sign::plus>(origin, i);

      REQUIRE(im);

      Vec correct = origin + basis.col(i);

      REQUIRE(gnorm(correct - *im) < 0.001);
    }
  }

  {
    Vec origin = basis * Vec::Constant(0.9);

    for (int i = 0; i < spatial_dims; i++) {
      REQUIRE(!grid.gen_image<Sign::plus>(origin, i));

      auto im = grid.gen_image<Sign::minus>(origin, i);

      REQUIRE(im);

      Vec correct = origin - basis.col(i);

      REQUIRE(gnorm(correct - *im) < 0.001);
    }
  }
}