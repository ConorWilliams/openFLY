// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/hessian.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace fly;

TEST_CASE("Hessian", "[system]") {
  //
  system::Hessian H;

  int n = 2;

  H.zero_for(n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      REQUIRE(H(i, j) == Mat::Zero());
    }
  }

  auto values = H.eigenvalues();

  for (auto const& elem : values) {
    REQUIRE(near(elem, 0.0));
  }

  H.zero_for(1);

  H(0, 0) = Mat::Constant(1.0);

  REQUIRE(near(H.eigenvalues().array().sum(), static_cast<double>(spatial_dims)));
}