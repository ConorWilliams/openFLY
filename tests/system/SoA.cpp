// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/SoA.hpp"

#include <catch2/catch_test_macros.hpp>
#include <utility>

#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("SoA", "[system]") {
  using namespace fly;
  using namespace system;

  SoA<Index, Position, PotentialGradient> arr;

  CHECK(arr.size() == 0);

  int n = 10;

  arr.destructive_resize(n);

  CHECK(arr.size() == n);

  for (int i = 0; i < n; i++) {
    arr(i_, i) = safe_cast<uint32_t>(i);
    arr(r_, i) = Position::matrix_t::Zero();
    arr(r_, i)[0] = i;
  }

  {
    unsigned int count = 0;

    for (auto&& elem : arr[i_]) {
      CHECK(elem == count++);
    }

    count = 0;

    for (auto&& elem : std::as_const(arr)[i_]) {
      CHECK(elem == count++);
    }

    count = 0;

    for (auto&& elem : arr[r_]) {
      if (count % spatial_dims == 0) {
        CHECK(near<double>(elem, count / spatial_dims));
      } else {
        CHECK(near<double>(elem, 0));
      }

      count++;
    }

    count = 0;

    for (auto&& elem : std::as_const(arr)[r_]) {
      if (count % spatial_dims == 0) {
        CHECK(near<double>(elem, count / spatial_dims));
      } else {
        CHECK(near<double>(elem, 0));
      }

      count++;
    }
  }

  arr[g_] = arr[r_];

  CHECK((arr[g_] == arr[r_]).all());

  //   /////////////////////////////////////

  SoA<Index&> view_p = arr;

  CHECK((view_p[i_] == arr[i_]).all());
}
