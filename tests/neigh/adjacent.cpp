// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/neigh/adjacent.hpp"

#include <catch2/catch_test_macros.hpp>

#include "libfly/system/atom.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("AdjacentCells", "[neigh]") {
  //
  if constexpr (fly::spatial_dims == 3) {
    //
    fly::neigh::detail::AdjacentCells cells({5, 5, 5});

    {
      auto neigh = cells[(1) + (1) * 5 + (1) * 5 * 5];

      REQUIRE(neigh.size() == 26);

      REQUIRE(neigh[0] == (1 - 1) + (1 - 1) * 5 + (1 - 1) * 5 * 5);
      REQUIRE(neigh[1] == (1 - 0) + (1 - 1) * 5 + (1 - 1) * 5 * 5);
      REQUIRE(neigh[2] == (1 + 1) + (1 - 1) * 5 + (1 - 1) * 5 * 5);

      REQUIRE(neigh[3] == (1 - 1) + (1 - 0) * 5 + (1 - 1) * 5 * 5);
      REQUIRE(neigh[4] == (1 - 0) + (1 - 0) * 5 + (1 - 1) * 5 * 5);
      REQUIRE(neigh[5] == (1 + 1) + (1 - 0) * 5 + (1 - 1) * 5 * 5);

      REQUIRE(neigh[6] == (1 - 1) + (1 + 1) * 5 + (1 - 1) * 5 * 5);
      REQUIRE(neigh[7] == (1 - 0) + (1 + 1) * 5 + (1 - 1) * 5 * 5);
      REQUIRE(neigh[8] == (1 + 1) + (1 + 1) * 5 + (1 - 1) * 5 * 5);

      REQUIRE(neigh[9] == (1 - 1) + (1 - 1) * 5 + (1 - 0) * 5 * 5);
      REQUIRE(neigh[10] == (1 - 0) + (1 - 1) * 5 + (1 - 0) * 5 * 5);
      REQUIRE(neigh[11] == (1 + 1) + (1 - 1) * 5 + (1 - 0) * 5 * 5);

      REQUIRE(neigh[12] == (1 - 1) + (1 - 0) * 5 + (1 - 0) * 5 * 5);
      // Does not include the centre cell
      REQUIRE(neigh[13] == (1 + 1) + (1 - 0) * 5 + (1 - 0) * 5 * 5);

      REQUIRE(neigh[14] == (1 - 1) + (1 + 1) * 5 + (1 - 0) * 5 * 5);
      REQUIRE(neigh[15] == (1 - 0) + (1 + 1) * 5 + (1 - 0) * 5 * 5);
      REQUIRE(neigh[16] == (1 + 1) + (1 + 1) * 5 + (1 - 0) * 5 * 5);

      REQUIRE(neigh[17] == (1 - 1) + (1 - 1) * 5 + (1 + 1) * 5 * 5);
      REQUIRE(neigh[18] == (1 - 0) + (1 - 1) * 5 + (1 + 1) * 5 * 5);
      REQUIRE(neigh[19] == (1 + 1) + (1 - 1) * 5 + (1 + 1) * 5 * 5);

      REQUIRE(neigh[20] == (1 - 1) + (1 - 0) * 5 + (1 + 1) * 5 * 5);
      REQUIRE(neigh[21] == (1 - 0) + (1 - 0) * 5 + (1 + 1) * 5 * 5);
      REQUIRE(neigh[22] == (1 + 1) + (1 - 0) * 5 + (1 + 1) * 5 * 5);

      REQUIRE(neigh[23] == (1 - 1) + (1 + 1) * 5 + (1 + 1) * 5 * 5);
      REQUIRE(neigh[24] == (1 - 0) + (1 + 1) * 5 + (1 + 1) * 5 * 5);
      REQUIRE(neigh[25] == (1 + 1) + (1 + 1) * 5 + (1 + 1) * 5 * 5);
    }

    // Test face cell

    {
      auto neigh = cells[((2) + (2) * 5 + (0) * 5 * 5)];

      REQUIRE(neigh.size() == 9 + 8);
    }

    // Test Edge cell

    {
      auto neigh = cells[((2) + (0) * 5 + (0) * 5 * 5)];

      REQUIRE(neigh.size() == 11);
    }

    // Test corner cell

    {
      auto neigh = cells[((0) + (0) * 5 + (0) * 5 * 5)];

      REQUIRE(neigh.size() == 7);
    }
  }
}