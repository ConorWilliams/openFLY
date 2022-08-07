

#include "libfly/neighbour/adjacent.hpp"

#include <catch2/catch_test_macros.hpp>

#include "libfly/system/atom.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("AdjacentCells", "[neighbour]") {
  //
  if constexpr (fly::spatial_dims == 3) {
    //
    fly::neighbour::AdjacentCells cells({5, 5, 5});

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