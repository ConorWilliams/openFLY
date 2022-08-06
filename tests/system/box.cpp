// clang-format off
#include <fmt/core.h>
#include "libfly/system/box.hpp"

#include <catch2/catch_test_macros.hpp>

#include "libfly/system/boxes/orthorhombic.hpp"

TEST_CASE("AdjacentCells", "[system]") {
  //
  fly::system::Orthorhombic box{{10, 10, 10}, {true, true, true}};

  auto grid = box.make_grid(3);

  using V = fly::Vec<fly::floating>;

  REQUIRE(grid.cell_idx(box.canon_image(V{0, 0, 0}) + grid.cell_offset()) == 1 + 1 * 5 + 1 * 5 * 5);
  REQUIRE(grid.cell_idx(box.canon_image(V{0, 0, 0}) + grid.cell_offset()) == 1 + 1 * 5 + 1 * 5 * 5);

  REQUIRE(grid.cell_idx(box.canon_image(V{5, 5, 5}) + grid.cell_offset()) == 2 + 2 * 5 + 2 * 5 * 5);

  REQUIRE(grid.cell_idx(box.canon_image(V{9.999, 9.999, 9.999}) + grid.cell_offset()) == 3 + 3 * 5 + 3 * 5 * 5);

  fly::system::AdjacentCells cells(grid.shape());

  //   // Test inner cell

 

  {
    auto neigh = cells[(1) + (1) * 5 + (1) * 5 * 5];

    for(auto elem : neigh){
     fmt::print("{}\n", elem);
}
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