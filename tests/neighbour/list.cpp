// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/neighbour/list.hpp"

#include <omp.h>

#include <catch2/catch_test_macros.hpp>
#include <random>

#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

using namespace fly;

TEST_CASE("List", "[neighbour]") {
  //
  REQUIRE(true);
  //
}

// struct Neigh {
//   std::size_t i;
//   Vec3<floating> dr;
// };

// void slow_neigh_list(std::vector<std::vector<Neigh>>& nl, SimCell const& atoms, floating r_cut) {
//   nl.resize(atoms.size());

//   for (std::size_t i = 0; i < atoms.size(); i++) {
//     nl[i].clear();

//     for (std::size_t j = 0; j < atoms.size(); j++) {
//       if (i != j) {
//         Vec3<floating> dr = atoms.min_image(atoms(Position{}, i), atoms(Position{}, j));

//         if (norm(dr) < r_cut) {
//           nl[i].push_back({j, dr});
//         }
//       }
//     }
//     std::sort(nl[i].begin(), nl[i].end(), [](auto const& a, auto const& b) { return a.i < b.i; });
//   }
// }

// TEST_CASE("Neighbour list speed testing") {
//   SimCell atoms({{1, 2, 1}, {true, true, true}});

//   randomise(atoms, 10'000 * atoms.extents().prod());

//   fmt::print("num atoms is {}\n", atoms.size());

//   floating r_cut = 0.1;

//   {
//     neighbour::List neigh;

//     neigh.rebuild(atoms, r_cut);  // Warm up + alloc

//     timeit("Fast", [&] { neigh.rebuild(atoms, r_cut); });

//     int x = 0;
//     int y = 0;

//     timeit("Counting", [&] {
//       for (std::size_t i = 0; i < atoms.size(); i++) {
//         //
//         x++;
//         neigh.for_neighbours(i, [&](std::size_t, floating, Vec3<floating> const&) { y++; });
//       }
//     });

//     fmt::print("num neigh = {}\n", (double)y / x);
//   }

//   {
//     std::vector<std::vector<Neigh>> nl;

//     slow_neigh_list(nl, atoms, r_cut);  // Warm up + alloc

//     timeit("Slow", [&] { slow_neigh_list(nl, atoms, r_cut); });
//   }
// }

void test(system::Box const& box, system::SoA<Position const&> atoms, double r_cut, int num_threads) {
  //
  neighbour::List neigh(box, r_cut);

  timeit("rebuild", [&] { neigh.rebuild(atoms, num_threads); });

  //   neigh.rebuild(atoms, num_threads);

  //   std::vector<std::vector<Neigh>> nl;

  //   slow_neigh_list(nl, atoms, r_cut);

  //   for (std::size_t i = 0; i < atoms.size(); i++) {
  //     //
  //     std::vector<Neigh> nl2;

  //     neigh.for_neighbours(i, [&](std::size_t n, Vec3<floating> const& dr) { nl2.push_back({neigh.image_to_real(n), dr}); });

  //     std::sort(nl2.begin(), nl2.end(), [](auto const& a, auto const& b) { return a.i < b.i; });

  //     // Test same number of neighbours
  //     REQUIRE(nl2.size() == nl[i].size());

  //     for (std::size_t j = 0; j < nl2.size(); j++) {
  //       // Test all neighbours have the same index
  //       REQUIRE(nl2[j].i == nl[i][j].i);
  //       // Test all neighbours have the same minimum image positions
  //       REQUIRE(norm(nl2[j].dr - nl[i][j].dr) < 0.001);
  //     }
  //   }
}

TEST_CASE("List fuzz-testing", "[neighbour]") {
  //

  //
  fly::Xoshiro gen({1, 2, 3, 4});
  std::uniform_real_distribution<double> dis(0, 1);
  std::uniform_int_distribution<Eigen::Index> idis(100, 1000);

  auto rand = [&dis, &gen]() { return dis(gen); };

  for (int i = 0; i < 1; i++) {
    // Random simulation box
    fly::system::Box box = [&]() {
      //
      Mat basis = Mat::Ones() + 9 * Mat::NullaryExpr(rand);

      basis = basis.triangularView<Eigen::Upper>();

      if (true) {
        basis = basis.triangularView<Eigen::Lower>();
      }

      fly::system::Box tmp{basis, fly::Arr<bool>::NullaryExpr(rand)};

      if (true) {
        CHECK(tmp.holding<system::Orthorhombic>());
      }

      return tmp;
    }();

    system::SoA<Position> cell(idis(gen) * 20);

    for (int j = 0; j < cell.size(); j++) {
      cell(r_, j) = box.basis() * fly::Vec::NullaryExpr(rand);
    }

    for (std::size_t j = 0; j < 10; j++) {
      test(box, cell, 0.25 + 0.2499 * dis(gen), 1);
      //   test(box, cell, 0.25 + 0.2499 * dis(gen), omp_get_max_threads());
    }
  }
}