// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/io/gsd.hpp"

#include <catch2/catch_test_macros.hpp>

#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"

TEST_CASE("FileGSD", "[io]") {
  //
  using namespace fly;

  {  // Write

    // io::FileGSD file("FileGSD_test.gsd", io::create);

    // system::Box box(Mat::Identity(), Arr<bool>::Constant(true));

    // system::TypeMap<Index> map(2);

    // map.set(0, "Fe", 6);
    // map.set(1, "H", 1);

    // CHECK(map.get(0, tp_) == "Fe");
    // // CHECK(map.get(0, d_) == 6);

    // CHECK(map.get(1, tp_) == "H");
    // // CHECK(map.get(1, d_) == 1);

    // system::SoA<Position> atom(4);

    // atom(r_, 0) = Vec{0, 0, 0};
    // atom(r_, 1) = Vec{1, 0, 0};
    // atom(r_, 2) = Vec{0, 1, 0};
    // atom(r_, 3) = Vec{0, 0, 1};

    // file.commit([&] {
    //   file.write(box);
    //   file.write() file.write(r_, atom);
    // });
  }
}
