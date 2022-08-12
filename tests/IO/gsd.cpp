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

    system::Box box(Mat::Identity(), Arr<bool>::Constant(true));

    system::TypeMap<Index> map(2);

    map.set(0, "Fe", 6);
    map.set(1, "H", 1);

    system::Supercell cell = system::make_supercell<Position>(box, map, 4);

    cell(r_, 0) = Vec{0, 0, 0};

    cell(r_, 1) = Vec{1, 0, 0};
    cell(r_, 2) = Vec{0, 1, 0};
    cell(r_, 3) = Vec{0, 0, 1};

    cell(id_, 0) = 0;

    cell(id_, 1) = 1;
    cell(id_, 2) = 1;
    cell(id_, 3) = 1;


  {  // Write

    io::FileGSD file("FileGSD_test.gsd", io::create);

    file.commit([&] {
      file.write(cell.box());
      file.write(cell.map());

      file.write(id_, cell);
      file.write(r_, cell);
    });
    
  }
}
