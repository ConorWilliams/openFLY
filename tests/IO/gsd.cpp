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
#include <cstdint>
#include <cstdlib>
#include <string_view>

#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

//

using namespace fly;

TEST_CASE("BinaryFile", "[io]") {
  //

  system::Box box(Mat::Identity(), Arr<bool>::Constant(true));

  system::TypeMap<Index> map(2);

  map.set(0, "Fe", 6);
  map.set(1, "H", 1);

  system::Supercell cell = system::make_supercell<Position>(box, map, 4);

  cell(r_, 0) = Vec::Constant(0);

  cell(r_, 1) = Vec::Constant(1);
  cell(r_, 2) = Vec::Constant(2);
  cell(r_, 3) = Vec::Constant(3);

  cell(id_, 0) = 0;

  cell(id_, 1) = 1;
  cell(id_, 2) = 1;
  cell(id_, 3) = 1;

  Position::array_t initial = cell[r_];

  {  // Write

    io::BinaryFile file("BinaryFile_test.gsd", io::create);

    file.commit([&] {
      file.write(cell.box());
      file.write(cell.map());
      file.write("particles/N", safe_cast<std::uint32_t>(cell.size()));
      file.write(id_, cell);

      file.write(r_, cell);
    });

    cell[r_] += 1;

    file.commit([&] {
      file.write("particles/N", safe_cast<std::uint32_t>(cell.size()));
      file.write(r_, cell);
    });
  }

  {  // Read
    io::BinaryFile file("BinaryFile_test.gsd", io::read_only);

    CHECK(file.n_frames() == 2);

    system::TypeMap out_map = file.read_map<Index>(0);

    CHECK(out_map.num_types() == map.num_types());

    for (uint32_t i = 0; i < map.num_types(); i++) {
      CHECK(out_map.get(i, Type{}) == map.get(i, Type{}));
      CHECK(out_map.get(i, Index{}) == map.get(i, Index{}));
    }

    system::Box out_box = file.read_box(0);

    Mat diff = out_box.basis() - box.basis();

    CHECK(gnorm(diff) < 0.0001);

    auto out_cell = system::make_supercell<Position>(out_box, out_map, file.read<std::uint32_t>(0, "particles/N"));

    file.read_to(0, id_, out_cell);
    file.read_to(0, r_, out_cell);

    CHECK((out_cell[id_] == cell[id_]).all());
    CHECK((out_cell[r_] == initial).all());

    file.read_to(1, r_, out_cell);

    CHECK((out_cell[r_] == out_cell[r_]).all());
  }
}
