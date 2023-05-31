// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/typemap.hpp"

#include <catch2/catch.hpp>

#include "libfly/io/gsd.hpp"
#include "libfly/system/property.hpp"

TEST_CASE("TypeMap", "[system]") {
  //
  using namespace fly;

  system::TypeMap<Index> map(2);

  CHECK(map.num_types() == 2);

  map.set(1, i_, 99u);
  CHECK(map.get(1, i_) == 99u);

  map.set(1, tp_, Type::matrix_t::Zero());
  CHECK((map.get<Type>(1, tp_) == 0).all());

  map.set(1, tp_, "Woooop");
  CHECK(map.get(1, tp_) == "Woooop");

  map.set(0, "Fe", 2u);
  CHECK(map.get(0, tp_) == "Fe");
  CHECK(map.get(0, i_) == 2u);
}
