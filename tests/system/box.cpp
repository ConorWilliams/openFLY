// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/system/box.hpp"

#include <catch2/catch_test_macros.hpp>

#include "libfly/utility/core.hpp"

TEST_CASE("Box::Box", "[system]") {
  //
  using namespace fly;

  {  //
    system::Box box(Mat::Identity(), Arr<bool>::Constant(true));

    CHECK(box.holding<system::Orthorhombic>());
  }

  {  //
    system::Box box(Mat::Ones().triangularView<Eigen::Upper>(), Arr<bool>::Constant(true));

    CHECK(box.holding<system::Triclinic>());
  }
}
