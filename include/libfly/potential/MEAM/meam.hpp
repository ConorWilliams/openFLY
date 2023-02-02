#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <memory>

#include "KIM/meam_c.hpp"
#include "KIM/meam_spline.hpp"
#include "KIM/meam_sw_spline.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file meam.hpp
 *
 * @brief MEAM potential implementation.
 */

namespace fly::potential {

  /**
   * @brief MEAM potential child class.
   */
  class MEAM {
  public:
    /**
     * @brief Construct a new MEAM object.
     *
     * @param map TypeMap of expected types to find in data.
     */
    MEAM(system::TypeMap<> const& map, std::string const&) {}

  private:
  };

}  // namespace fly::potential