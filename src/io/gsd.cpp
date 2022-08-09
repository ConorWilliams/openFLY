// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/io/gsd.hpp"

#include <fmt/core.h>

#include <cstring>
#include <stdexcept>
#include <string_view>

#include "external/gsd.h"

namespace fly::io {

  /**
   * @brief Process and throw a GSD error
   */
  void throw_err(gsd_error retval, std::string_view fname) {
    switch (retval) {
      case GSD_SUCCESS:
        return;
      case GSD_ERROR_IO:
        throw std::runtime_error(fmt::format("GSD: {} - {}", std::strerror(errno), fname));
      case GSD_ERROR_INVALID_ARGUMENT:
        throw std::runtime_error(fmt::format("GSD: Invalid argument - {}", fname));
      case GSD_ERROR_NOT_A_GSD_FILE:
        throw std::runtime_error(fmt::format("GSD: Not a GSD file - {}", fname));
      case GSD_ERROR_INVALID_GSD_FILE_VERSION:
        throw std::runtime_error(fmt::format("GSD: Invalid GSD file version - {}", fname));
      case GSD_ERROR_FILE_CORRUPT:
        throw std::runtime_error(fmt::format("GSD: File corrupt - {}", fname));
      case GSD_ERROR_MEMORY_ALLOCATION_FAILED:
        throw std::runtime_error(fmt::format("GSD: Memory allocation failed - {}", fname));
      case GSD_ERROR_NAMELIST_FULL:
        throw std::runtime_error(fmt::format("GSD: Namelist full - {}", fname));
      case GSD_ERROR_FILE_MUST_BE_WRITABLE:
        throw std::runtime_error(fmt::format("GSD: File must be writeable - {}", fname));
      case GSD_ERROR_FILE_MUST_BE_READABLE:
        throw std::runtime_error(fmt::format("GSD: File must be readable - {}", fname));
      default:
        throw std::runtime_error(fmt::format("GSD: Unknown error - {}", fname));
    }
  }

  void test() {
    gsd_handle file_handle;

    gsd_create_and_open(&file_handle, "test.gsd", "openFLY", "hoomd", gsd_make_version(1, 4), GSD_OPEN_READWRITE, 0);
  }

}  // namespace fly::io