#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <string_view>

#include "libfly/utility/current_function.hpp"

/**
 * \file asserts.hpp
 *
 * @brief Defines a set of \<cassert\> like macros.
 */

namespace fly::detail {

  /**
   * @brief Kill program printing diagnostics + stacktrace, should not be inlined to reduce code size.
   */
  [[noreturn]] void assert_handler(std::string_view expr, std::string_view msg, std::string_view file, long line,
                                   std::string_view func);

}  // namespace fly::detail

/**
 * @brief Use like std \c assert(expr) but with an error message, not disabled by \c NDEBUG.
 */
#define VERIFY(expr, msg)                                                                \
  do {                                                                                   \
    if (!(expr)) {                                                                       \
      fly::detail::assert_handler(#expr, msg, __FILE__, __LINE__, FLY_CURRENT_FUNCTION); \
    }                                                                                    \
  } while (false)

#ifndef NDEBUG

/**
 * @brief Use like std \c assert(expr) but with an error message.
 */
#  define ASSERT(expr, msg) VERIFY(expr, msg)

#else

/**
 * @brief Use like std \c assert(expr) but with an error message.
 */
#  define ASSERT(expr, msg) \
    do {                    \
    } while (false)

#endif  // !NDEBUG
