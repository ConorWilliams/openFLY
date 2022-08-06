#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
