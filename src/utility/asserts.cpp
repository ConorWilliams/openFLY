// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: MPL-2.0

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "libfly/utility/asserts.hpp"

#include <fmt/core.h>

#include <exception>
#include <stdexcept>
#include <string_view>

namespace otf::detail {

  struct libatom_unrecoverable : std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  [[noreturn]] void assert_handler(std::string_view expr, std::string_view msg, std::string_view file, long line,
                                   std::string_view func) {
    //
    fmt::print(stderr, "In {}:{}\n", file, line);
    fmt::print(stderr, "From function: {}\n", func);
    fmt::print(stderr, "Assertion \"{}\" failied!\n", expr);
    fmt::print(stderr, "With message: {}\n", msg);

    throw libatom_unrecoverable("libatom encountered an unrecoverable error");
  }

}  // namespace otf::detail