// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "libfly/utility/asserts.hpp"

#include <fmt/core.h>

#include <exception>
#include <stdexcept>
#include <string_view>

namespace fly::detail {

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

}  // namespace fly::detail