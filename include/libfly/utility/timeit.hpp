#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <fmt/chrono.h>
#include <fmt/core.h>

#include <chrono>
#include <functional>
#include <string_view>
#include <type_traits>
#include <utility>

#include "libfly/utility/core.hpp"

/**
 *
 * \file timeit.hpp
 *
 * @brief Timing utilities.
 *
 */

namespace fly {

  /**
   * @brief Transparent function wrapper that measures the execution time of a function.
   *
   * The execution time is printed to stdout. Garantees RVO.
   *
   * @param name Name of function being called, also printed to stdout.
   * @param f Function call.
   * @param args Arguments to call \c f with.
   * @return std::invoke_result_t<F&&, Args&&...> The result of calling \c f with \c args... .
   */
  template <typename F, typename... Args> std::invoke_result_t<F&&, Args&&...> timeit(std::string_view name, F&& f, Args&&... args) {
    //
    auto start = std::chrono::steady_clock::now();

    Defer _ = [&]() noexcept {
      //
      using namespace std::chrono;

      auto elapsed = steady_clock::now() - start;

      auto sec = duration_cast<seconds>(elapsed);

      elapsed -= sec;

      auto mil = duration_cast<milliseconds>(elapsed);

      elapsed -= mil;

      auto mic = duration_cast<microseconds>(elapsed);

      elapsed -= mic;

      auto nan = duration_cast<nanoseconds>(elapsed);

      fmt::print("Timing \"{}\" {:>4} {:>5} {:>5} {:>5}\n", name, sec, mil, mic, nan);
    };

    if constexpr (std::is_void_v<std::invoke_result_t<F&&, Args&&...>>) {
      std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    }
  }

}  // namespace fly
