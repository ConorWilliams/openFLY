#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

//

#include <Eigen/Core>
#include <Eigen/Dense>

//

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

/**
 * \file core.hpp
 *
 * @brief Miscellaneous utilities.
 *
 * \rst
 * .. todo::
 *    refactor to trailing return type syntax.
 * \endrst
 */

/**
 * @brief Specify the number of spatial dimensions at compile time, defaults
 * to 3.
 *
 * \rst
 *
 * .. _`configure FLY_SPATIAL_DIMS`:
 *
 * To customize, at the compilation :ref:`configuration <Compiling openFLY>`
 * step append:
 *
 * .. code:: console
 *
 *     -DCMAKE_CXX_FLAGS='-DFLY_SPATIAL_DIMS=<number>'
 *
 * replacing ``<number>`` with the number of dimensions.
 *
 * \endrst
 */
#ifndef FLY_SPATIAL_DIMS
#  define FLY_SPATIAL_DIMS 3
#endif  // !FLY_SPATIAL_DIMS

namespace fly {

  // ------------------- Error Handling ---------------- //

  /**
   * @brief libFLY's catchable error type.
   */
  struct RuntimeError : std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * @brief Utility to make a RuntimeError object using the `{fmt}` library to
   * format the error message.
   *
   * @param fmt Format string.
   * @param args Arguments to forward to format string.
   * @return RuntimeError containing the formatted error message.
   */
  template <typename... Args>
  RuntimeError error(fmt::format_string<Args...> fmt, Args &&...args) {
    try {
      return RuntimeError{fmt::format(fmt, std::forward<Args>(args)...)};
    } catch (fmt::format_error const &err) {
      return RuntimeError{fmt::format("Failed to format '{}' with err: {}", fmt, err.what())};
    } catch (...) {
      return RuntimeError("Error during error handling");
    }
  }

  /**
   * @brief Utility to check condition and throw RuntimeError if condition is
   * false.
   *
   * Forwards ``fmt`` and ``args`` to fly::error().
   */
  template <typename... Args>
  void verify(bool condition, fmt::format_string<Args...> fmt, Args &&...args) {
    if (!condition) {
      throw error(std::move(fmt), std::forward<Args>(args)...);
    }
  }

#ifndef NDEBUG

  namespace detail {
    constexpr std::string_view file_name(std::string_view path) {
      if (auto k = path.find_last_of("/\\"); k != std::string_view::npos) {
        path.remove_prefix(k);
      }
      return path;
    }
  }  // namespace detail

/**
 * @brief Use like fly::verify but disabled if NDEBUG defined..
 */
#  define ASSERT(expr, string_literal, ...)                                                                          \
    do {                                                                                                             \
      if (constexpr std::string_view fname = fly::detail::file_name(__FILE__); !(expr)) {                            \
        throw fly::error("ASSERT \"{}\" failed in ...{}:{} | " string_literal, #expr, fname, __LINE__, __VA_ARGS__); \
      }                                                                                                              \
    } while (false)

#else

/**
 * @brief Use like std \c assert(expr) but with an error message.
 */
#  define ASSERT(...) \
    do {              \
    } while (false)

#endif  // !NDEBUG

  // ------------------- Defines, variables, etc. ---------------- //

  /**
   * @brief The number of spatial dimensions.
   *
   * \rst
   *
   * Configurable using the :ref:`FLY_SPATIAL_DIMS <configure FLY_SPATIAL_DIMS>`
   * macro.
   *
   * \endrst
   */
  inline constexpr int spatial_dims = FLY_SPATIAL_DIMS;

  static_assert(spatial_dims >= 2, "libFLY is not optimal for 1D simulations");

  /**
   * @brief Shorthand for creating an \c Eigen::Vector of doubles length \c
   * spatial_dims
   */
  using Vec = Eigen::Vector<double, spatial_dims>;

  /**
   * @brief Shorthand for creating an \c spatial_dims x \c spatial_dims \c
   * Eigen::Matrix of doubles.
   */
  using Mat = Eigen::Matrix<double, spatial_dims, spatial_dims>;

  /**
   * @brief Shorthand for creating an \c spatial_dims x 1 \c Eigen::Array.
   *
   * @tparam T The scalar type of the \c Eigen::Array
   */
  template <typename T>
  using Arr = Eigen::Array<T, spatial_dims, 1>;

  /**
   * @brief Strongly typed +/- sign.
   */
  enum class Sign : int {
    plus = 1,
    minus = -1,
  };

  namespace detail {

    template <typename T, typename...>
    struct First {
      using type = T;
    };

  }  // namespace detail

  /**
   * @brief Extracts the first type from a parameter pack.
   */
  template <typename... T>
  using first_t = typename detail::First<T...>::type;

  /**
   * @brief Non-deducable \c false for use in \c static_assert.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/utility/static_assert_false.cpp
   *    :code:
   *
   * \endrst
   */
  template <typename...>
  inline constexpr bool always_false = false;

  /**
   * @brief Strip all reference and const qualifications from \c T.
   *
   * @tparam T The type to strip ref/const qualifications from.
   */
  template <typename T>
  using remove_cref_t = std::remove_const_t<std::remove_reference_t<T>>;

  // ------------------- Small functions ---------------- //

  namespace detail {
    // C++20 functions see: https://en.cppreference.com/w/cpp/utility/intcmp
    template <class T, class U>
    constexpr bool cmp_less(T t, U u) noexcept {
      using UT = std::make_unsigned_t<T>;
      using UU = std::make_unsigned_t<U>;
      if constexpr (std::is_signed_v<T> == std::is_signed_v<U>)
        return t < u;
      else if constexpr (std::is_signed_v<T>)
        return t < 0 ? true : UT(t) < u;
      else
        return u < 0 ? false : t < UU(u);
    }

    template <class T, class U>
    constexpr bool cmp_greater(T t, U u) noexcept {
      return cmp_less(u, t);
    }

    template <class T, class U>
    constexpr bool cmp_less_equal(T t, U u) noexcept {
      return !cmp_greater(t, u);
    }

    template <class T, class U>
    constexpr bool cmp_greater_equal(T t, U u) noexcept {
      return !cmp_less(t, u);
    }

  }  // namespace detail

  /**
   * @brief Cast integral types asserting that conversion is lossless.
   *
   * Perform a `static_cast` from type `T` to `R` with bounds checking in debug
   * builds.
   *
   * Only SFINE enabled for integral types.
   *
   * @tparam R Target type to cast to.
   */
  template <typename R, typename T>
  constexpr auto safe_cast(T x) -> std::enable_if_t<std::is_integral_v<R> && std::is_integral_v<T>, R> {
    //
    static_assert(std::numeric_limits<unsigned int>::min() == 0);

    auto constexpr R_max = std::numeric_limits<R>::max();
    auto constexpr R_min = std::numeric_limits<R>::min();

    auto constexpr T_max = std::numeric_limits<T>::max();
    auto constexpr T_min = std::numeric_limits<T>::min();

    if constexpr (detail::cmp_less(R_max, T_max)) {
      ASSERT(detail::cmp_less_equal(x, R_max), "Could not cast '{}' to type R with R_max={}", x, R_max);
    }

    if constexpr (detail::cmp_greater(R_min, T_min)) {
      ASSERT(detail::cmp_greater_equal(x, R_min), "Could not cast '{}' to type R with R_min={}", x, R_min);
    }

    return static_cast<R>(x);
  }

  // ------------------- Math functions ---------------- //

  /**
   * @brief Test if two floating point numbers are within 0.01% of each other.
   *
   * Only SFINE enabled if T is floating point.
   */
  template <typename T>
  constexpr auto near(T a, T b) -> std::enable_if_t<std::is_floating_point_v<T>, bool> {
    return std::abs(a - b) <= 0.0001 * std::max(std::abs(a), std::abs(b));
  }

  /**
   * @brief Compute the shifted prefix product of an array.
   *
   * \rst
   *
   * If the inputs are :math:`x_i` the outputs are:
   *
   * .. math::
   *
   *    y_1 & = 1 \\
   *    y_2 & = 1 \times x_1 \\
   *    y_3 & = 1 \times x_1 \times x_2 \\
   *        & \vdots \\
   *    y_n & = 1 \times x_1\times x_{2} \times \dots \times x_{n-1}
   *
   * This is a common precursor operation when computing cell indices.
   *
   * \endrst
   *
   */
  template <typename T, int N>
  Eigen::Array<T, N, 1> product_scan(Eigen::Array<T, N, 1> x) {
    T prod = 1;
    for (Eigen::Index i = 0; i < x.size(); i++) {
      prod *= std::exchange(x[i], prod);
    }
    return x;
  }

  /**
   * @brief Compute integer powers of arithmetic types at compile time.
   *
   * Only SFINE enabled if base is an arithmetic type.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \text{x}^{\text{Exp}}
   *
   * \endrst
   */
  template <std::size_t Exp, typename T>
  constexpr auto ipow(T x) -> std::enable_if_t<std::is_arithmetic_v<T>, T> {
    if constexpr (Exp == 0) {
      return T(1);
    }
    if constexpr (Exp == 1) {
      return x;
    }
    if constexpr (Exp % 2 == 0) {
      return ipow<Exp / 2>(x) * ipow<Exp / 2>(x);
    } else {
      return ipow<Exp - 1>(x) * x;
    }
  }

  /**
   * @brief Generalised dot-product between two \c Eigen::Array objects.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \sum_{i} \sum_{j} a_{ij} b_{ij}
   *
   * \endrst
   */
  template <typename E1, typename E2>
  constexpr auto gdot(E1 const &a, E2 const &b) {
    return (a.array() * b.array()).sum();
  }

  /**
   * @brief Generic squared Frobenius norm of an \c Eigen::Array.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \sum_{i} \sum_{j} r^2_{ij}
   *
   * \endrst
   */
  template <typename E>
  constexpr auto gnorm_sq(E const &r) {
    return (r.array() * r.array()).sum();
  }

  /**
   * @brief Generic Frobenius norm of an \c Eigen::Array.
   *
   * \rst
   *
   * Computes:
   *
   * .. math::
   *
   *    \sqrt{\sum_{i} \sum_{j} r^2_{ij}}
   *
   * \endrst
   */
  template <typename E>
  constexpr auto gnorm(E &&expr) {
    return std::sqrt(gnorm_sq(std::forward<E>(expr)));
  }

  /**
   * @brief Compute the unit normal of a hyperplane through a set of points.
   *
   * \rst
   *
   * Only SFINE enabled for fixed-size square matrices.
   *
   * See: `StackExchange
   * <https://math.stackexchange.com/questions/2301110/fastest-way-to-find-equation-of-hyperplane>`_
   *
   * \endrst
   *
   * @param P An NxN fixed-size matrix.
   * @return Eigen::Vector<Scalar, N> The unit normal of the hyperplane passing
   * through the column vectors of ``P``.
   */
  template <typename Scalar, int N>
  auto hyperplane_normal(Eigen::Matrix<Scalar, N, N> const &P) -> std::enable_if_t<N != Eigen::Dynamic, Eigen::Vector<Scalar, N>> {
    //
    Eigen::Matrix<Scalar, N, N + 1> H = Eigen::Matrix<Scalar, N, N + 1>::Ones();

    H.template topLeftCorner<N, N>() = P.transpose();

    Eigen::FullPivLU<decltype(H)> lu(H);

    if (auto kd = lu.dimensionOfKernel(); kd != 1) {
      throw error("Points passed to hyperplane_normal are linearly dependant");
    }

    ASSERT(near(lu.kernel()(N, 0), 0.0), "Homogeneous coordinate of the kernel = {} but should be zero!", lu.kernel()(N, 0));

    return lu.kernel().template topLeftCorner<N, 1>().normalized();
  }

  // ------------------- Classes ---------------- //

  /**
   * @brief Basic implementation of a Golang like defer.
   *
   * Only SFINE enabled if callable is noexcept.
   *
   * \tparam F The invocable's type, this **MUST** be deducted through CTAD by the
   * deduction guide.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/utility/defer.cpp
   *    :code:
   *
   * \endrst
   */
  template <class F, typename = std::enable_if_t<std::is_nothrow_invocable_v<F &&>>>
  class [[nodiscard]] Defer {
  public:
    /**
     * @brief Construct a new Defer object.
     *
     * @param f Invocable forwarded into object and invoked by destructor.
     */
    constexpr Defer(F &&f) : m_f(std::forward<F>(f)) {}

    Defer(const Defer &) = delete;
    Defer(Defer &&other) = delete;
    Defer &operator=(const Defer &) = delete;
    Defer &operator=(Defer &&) = delete;

    /**
     * @brief Call the invocable.
     */
    ~Defer() noexcept { std::invoke(std::forward<F>(m_f)); }

  private:
    F m_f;
  };

  /**
   * @brief Forwarding deduction guide.
   */
  template <typename F>
  Defer(F &&) -> Defer<F>;

  // ------------------- Timing ---------------- //

  /**
   * @brief Transparent function wrapper that measures the execution time of a
   * function.
   *
   * The execution time is printed to stdout. Garantees RVO.
   *
   * @param name Name of function being called, also printed to stdout.
   * @param f Function call.
   * @param args Arguments to call \c f with.
   * @return std::invoke_result_t<F&&, Args&&...> The result of calling \c f with
   * \c args... .
   */
  template <typename F, typename... Args>
  std::invoke_result_t<F &&, Args &&...> timeit(std::string_view name, F &&f, Args &&...args) {
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

    if constexpr (std::is_void_v<std::invoke_result_t<F &&, Args &&...>>) {
      std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    }
  }

}  // namespace fly
