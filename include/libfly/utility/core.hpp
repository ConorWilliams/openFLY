#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <fmt/core.h>
#include <fmt/format.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "libfly/utility/asserts.hpp"  // API specifies we must include this

/**
 * \file core.hpp
 *
 * @brief Miscellaneous utilities.
 */

/**
 * @brief Specify the number of spatial dimensions at compile time, defaults to 3.
 *
 * \rst
 *
 * .. _`configure FLY_SPATIAL_DIMS`:
 *
 * To customize, at the compilation :ref:`configuration <Compiling openFLY>` step append:
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

  /**
   * @brief The number of spatial dimensions.
   *
   * \rst
   *
   * Configurable using the :ref:`FLY_SPATIAL_DIMS <configure FLY_SPATIAL_DIMS>` macro.
   *
   * \endrst
   */
  inline constexpr int spatial_dims = FLY_SPATIAL_DIMS;

  static_assert(spatial_dims >= 2, "libFLY is not optimal for 1D simulations");

  /**
   * @brief The maximum atomic number that any atom can have.
   */
  inline constexpr int max_atomic_num = 111;

  /**
   * @brief Shorthand for creating an \c Eigen::Vector of length \c spatial_dims
   *
   * @tparam T The scalar type of the \c Eigen::Vector
   */
  template <typename T>
  using Vec = Eigen::Vector<T, spatial_dims>;

  /**
   * @brief Shorthand for creating an \c spatial_dims x 1 \c Eigen::Array.
   *
   * @tparam T The scalar type of the \c Eigen::Array
   */
  template <typename T>
  using Arr = Eigen::Array<T, spatial_dims, 1>;

  /**
   * @brief Shorthand for creating an \c spatial_dims x \c spatial_dims \c Eigen::Matrix.
   *
   * @tparam T The scalar type of the \c Eigen::Matrix
   */
  template <typename T>
  using Mat = Eigen::Matrix<T, spatial_dims, spatial_dims>;

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
  template <typename... Ts>
  using first_t = typename detail::First<Ts...>::type;

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
   * @brief Bounds checked integer cast in debug builds.
   *
   * Perform a `static_cast` from type `T` to `R` with bounds checking in debug builds.
   *
   * @tparam R Target type to cast to.
   */
  template <typename R, typename T, typename = std::enable_if_t<std::is_integral_v<R> && std::is_integral_v<T>>>
  R safe_cast(T x) {
    //
    static_assert(std::numeric_limits<unsigned int>::min() == 0);

    auto constexpr R_max = std::numeric_limits<R>::max();
    auto constexpr R_min = std::numeric_limits<R>::min();

    auto constexpr T_max = std::numeric_limits<T>::max();
    auto constexpr T_min = std::numeric_limits<T>::min();

    if constexpr (detail::cmp_less(R_max, T_max)) {
      ASSERT(detail::cmp_less_equal(x, R_max), fmt::format("Could not cast '{}' to desired type", x));
    }

    if constexpr (detail::cmp_greater(R_min, T_min)) {
      ASSERT(detail::cmp_greater_equal(x, R_min), fmt::format("Could not cast '{}' to desired type", x));
    }

    return static_cast<R>(x);
  }

  /**
   * @brief libFLY's catchable error type.
   */
  struct RuntimeError : std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * @brief Utility to build a RuntimeError using `{fmt}` to build the error message.
   *
   * @param fmt Format string.
   * @param args Arguments to forward to format string.
   * @return RuntimeError containing the formatted error message.
   */
  template <typename... T, typename... Args>
  RuntimeError error(fmt::format_string<T...> fmt, Args&&... args) {
    return RuntimeError{fmt::format(fmt, std::forward<Args>(args)...)};
  }

  /**
   * @brief Force constant evaluation in C++17
   */
  template <auto Value>
  inline constexpr auto const_eval = Value;

  namespace detail {

    constexpr std::string_view file_name(std::string_view path) {
      if (auto k = path.find_last_of("/\\"); k != std::string_view::npos) {
        path.remove_prefix(k);
      }
      return path;
    }

  }  // namespace detail

#ifndef NDEBUG

/**
 * @brief Use like std \c assert(expr) but with a formatted error message.
 */
#  define XASSERT(expr, string_literal, ...)                                                                         \
    do {                                                                                                             \
      if (constexpr std::string_view fname = fly::detail::file_name(__FILE__); !(expr)) {                            \
        throw fly::error("ASSERT \"{}\" failed in ...{}:{} | " string_literal, #expr, fname, __LINE__, __VA_ARGS__); \
      }                                                                                                              \
    } while (false)

#else

/**
 * @brief Use like std \c assert(expr) but with an error message.
 */
#  define XASSERT(...) \
    do {               \
    } while (false)

#endif  // !NDEBUG

  /**
   * @brief Test if two floating point numbers are within 0.01% of each other.
   */
  template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
  constexpr bool near(T a, T b) {
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
    for (int i = 0; i < x.size(); i++) {
      prod *= std::exchange(x[i], prod);
    }
    return x;
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
  auto gdot(E1 const& a, E2 const& b) {
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
  auto gnorm_sq(E const& r) {
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
  auto gnorm(E&& expr) {
    return std::sqrt(gnorm_sq(std::forward<E>(expr)));
  }

  /**
   * @brief Compute the unit normal of a hyperplane through a set of points.
   *
   * \rst
   *
   * See: `StackExchange <https://math.stackexchange.com/questions/2301110/fastest-way-to-find-equation-of-hyperplane>`_
   *
   * \endrst
   *
   * @param P An NxN fixed-size matrix.
   * @return Eigen::Vector<Scalar, N> The unit normal of the hyperplane passing through the column vectors of ``P``.
   */
  template <typename Scalar, int N>
  std::enable_if_t<N != Eigen::Dynamic, Eigen::Vector<Scalar, N>> hyperplane_normal(Eigen::Matrix<Scalar, N, N> const& P) {
    //
    Eigen::Matrix<Scalar, N, N + 1> H = Eigen::Matrix<Scalar, N, N + 1>::Ones();

    H.template topLeftCorner<N, N>() = P.transpose();

    Eigen::FullPivLU<decltype(H)> lu(H);

    ASSERT(lu.dimensionOfKernel() == 1, "Points are linearly dependant!");

    ASSERT(near(lu.kernel()(N, 0), 0.0), "Homogeneous coordinate of the kernel should be zero!");

    return lu.kernel().template topLeftCorner<N, 1>().normalized();
  }

  /**
   * @brief Compute integer powers of arithmetic types at compile time.
   *
   * Only SFINE if base is an arithmetic type.
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
  constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> ipow(T x) {
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
   * @brief Basic implementation of a Golang like defer.
   *
   * Only SFINE if callable is noexcept.
   *
   * \tparam F The invocable's type, this **MUST** be deducted through CTAD by the deduction guide.
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
  template <class F, typename = std::enable_if_t<std::is_nothrow_invocable_v<F&&>>>
  class [[nodiscard]] Defer {
  public:
    /**
     * @brief Construct a new Defer object.
     *
     * @param f Invocable forwarded into object and invoked by destructor.
     */
    Defer(F&& f) : m_f(std::forward<F>(f)) {}

    Defer(const Defer&) = delete;
    Defer(Defer&& other) = delete;
    Defer& operator=(const Defer&) = delete;
    Defer& operator=(Defer&&) = delete;

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
  Defer(F&&) -> Defer<F>;

}  // namespace fly
