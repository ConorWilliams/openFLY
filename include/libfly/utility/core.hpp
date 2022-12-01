#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

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
#include <cereal/types/vector.hpp>
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
#include <variant>
#include <vector>

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
   * @param args Format arguments (forwarded to ``fmt::format``).
   * @return A fly::RuntimeError containing the formatted error message as by ``fmt::format(fmt, args...)``.
   */
  template <typename... Args>
  auto error(fmt::format_string<Args...> fmt, Args &&...args) -> RuntimeError {
    try {
      return RuntimeError{fmt::format(fmt, std::forward<Args>(args)...)};
    } catch (fmt::format_error const &err) {
      return RuntimeError{fmt::format("Failed to format '{}' with err: {}", fmt, err.what())};
    } catch (...) {
      return RuntimeError("Error during error handling");
    }
  }

  /**
   * @brief Utility to check ``cond`` and throw RuntimeError if condition is ``false``.
   *
   * \rst
   *
   * Shorthand for:
   *
   * .. code::
   *
   *     if (!cond) {
   *         throw fly::error(fmt, args...);
   *     }
   *
   * \endrst
   *
   * @param cond Condition to check.
   * @param fmt Format string.
   * @param args Format arguments (forwarded to ``fmt::format``).
   *
   * @return void.
   */
  template <typename... Args>
  constexpr auto verify(bool cond, fmt::format_string<Args...> fmt, Args &&...args) -> void {
    if (!cond) {
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
 * @brief Use like fly::verify() but disabled if ``NDEBUG`` defined.
 */
#  define ASSERT(expr, string_literal, ...)                                                                          \
    do {                                                                                                             \
      if (constexpr std::string_view fname = fly::detail::file_name(__FILE__); !(expr)) {                            \
        throw fly::error("ASSERT \"{}\" failed in ...{}:{} | " string_literal, #expr, fname, __LINE__, __VA_ARGS__); \
      }                                                                                                              \
    } while (false)

#else

/**
 * @brief Use like fly::verify() but disabled if ``NDEBUG`` defined.
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
   * Must be greater than or equal to 2. Configurable using the :ref:`FLY_SPATIAL_DIMS <configure FLY_SPATIAL_DIMS>` macro.
   *
   * \endrst
   */
  inline constexpr int spatial_dims = FLY_SPATIAL_DIMS;

  static_assert(spatial_dims >= 2, "libFLY is not optimal for 1D simulations");

  /**
   * @brief Shorthand for creating an \c Eigen::Vector of doubles of length \c spatial_dims.
   */
  using Vec = Eigen::Vector<double, spatial_dims>;

  /**
   * @brief Shorthand for creating an \c spatial_dims x \c spatial_dims \c Eigen::Matrix of doubles.
   */
  using Mat = Eigen::Matrix<double, spatial_dims, spatial_dims>;

  /**
   * @brief Shorthand for creating an \c spatial_dims x 1 \c Eigen::Array.
   *
   * @tparam T The scalar type of the \c Eigen::Array.
   */
  template <typename T>
  using Arr = Eigen::Array<T, spatial_dims, 1>;

  /**
   * @brief Strongly typed +/- sign.
   */
  enum class Sign : int {
    plus = 1,    ///< Representing +1 or the positive direction.
    minus = -1,  ///< Representing -1 or the negative direction.
  };

  // ------------------- Meta --------------------- //

  namespace detail {

    template <typename From, typename To, typename = void>
    struct is_narrowing_conversion_impl : std::true_type {};

    template <typename From, typename To>
    struct is_narrowing_conversion_impl<From, To, std::void_t<decltype(To{std::declval<From>()})>> : std::false_type {};

  }  // namespace detail

  /**
   * @brief Check if a numerical conversion is narrowing.
   *
   * Leverages brace initialisation, adapted from: https://stackoverflow.com/a/67603594
   *
   * @tparam From The source type.
   * @tparam To The target type.
   */
  template <typename From, typename To>
  inline constexpr bool is_narrowing_conversion_v = detail::is_narrowing_conversion_impl<From, To>::value;

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
   * @brief Non-deducible \c false for use in \c static_assert.
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

    // See https://en.cppreference.com/w/cpp/experimental/is_detected

    template <class Default, class AlwaysVoid, template <class...> class Fn, class... Args>
    struct detector : std::false_type {
      using type = Default;
    };

    template <class Default, template <class...> class Fn, class... Args>
    struct detector<Default, std::void_t<Fn<Args...>>, Fn, Args...> : std::true_type {
      using type = Fn<Args...>;
    };
  }  // namespace detail

  /**
   * @brief Utility to detect if a meta function has substitution failure.
   *
   * False if ``Fn<Args...>`` triggers substitution failure, true otherwise.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/utility/detected.cpp
   *    :code:
   *
   * \endrst
   *
   * @tparam Fn A template template meta-function.
   * @tparam Args Meta arguments to invoke with meta function.
   */
  template <template <class...> class Fn, class... Args>
  inline constexpr bool is_detected_v = detail::detector<void, void, Fn, Args...>::value;

  /**
   * @brief Utility to get the result of a meta function or a default value.
   *
   * If ``Fn<Args...>`` triggers substitution failure returns ``Default`` otherwise, returns ``Fn<Args...>``.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/utility/detected.cpp
   *    :code:
   *
   * \endrst
   *
   * @tparam Default The type to return in case of substitution failure.
   * @tparam Fn A template template meta-function.
   * @tparam Args Meta arguments to invoke with meta function.
   */
  template <class Default, template <class...> class Fn, class... Args>
  using detected_or_t = typename detail::detector<Default, void, Fn, Args...>::type;

  namespace detail {

    //
    template <typename...>
    struct find : std::false_type {};

    template <typename T, typename... Ts>
    struct find<T, T, Ts...> : std::true_type {};

    template <typename T, typename U, typename... Ts>
    struct find<T, U, Ts...> : find<T, Ts...> {};

  }  // namespace detail

  /**
   * @brief Test if a parameter pack contains a type.
   *
   * @tparam Target Type to search for.
   * @tparam Candidates Types to search through.
   */
  template <typename Target, typename... Candidates>
  inline constexpr bool contains_v = detail::find<Target, Candidates...>::value;

  // ------------------- Small functions ---------------- //

  /**
   * @brief Conditional printing utility.
   *
   * If ``cond`` is true forwards ``fmt`` and ``args...`` to ``fmt::print``. Useful for printing debug messages.
   */
  template <typename... Args>
  auto dprint(bool cond, fmt::format_string<Args...> fmt, Args &&...args) -> void {
    if (cond) {
      fmt::print(fmt, std::forward<Args>(args)...);
    }
  }

  /**
   * @brief Utility to reverse argument order to ``std::visit``.
   *
   * @param v A ``std::variant`` to pass to the visitor.
   * @param f A callable that accepts every possible alternative in the variant ``v``.
   *
   * @return The result of calling ``std::visit(std::forward<F>(f), std::forward<V>(v))``.
   */
  template <typename V, typename F>
  auto visit(V &&v, F &&f) -> decltype(auto) {
    return std::visit(std::forward<F>(f), std::forward<V>(v));
  }

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
   * Perform a `static_cast` from type `T` to `R` with bounds checking in debug builds.
   *
   * @tparam T Type of input ``x``, must be integral.
   * @tparam R Target type to cast to, must be integral.
   *
   * @param x The value to cast to a new type.
   *
   * @return Exactly ``static_cast<R>(x)``.
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

  namespace detail {
    template <int N = 0, typename F, typename... Args, typename T>
    void template_for_impl(Arr<T> const &beg, Arr<T> const &end, F const &f, Args... args) {
      if constexpr (N == spatial_dims) {
        if constexpr (std::is_invocable_v<F const &, Args...>) {
          std::invoke(f, args...);
        } else if constexpr (std::is_invocable_v<F const &, Arr<T>>) {
          std::invoke(f, Arr<T>{args...});
        } else {
          static_assert(always_false<F>, "template_for()'s function argument 'f' not invokable with indices or Arr<...>");
        }
      } else {
        for (T i = beg[spatial_dims - 1 - N]; i < end[spatial_dims - 1 - N]; i++) {
          template_for_impl<N + 1>(beg, end, f, i, args...);
        }
      }
    }
  }  // namespace detail

  /**
   * @brief Invoke ``f`` with every combination of indexes between ``beg`` and ``end``.
   *
   * \rst
   *
   * Effectively expands to:
   *
   * .. code::
   *
   *    for(T n = beg[spatial_dims]; n < end[spatial_dims]; n++){
   *        // ... //
   *        for(T j = beg[1]; j < end[1]; j++){
   *            for(T i = beg[0]; i < end[0]; i++){
   *                f(i, j,..., n);
   *            }
   *        }
   *    }
   *
   * Additionally, this will detect if ``f`` is callable with an array of indexes, e.g. (in the notation from above):
   *
   * .. code::
   *
   *    f(Arr<int>{i, j, ..., n});
   *
   * \endrst
   *
   * @tparam T The scalar type of the input arrays.
   *
   * @param beg Array of loop start indexes.
   * @param end Array of loop end indexes.
   * @param f Invokable to call with every combination of indexes.
   */
  template <typename T, typename F>
  void template_for(Arr<T> const &beg, Arr<T> const &end, F const &f) {
    detail::template_for_impl(beg, end, f);
  }

  /**
   * @brief Get the signed size of a container.
   *
   * See https://en.cppreference.com/w/cpp/iterator/size
   *
   * @param c Container to find the size of.
   * @return The result of ``c.size()`` cast to an appropriate signed type.
   */
  template <typename C>
  constexpr auto ssize(C const &c) -> std::common_type_t<std::ptrdiff_t, std::make_signed_t<decltype(c.size())>> {
    using R = std::common_type_t<std::ptrdiff_t, std::make_signed_t<decltype(c.size())>>;
    return static_cast<R>(c.size());
  }

  // ------------------- Math functions ---------------- //

  /**
   * @brief Test if two floating point numbers are close.
   *
   * @tparam T Type of inputs, must be a floating point type.
   *
   * @param a First input.
   * @param b Second input.
   *
   * @param atol Tolerance of absolute difference between ``a`` and  ``b`` for them to be close.
   * @param ftol Tolerance of fractional difference between ``a`` and  ``b`` for them to be close.
   *
   * @return ``true`` if ``a`` and  ``b`` are within ``atol`` or fractionally within ``ftol`` of each other.
   */
  template <typename T>
  constexpr auto near(T a, T b, T atol = 1e-10, T ftol = 0.0001) -> std::enable_if_t<std::is_floating_point_v<T>, bool> {
    return std::abs(a - b) < atol || std::abs(a - b) <= ftol * std::max(std::abs(a), std::abs(b));
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
   * @tparam T Scalar type of the input array, must be arithmetic.
   * @tparam N Number of columns in the input array.
   *
   * @param x The input array.
   *
   * @return The shifted prefix product (see above) of ``x``.
   */
  template <typename T, int N>
  auto product_scan(Eigen::Array<T, N, 1> x) -> std::enable_if_t<std::is_arithmetic_v<T>, Eigen::Array<T, N, 1>> {
    T prod = 1;
    for (Eigen::Index i = 0; i < x.size(); i++) {
      prod *= std::exchange(x[i], prod);
    }
    return x;
  }

  /**
   * @brief Compute integer powers of arithmetic types at compile time.
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
   *
   * @tparam Exp The exponent.
   * @tparam T The type of ``x``, must be arithmetic.
   *
   * @param x The input parameter.
   *
   * @return ``x``^``Exp``.
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
   * @brief Generalised dot-product between two \c Eigen objects.
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
   *
   * @param a Array-like input.
   * @param b Array-like input.
   *
   * @return The generalised dot-product (see above) of ``a`` and ``b``.
   */
  template <typename E1, typename E2>
  constexpr auto gdot(E1 const &a, E2 const &b) {
    return (a.array() * b.array()).sum();
  }

  /**
   * @brief Generic squared Frobenius norm of an \c Eigen object.
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
   *
   * @param r Array-like input.
   *
   * @return The squared Frobenius norm (see above) of ``r``.
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
   *
   * @param r Array-like input.
   *
   * @return The Frobenius norm (see above) of ``expr``.
   */
  template <typename E>
  constexpr auto gnorm(E &&r) {
    return std::sqrt(gnorm_sq(std::forward<E>(r)));
  }

  /**
   * @brief Compute the unit normal of a hyperplane through a set of points.
   *
   * \rst
   *
   * See: `StackExchange
   * <https://math.stackexchange.com/questions/2301110/fastest-way-to-find-equation-of-hyperplane>`_
   *
   * \endrst
   *
   * @tparam Scalar The scalar type of the input matrix.
   * @tparam N The dimensions of the square-matrix ``p``, must not be ``Eigen::Dynamic``.
   *
   * @param P The square input-matrix, whose columns are the N-points the hyper plane will pass through.
   *
   * @return The unit normal of a hyperplane passing through the columns of ``P``.
   */
  template <typename Scalar, int N>
  auto hyperplane_normal(Eigen::Matrix<Scalar, N, N> const &P) -> std::enable_if_t<N != Eigen::Dynamic, Eigen::Vector<Scalar, N>> {
    //
    Eigen::Matrix<Scalar, N, N + 1> H = Eigen::Matrix<Scalar, N, N + 1>::Ones();

    H.template topLeftCorner<N, N>() = P.transpose();

    Eigen::FullPivLU<decltype(H)> lu(H);

    if (auto kd = lu.dimensionOfKernel(); kd != 1) {
      throw error("Dimension is kernel={}", kd);
    }

    Eigen::Matrix<Scalar, N + 1, 1> ker = lu.kernel();

    ASSERT(near(ker[N], 0.0), "Homogeneous coordinate of the kernel = {} but should be zero!", ker[N]);

    return ker.template head<N>().normalized();
  }

  // ------------------- Classes ---------------- //

  /**
   * @brief A wrapper around a ``std::vector``.
   *
   * Does bounds checking in debug and has a signed size_type.
   */
  template <typename T>
  class Vector : private std::vector<T> {
  private:
    using Base = std::vector<T>;
    using size_type = typename Base::size_type;

  public:
    /**
     * @brief Lib cereal serialization support.
     */
    template <class Archive>
    void serialize(Archive &archive) {
      archive(static_cast<Base &>(*this));
    }

    /**
     * @brief Construct a new empty Vector.
     */
    Vector() = default;

    /**
     * @brief Construct a new Vector object containing ``size`` default initialised elements.
     *
     * @param size The size of the new container.
     */
    explicit Vector(Eigen::Index size) : Base(safe_cast<std::size_t>(size)) {
      verify(size >= 0, "{} is not a valid size for a Vector", size);
    }

    /**
     * @brief Replaces the contents with ``count`` copies of value ``value``.
     *
     * @param count The new size of the container.
     * @param value The value to initialize elements of the container with.
     */
    auto assign(Eigen::Index count, T const &value) -> void {
      verify(count >= 0, "Cannot assign {} values", count);
      Base::assign(static_cast<size_type>(count), value);
    }

    /**
     * @brief Returns a reference to the element at specified location ``pos``.
     *
     * Bounds checking in DEBUG builds.
     *
     * @param pos Position of the element to return.
     *
     * @return A mutable reference to the element.
     */
    auto operator[](Eigen::Index pos) -> T & {
      ASSERT(pos >= 0 && pos < size(), "Position, {}, is OoB in Vector length {}", pos, size());
      return Base::operator[](static_cast<std::size_t>(pos));
    }

    /**
     * @brief Returns a reference to the element at specified location ``pos``.
     *
     * Bounds checking in DEBUG builds.
     *
     * @param pos Position of the element to return.
     *
     * @return A constant reference to the element.
     */
    auto operator[](Eigen::Index pos) const -> T const & {
      ASSERT(pos >= 0 && pos < size(), "Position, {}, is OoB in Vector length {}", pos, size());
      return Base::operator[](static_cast<std::size_t>(pos));
    }

    using Base::begin;
    using Base::end;

    /**
     * @brief Fetch the number of elements in the vector.
     *
     * @return The number of elements in the Vector.
     */
    auto size() const -> Eigen::Index { return safe_cast<Eigen::Index>(Base::size()); }

    /**
     * @brief Clears the contents.
     */
    auto clear() -> void { Base::clear(); }

    /**
     * @brief Resizes the container to contain count elements.
     *
     * If the current size is greater than count, the container is reduced to its first count elements. If the current size is less
     * than count, additional default-inserted elements are appended.
     *
     * @param count New size of the container
     */
    auto resize(Eigen::Index count) -> void {
      verify(count >= 0, "Cannot assign {} values", count);
      Base::resize(static_cast<std::size_t>(count));
    }

    using Base::emplace_back;
    using Base::push_back;
  };

  /**
   * @brief Basic implementation of a Golang like defer.
   *
   *
   * \tparam F The nullary invocable's type, this **MUST** be deducted through CTAD by the deduction guide and it must be ``noexcept``
   * callable.
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
     * @param f Nullary invocable forwarded into object and invoked by destructor.
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
   * @param f Invocable to call.
   * @param args Arguments to invoke \c f with.
   * @return The result of invoking \c f with \c args... .
   */
  template <typename F, typename... Args>
  auto timeit(std::string_view name, F &&f, Args &&...args) -> std::invoke_result_t<F &&, Args &&...> {
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

      fmt::print("{} | {:>4} {:>5} {:>5} {:>5}\n", name, sec, mil, mic, nan);
    };

    if constexpr (std::is_void_v<std::invoke_result_t<F &&, Args &&...>>) {
      std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    }
  }

}  // namespace fly
