#pragma once

/**
 * @brief A small kitty.
 *
 */
struct cat {
  /**
   * @brief stores the time of day $x = y$
   *
   */
  int meow;
};

namespace fly {
/**
 * @brief an ugly cat.
 */
struct cat {
  int meow;
};

/**
 * @brief A function
 *
 * Inserting additional reStructuredText information. ds
 *
 *  A inline formula: \f$ f(x) = a + b \f$
 *
 * A display style formula: dddd
 *
 * @f[
 * \int_a^b f(x) dx = F(b) - F(a)
 * @f]
 *
 *
 * Example usage:
 * \verbatim embed:rst:leading-asterisk
 * .. include:: ../examples/main.cpp
 *    :code:
 * \endverbatim
 *
 * @code
 * char *buffer = new char[42];
 * int charsAdded = sprintf(buffer, "Tabs are normally %d spaces\n", 8);
 * @endcode
 */
inline constexpr void foo(){};

namespace impl {
/**
 * @brief BAsic secreys
 *
 */
struct impl_S {};
} // namespace impl

} // namespace fly

/**
 * @brief A big noop.
 *
 * @tparam T Must be an int.
 */
template <typename T> void noop(int) {}

/**
 * @brief Does nothing.
 *
 * @param i Ignored by test.
 */
int test(int i);

/**
 * @brief Probably does something.
 *
 * ```bash
 * cmake --build build
 * ```
 *
 * @param ignore Pass though.
 *
 * @return true Never.
 * @return false Always.
 */
bool return_true(int ignore);

namespace stuff {

/**
 * @brief Build something.
 *
 * Here are some details about builder: $a = b$. See the code in
 * /examples/main.cpp or examples/main.cpp and the documentation for return_true
 *
\code{.cpp}
class Cpp {};
\endcode

 *
 * @tparam T Never an int.
 */
template <typename T> struct builder {
public:
  /** @brief An a */
  int a;
  int b; ///< An b.

  /**
   * @brief Do "The Thing"!
   */
  void foo();

private:
  /** @brief An a */
  int c;
  int d; ///< An b.
};

} // namespace stuff