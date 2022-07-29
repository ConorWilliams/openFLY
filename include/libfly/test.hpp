#pragma once

/**
 * @brief A small kitty.
 *
 */
struct cat {
  /**
   * @brief stores the time of day
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
 * Example of math
 *
 *
 *
 * Example of an example:
 *
 * \verbatim embed:rst:leading-asterisk
 * .. include:: ../examples/main.cpp
 *    :code:
 * \endverbatim
 *
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
 * Here are some details about builder: See the code in
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