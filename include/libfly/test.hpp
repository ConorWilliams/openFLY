#pragma once

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
 * Here are some details about builder: $a = a \times b$.
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