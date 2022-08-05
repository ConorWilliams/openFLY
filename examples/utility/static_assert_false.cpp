#include "libfly/utility/core.hpp"

template <typename T> struct referance_banned {};

template <typename T> struct referance_banned<T&> {
  // This static_assert will always fail but the compiler cannot deduce this
  // until it is instanciated as fly::always_false<T> could be specialised
  // at any point.
  static_assert(fly::always_false<T>, "T cannot be an l-value referances!");
};
