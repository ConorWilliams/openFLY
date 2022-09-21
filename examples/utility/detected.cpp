#include <utility>

#include "libfly/utility/core.hpp"

template <typename T>
// Will error if T does not have foo member.
using has_foo = decltype(std::declval<T>().foo());

// If T supports the .foo() method call and return it, otherwise throw an error.
template <typename T>
auto foo_or_throw(T any) -> fly::detected_or_t<void, has_foo, T> {
  if constexpr (fly::is_detected_v<has_foo, T>) {
    return any.foo();
  } else {
    throw fly::error("Type 'T' does not have a foo method");
  }
}