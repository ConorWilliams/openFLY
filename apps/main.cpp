
#include <fmt/core.h>

#include <array>
#include <ctime>
#include <iostream>
#include <utility>
#include <variant>

#include "libfly/io/gsd.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/utility/core.hpp"

namespace detail {

  // See https://en.cppreference.com/w/cpp/experimental/is_detected

  template <class Default, class AlwaysVoid, template <class...> class Op, class... Args>
  struct detector : std::false_type {
    using type = Default;
  };

  template <class Default, template <class...> class Op, class... Args>
  struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> : std::true_type {
    using type = Op<Args...>;
  };
}  // namespace detail

/**
 * @brief Utility to detect the presence
 *
 * @tparam Op
 * @tparam Args
 */
template <template <class...> class Op, class... Args>
inline constexpr bool is_detected_v = detail::detector<int, void, Op, Args...>::value;

using namespace fly;

// Implementations of potential are concrete

struct A {
  void foo(system::SoA<Frozen const&>){};
};

struct B {
  void foo(system::SoA<Position const&>){};
};

class Pot {
private:
  template <typename Ptr, typename... Args>
  using foo_callable_with = decltype(std::declval<Ptr>().foo(std::declval<Args>()...));

public:
  template <typename... U, typename... V>
  void foo(system::Supercell<system::TypeMap<U...>, V...> const& cell) {
    //
    visit(pots, [&](auto& x) {
      fmt::print("Detected={}\n", is_detected_v<foo_callable_with, decltype(x), decltype(cell.soa())>);

      //   x.foo(cell.soa());
      //
    });

    // if constexpr (is_detected_v<foo_impl_callable_with, T*, decltype(any.soa())>) {
    //   return parent()->foo_impl(any.soa());
    // } else {
    //   throw error("Supercell does not supply properties required for this potential");
    // }
  }

private:
  std::variant<A, B> pots;
};

int main() {
  //

  system::Supercell<system::TypeMap<>, Mass> cell({Mat::Identity(), Arr<bool>::Constant(true)}, system::TypeMap<>(2), 100);

  Pot p;

  p.foo(cell);

  /**
   * Supercell = SoA + Box + Map
   */

  fmt::print("Done\n");

  return 0;
}
