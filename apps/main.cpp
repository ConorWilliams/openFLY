
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

using namespace fly;

// Implementations of potential are concrete

struct A {
  void foo(system::SoA<Mass const&>){};
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
