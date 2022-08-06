#include "libfly/utility/core.hpp"

#include <catch2/catch_test_macros.hpp>
#include <type_traits>

static_assert(std::is_same_v<fly::first_t<int, void, float>, int>);

static_assert(!fly::always_false<int>);

static_assert(!fly::always_false<int, float>);

TEST_CASE("ipow", "[utility]") {
  //
  CHECK(fly::ipow<0>(10) == 1);
  CHECK(fly::ipow<1>(10) == 10);
  CHECK(fly::ipow<2>(10) == 100);
  CHECK(fly::ipow<3>(10) == 1000);
  CHECK(fly::ipow<4>(10) == 10000);
  CHECK(fly::ipow<5>(10) == 100000);
}

TEST_CASE("defer", "[utility]") {
  //
  int i = 0;

  {
    fly::Defer _ = [&i]() noexcept { i = 1; };

    CHECK(i == 0);
  }

  CHECK(i == 1);

  try {
    fly::Defer _ = [&i]() noexcept { i = 2; };
    CHECK(i == 1);
    throw "!";
  } catch (...) {
    CHECK(i == 2);
  }
}