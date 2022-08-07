#include "libfly/utility/core.hpp"

#include <catch2/catch_test_macros.hpp>
#include <type_traits>

static_assert(std::is_same_v<fly::first_t<int, void, float>, int>);

static_assert(!fly::always_false<int>);

static_assert(!fly::always_false<int, float>);

TEST_CASE("hyperplane_normal", "[utility]") {
  //
  using namespace fly;
  {
    Eigen::Matrix<double, 4, 4> points{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 0},
    };

    auto n = hyperplane_normal(points);

    REQUIRE(norm(n - Eigen::Vector<double, 4>{0, 0, 0, 1}) < 0.001);
  }

  {
    Eigen::Matrix<double, 3, 3> points{
        {-1, 0, 0},
        {+1, 0, 0},
        {+1, 2, 1},
    };

    auto n = hyperplane_normal(points);

    REQUIRE(norm(n - Eigen::Vector<double, 3>{1 / std::sqrt(2.), 1 / std::sqrt(2.), 0}) < 0.001);
  }
}

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