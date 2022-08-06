#include "libfly/system/boxes/orthorhombic.hpp"

#include <catch2/catch_test_macros.hpp>
#include <random>

TEST_CASE("Orthorhombic::mini_image_norm", "[system]") {
  //
  using namespace fly;
  using namespace fly::system;

  Orthorhombic box{Arr<floating>::Constant(10), Arr<bool>::Constant(true)};

  Vec<floating> a = Vec<floating>::Constant(1);
  Vec<floating> b = Vec<floating>::Constant(9);

  Vec<floating> m = box.min_image(a, b);

  Vec<floating> x = Vec<floating>::Constant(-2);

  REQUIRE(norm(m - x) < 0.001);
}

TEST_CASE("Orthorhombic::canon_image", "[system]") {
  //
  using namespace fly;

  std::mt19937 gen(33);
  std::uniform_real_distribution<floating> dis(0, 1);

  auto vrand = [&] { return Arr<floating>::NullaryExpr([&]() { return dis(gen); }); };

  for (std::size_t i = 0; i < 100'000; i++) {
    //
    Arr<floating> extents = vrand() + 1;
    Arr<bool> periodic = vrand() < .5;

    fly::system::Orthorhombic box{extents, periodic};

    // Random points inside simbox
    Vec<floating> a = vrand() * extents;
    Vec<floating> b = vrand() * extents;

    // Displace by integral number of random extents in each periodic periodic direction
    Vec<floating> b_prime = periodic.select(b.array() + (10 * vrand()).floor() * extents, b);

    // Check in same position
    REQUIRE(std::abs(norm(box.canon_image(b_prime) - a) - fly::norm(a - b)) < 0.001);
  }
}

TEST_CASE("Orthorhombic::Grid", "[system]") {
  //
  fly::system::Orthorhombic box{{10, 10, 10}, {true, true, true}};

  auto grid = box.make_grid(3);

  using V = fly::Vec<fly::floating>;

  REQUIRE(grid.cell_idx(box.canon_image(V{0, 0, 0}) + grid.cell_offset()) == 1 + 1 * 5 + 1 * 5 * 5);
  REQUIRE(grid.cell_idx(box.canon_image(V{0, 0, 0}) + grid.cell_offset()) == 1 + 1 * 5 + 1 * 5 * 5);

  REQUIRE(grid.cell_idx(box.canon_image(V{5, 5, 5}) + grid.cell_offset()) == 2 + 2 * 5 + 2 * 5 * 5);

  REQUIRE(grid.cell_idx(box.canon_image(V{9.999, 9.999, 9.999}) + grid.cell_offset()) == 3 + 3 * 5 + 3 * 5 * 5);
}