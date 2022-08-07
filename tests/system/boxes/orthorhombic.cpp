#include "libfly/system/boxes/orthorhombic.hpp"

#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <random>

#include "libfly/system/atom.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("Orthorhombic::min_image", "[system]") {
  //
  using namespace fly;
  using namespace fly::system;

  Orthorhombic box{Arr<Position::scalar_t>::Constant(10), Arr<bool>::Constant(true)};

  Vec<Position::scalar_t> a = Vec<Position::scalar_t>::Constant(1);
  Vec<Position::scalar_t> b = Vec<Position::scalar_t>::Constant(9);

  Vec<Position::scalar_t> m = box.min_image(a, b);

  Vec<Position::scalar_t> x = Vec<Position::scalar_t>::Constant(-2);

  REQUIRE(norm(m - x) < 0.001);
}

TEST_CASE("Orthorhombic::canon_image", "[system]") {
  //
  using namespace fly;

  std::mt19937 gen(33);
  std::uniform_real_distribution<Position::scalar_t> dis(0, 1);

  auto vrand = [&] { return Arr<Position::scalar_t>::NullaryExpr([&]() { return dis(gen); }); };

  for (std::size_t i = 0; i < 100'000; i++) {
    //
    Arr<Position::scalar_t> extents = vrand() + 1;
    Arr<bool> periodic = vrand() < .5;

    fly::system::Orthorhombic box{extents, periodic};

    // Random points inside simbox
    Vec<Position::scalar_t> a = vrand() * extents;
    Vec<Position::scalar_t> b = vrand() * extents;

    // Displace by integral number of random extents in each periodic periodic direction
    Vec<Position::scalar_t> b_prime = periodic.select(b.array() + (10 * vrand()).floor() * extents, b);

    // Check in same position
    REQUIRE(std::abs(norm(box.canon_image(b_prime) - a) - fly::norm(a - b)) < 0.001);
  }
}

TEST_CASE("Orthorhombic::Grid::cell_idx", "[system]") {
  //
  using namespace fly;

  system::Orthorhombic box{Arr<Position::scalar_t>::Constant(10), Arr<bool>::Constant(true)};

  auto grid = box.make_grid(3);

  int A = grid.cell_idx(Vec<Position::scalar_t>::Constant(0));
  int B = grid.cell_idx(Vec<Position::scalar_t>::Constant(5));
  int C = grid.cell_idx(Vec<Position::scalar_t>::Constant(9.5));

  if constexpr (spatial_dims == 3) {
    REQUIRE(A == 1 + 1 * 5 + 1 * 5 * 5);
    REQUIRE(B == 2 + 2 * 5 + 2 * 5 * 5);
    REQUIRE(C == 3 + 3 * 5 + 3 * 5 * 5);
  }

  if constexpr (spatial_dims == 2) {
    REQUIRE(A == 1 + 1 * 5);
    REQUIRE(B == 2 + 2 * 5);
    REQUIRE(C == 3 + 3 * 5);
  }
}

TEST_CASE("Orthorhombic::Grid::gen_image", "[system]") {
  using namespace fly;

  system::Orthorhombic box{Arr<Position::scalar_t>::Constant(10), Arr<bool>::Constant(true)};

  system::Orthorhombic::Grid grid = box.make_grid(3);

  {
    std::optional im = grid.gen_image<Sign::plus>(Position::matrix_t::Constant(5), 0);

    REQUIRE(!im);
  }

  {
    Position::matrix_t origin = Position::matrix_t::Constant(0);

    origin[0] += 0.2;

    REQUIRE(!grid.gen_image<Sign::minus>(origin, 0));

    auto im = grid.gen_image<Sign::plus>(origin, 0);

    REQUIRE(im);

    Position::matrix_t correct = origin;

    correct[0] += 10;

    REQUIRE(norm(correct - *im) < 0.001);
  }

  {
    Position::matrix_t origin = Position::matrix_t::Constant(0);

    origin[0] = 10 - 0.2;

    REQUIRE(!grid.gen_image<Sign::plus>(origin, 0));

    auto im = grid.gen_image<Sign::minus>(origin, 0);

    REQUIRE(im);

    Position::matrix_t correct = origin;

    correct[0] -= 10;

    REQUIRE(norm(correct - *im) < 0.001);
  }
}