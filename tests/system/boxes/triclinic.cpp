

#include "libfly/system/boxes/triclinic.hpp"

#include <Eigen/src/Core/Diagonal.h>
#include <Eigen/src/Core/DiagonalMatrix.h>

#include <catch2/catch_test_macros.hpp>
#include <ios>
#include <iostream>
#include <random>

#include "libfly/system/atom.hpp"
#include "libfly/utility/core.hpp"

TEST_CASE("Triclinic::canon_image", "[system]") {
  //
  using namespace fly;

  std::mt19937 gen(33);
  std::uniform_real_distribution<Position::scalar_t> dis(0, 1);

  // Random numbers on interval 0 to 1
  auto vrand = [&] { return Position::matrix_t::NullaryExpr([&]() { return dis(gen); }); };

  using M = Mat<Position::scalar_t>;

  for (std::size_t i = 0; i < 1'000; i++) {
    //
    M basis = [&]() {
      M tmp = M::NullaryExpr([&]() { return dis(gen); });
      tmp.array() += 1;
      tmp = tmp.triangularView<Eigen::Upper>();
      return tmp;
    }();

    Arr<bool> periodic = vrand().array() > .5;

    fly::system::Triclinic box{basis, periodic};

    // // Random points inside simbox
    Position::matrix_t a = basis * vrand();
    Position::matrix_t b = basis * vrand();

    Position::matrix_t integer_randoms = (10 * vrand()).array().floor();

    Eigen::DiagonalMatrix<Position::scalar_t, spatial_dims, spatial_dims> diag(integer_randoms);

    M offsets = basis * diag;

    Position::matrix_t b_prime = b;

    for (int j = 0; j < 3; j++) {
      if (periodic[j]) {
        b_prime += offsets.col(j);
      }
    }

    Position::matrix_t b_undo = box.canon_image(b_prime);

    // // Check in same position
    REQUIRE(std::abs(norm(b_undo - a) - fly::norm(a - b)) < 0.001);
  }
}