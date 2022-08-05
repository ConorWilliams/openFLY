#include "libfly/system/boxes/orthorhombic.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Orthorhombic::mini_image_norm", "[system]") {
  //
  using namespace fly;
  using namespace fly::system;

  Orthorhombic box{Arr<floating>::Constant(10), Arr<bool>::Constant(true)};

  Arr<floating> a = Arr<floating>::Constant(1);
  Arr<floating> b = Arr<floating>::Constant(9);

  Vec<floating> m = box.min_image(a, b);

  Vec<floating> x = Vec<floating>::Constant(-2);

  REQUIRE(norm(m - x) < 0.001);
}
