
#include <array>
#include <iostream>

#include "libfly/io/gsd.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/utility/core.hpp"

int main() {
  //
  using namespace fly;

  io::FileGSD file("build/test.gsd", io::create);

  if (true) {
    system::Box box(Mat::Identity(), Arr<bool>::Constant(true));

    // fmt::print("test {}", Mat<int>{});

    box.canon_image(Vec{10, 10, 10});

    system::SoA<Position> atom(4);

    atom(r_, 0) = Vec{0, 0, 0};
    atom(r_, 1) = Vec{1, 0, 0};
    atom(r_, 2) = Vec{0, 1, 0};
    atom(r_, 3) = Vec{0, 0, 1};

    for (int i = 0; i < 10; i++) {
      timeit("dump all", [&] { file.dump(box, atom, atom); });

      atom[r_](Eigen::lastN(3 * 3)) += 0.1;
    }

    // file.write(box, atom, r_, v_);

  } else {
    //
  }

  return 0;
}
