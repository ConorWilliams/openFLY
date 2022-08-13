
#include <array>
#include <ctime>
#include <iostream>

#include "libfly/io/gsd.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

int main() {
  //
  using namespace fly;

  io::BinaryFile file("build/test.gsd", io::read_write);

  if (false) {
    //
    file.clear();

    system::Box box(Mat::Identity(), Arr<bool>::Constant(true));

    // fmt::print("test {}", Mat<int>{});

    // box.canon_image(Vec{10, 10, 10});

    system::SoA<Position> atom(4);

    atom(r_, 0) = Vec{0, 0, 0};
    atom(r_, 1) = Vec{1, 0, 0};
    atom(r_, 2) = Vec{0, 1, 0};
    atom(r_, 3) = Vec{0, 0, 1};

    timeit("Write", [&] {
      file.commit([&] {
        file.write(box);
        file.write(r_, atom);
      });
    });

  } else {
    system::SoA<Position> atom(4);

    timeit("Read", [&] { file.read_to(0, r_, atom); });

    file.read_to(0, r_, atom);

    verify(atom(r_, 1) == Vec{1, 0, 0}, "Oops its {}", atom(r_, 1));
  }

  return 0;
}
