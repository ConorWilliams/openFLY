#include <array>
#include <iostream>

#include "libfly/io/gsd.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/timeit.hpp"

int main() {
  //
  using namespace fly;

  io::FileGSD file("build/test.gsd", io::create);

  if (true) {
    system::Box box(Mat<double>::Identity(), Arr<bool>::Constant(true));

    system::SoA<Position> atom(1);

    atom(r_, 0)[0] = 0;
    atom(r_, 0)[1] = 0;
    atom(r_, 0)[2] = 0;

    file.dump(box, atom);

  } else {
  }

  //   fly::system::Box box_read;

  //   file.load(0, box_read);

  //   file.clear();

  //   std::cout << "n " << file.n_frames() << std::endl;

  //   ASSERT(box_read == box, "same?");

  //   ASSERT(box.holding<fly::system::Orthorhombic>(), "");

  //   fly::timeit("dump box 1", [&] { file.dump(box); });

  //   fly::timeit("dump box 2", [&] { file.dump(box); });

  //   std::cout << "n " << file.n_frames() << std::endl;

  //   std::array<double, 3> tmp;

  //   file.dump_span("", 2, tmp);

  //   fly::safe_cast<unsigned int>(-1);

  return 0;
}