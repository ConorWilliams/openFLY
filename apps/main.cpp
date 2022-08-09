#include <iostream>

#include "libfly/io/gsd.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/timeit.hpp"

int main() {
  //

  fly::io::FileGSD file("build/test.gsd", fly::io::read_write);

  fly::system::Box box_read(2 * fly::Mat<double>::Identity(), fly::Arr<bool>::Constant(true));

  file.load(0, box_read);

  file.clear();

  std::cout << "n " << file.n_frames() << std::endl;

  fly::system::Box box(fly::Mat<double>::Identity(), fly::Arr<bool>::Constant(true));

  ASSERT(box_read == box, "same?");

  ASSERT(box.holding<fly::system::Orthorhombic>(), "");

  fly::timeit("dump box 1", [&] { file.dump(box); });

  std::cout << "n " << file.n_frames() << std::endl;

  return 0;
}