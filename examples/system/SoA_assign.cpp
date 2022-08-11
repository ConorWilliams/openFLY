
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"  //< Property

void SoA_assign() {
  // Define a property to represent position that is a vector of 3 doubles.
  struct xyz : fly::system::Property<double, 3> {};

  fly::system::SoA<xyz> a(10);
  fly::system::SoA<xyz> b(10);

  a[xyz{}] = 1;
  b[xyz{}] = 2;

  fly::system::SoA<xyz &> view_a = a;
  fly::system::SoA<xyz &> view_b = b;
  fly::system::SoA<xyz &> view_b2 = b;

  view_a = view_b;                 // This does NOT set a[xyz{}] to  b[xyz{}]
  view_a = std::move(view_b);      // Neither does this
  view_a[xyz{}] = view_b2[xyz{}];  // This does :)
}