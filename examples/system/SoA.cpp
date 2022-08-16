
#include "libfly/system/SoA.hpp"

#include "libfly/system/property.hpp"

// Define a property to represent spin that is a scalar of type bool.
struct spin : fly::system::Property<bool, 1, 1, Eigen::Matrix> {};

// Define a property to represent position that is a vector of 3 doubles.
// Note there is a built-in property for this (Position).
struct xyz : fly::system::Property<double, 3, 1, Eigen::Matrix> {};

// Define a property to represent the inertia tensor that is a 3x3 matrix of doubles.
struct I : fly::system::Property<double, 3, 3, Eigen::Matrix> {};

// Use reference property -> pass by value.
void zero_xyz(fly::system::SoA<xyz &> view_xyz) {
  // Zero the positions
  view_xyz[xyz{}] = 0;
}

void example_SoA() {
  // SoA of 10 uninitialized atoms.
  fly::system::SoA<spin, xyz, I> atoms(10);

  // Set the spin of the fist atom.
  atoms(spin{}, 0) = true;

  // Implicitly sliced at call.
  zero_xyz(atoms);
}