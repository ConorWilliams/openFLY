
#include "libfly/system/atom.hpp"

#include "libfly/system/property.hpp"

void example_atom() {
  // Define a property to represent spin that is a scalar of type bool.
  struct spin : fly::system::Property<bool> {};

  // Define a property to represent position that is a vector of 3 doubles.
  // Note there is a built-in property for this (Position).
  struct xyz : fly::system::Property<double, 3> {};

  // Define a property to represent the inertia tensor that is a 3x3 matrix of doubles.
  struct I : fly::system::Property<double, 3, 3> {};

  // Create an atom
  fly::system::Atom<spin, xyz, I> atom{
      true,
      {1, 2, 3},
      I::matrix_t::Zero(),  // I::matrix_t is Eigen::Matrix<double, 3, 3>
  };

  //   Access and set the spin property
  atom[spin{}] = false;
}