
#include "libfly/system/atom.hpp"

#include <Eigen/src/Core/Matrix.h>

void example_atom() {
  // Define a member to represent spin that is a scalar of type bool.
  struct spin : fly::system::MemTag<bool> {};

  // Define a member to represent position that is a vector of 3 doubles.
  struct xyz : fly::system::MemTag<double, 3> {};

  // Define a member to represent the inertia tensor that is a 3x3 matrix of doubles.
  struct I : fly::system::MemTag<double, 3, 3> {};

  // Create an atom
  fly::system::Atom<spin, xyz, I> atom{
      true,
      {1, 2, 3},
      I::matrix_t::Zero(),  // I::matrix_t is Eigen::Matrix<double, 3, 3>
  };

  //   Access and set the spin property
  atom[spin{}] = false;
}