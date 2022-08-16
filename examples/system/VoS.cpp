
#include "libfly/system/VoS.hpp"

#include "libfly/system/atom.hpp"
#include "libfly/system/property.hpp"

void example_VoS() {
  // Define a property to represent position that is a vector of 3 doubles.
  // Note there is a built-in property for this (Position).
  struct xyz : fly::system::Property<double, 3, 1, Eigen::Matrix> {};

  fly::system::VoS<xyz> atoms;  // Empty VoS

  atoms.emplace_back({1, 2, 3});  // Add an atom at position {1, 2, 3}

  fly::system::Atom<xyz>& first_atom = atoms[0];  // Reference to first atom.

  atoms.push_back(first_atom);  // Make a second atom a copy of the first.
}