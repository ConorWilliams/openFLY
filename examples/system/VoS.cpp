
#include "libfly/system/VoS.hpp"

#include "libfly/system/atom.hpp"  //< MemTag, atom

void example_VoS() {
  // Define a member to represent position that is a vector of 3 doubles.
  struct xyz : fly::system::MemTag<double, 3> {};

  fly::system::VoS<xyz> atoms;  // Empty VoS

  atoms.emplace_back({1, 2, 3});  // Add an atom at position {1, 2, 3}

  fly::system::Atom<xyz>& first_atom = atoms[0];  // Referance to first atom.

  atoms.push_back(first_atom);  // Make a second atom a copy of the first.
}