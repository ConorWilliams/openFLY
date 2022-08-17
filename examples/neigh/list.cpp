
#include "libfly/neigh/list.hpp"

#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

#define MU [[maybe_unused]]  // Silence unused variable warning.

namespace fs = fly::system;

void example_list(fs::Supercell<fs::TypeMap<>, fly::Position> const& cell) {
  //
  constexpr double r_cut = 0.5;

  fly::neigh::List nl(cell.box(), r_cut);

  int sum = 0;

  for (int i = 0; i < nl.size(); i++) {
    nl.for_neighbours(i, r_cut, [&](MU auto n, MU auto r, MU auto const& dr) {
      // Here "n"  is the index of the neighbour,
      //      "r"  is the minimum-image (MI) distance to "n" and
      //      "dr" is the position vector connecting "i" to the MI of "n".

      sum += 1;
    });
  }

  // Average number of neighbours is sum / nl.size() //
}