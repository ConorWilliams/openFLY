
#include "libfly/env/geometry.hpp"

#include "libfly/system/VoS.hpp"
#include "libfly/system/property.hpp"

using Set = fly::system::VoS<fly::Position, fly::Colour>;

template <typename... Args>
auto work(Args const&...) -> void {
  // Artificial work.
}

auto count_sym(Set const& ref, double delta) -> int {
  //
  Set mut = ref;  // Make a copy of the reference set

  int count = 0;

  // Explore permutations fixing the first atom.
  fly::env::for_equiv_perms(mut, ref, delta, 1, [&](fly::Mat const& O, double rmsd) {
    // If this lambda is called then "mut" has been permuted into an equivalent
    // permutation. "O" is the matrix that maps "mut" onto "ref" such that that:
    //
    //     "rmsd" == grmsd(O, mut, ref) and "rmsd" < delta.

    // We could now do something with "rmsd" and "O".

    work(O, rmsd);  // Suppresses warnings

    // We MUST NOT mutate "ref" or "mut".

    // If we "return true" then exploration of permutations terminates.
    // If we "return false" then exploration of permutations continues.

    // Here we are counting all the permutations so we will return true.
    count++;

    return true;
  });

  return count;  // The number of approximate symmetries that ref has.
}