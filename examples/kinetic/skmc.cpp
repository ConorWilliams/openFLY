#include "libfly/kinetic/skmc.hpp"

using namespace fly;

void example_skmc(system::Supercell<system::TypeMap<Mass>, Position, Frozen> cell,
                  potential::Generic potential_func) {
  // Construct am SKMC object with some non-default options.
  kinetic::SKMC runner = { 
      { 
          .debug = true,
          .opt_cache = {
              .debug = true,
              .opt_basin = {
                  .debug = true,
                  .temp = 500, ///< Set the temperature of the simulation
              },
              .opt_sb = {
                  .debug = true,
              },
          },
          .opt_master = {},
      },
      cell.box(), 
      cell.map(),
      { {}, cell.box() },
      potential_func,
      { {}, {}, cell.box() },
  };

  // Run an OLKMC simulation using the runner object.
  runner.skmc(cell,
              omp_get_max_threads(),
              [&](double time,                         ///< Total time just after system at post.
                  system::SoA<Position const &> pre,   ///< State just before mech applied.
                  double E0,                           ///< Energy of the system in state pre.
                  int atom,                            ///< Index of central atom of mechanism.
                  env::Mechanism const &mech,          ///< Chosen mechanism.
                  system::SoA<Position const &> post,  ///< Final state of system after this iteration/mech.
                  double Ef                            ///< Energy of system in state post.
              ) {
                // Called every step of the simulation.
                // Can output pre/post do processing etc.
                // Return true to stop and false to continue simulation.
                return false;
              });
}