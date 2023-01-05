

#include <fmt/core.h>
#include <omp.h>

//

#include <algorithm>
#include <array>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/kinetic/basin.hpp"
#include "libfly/kinetic/skmc.hpp"
#include "libfly/kinetic/superbasin.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle//find.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/saddle/perturb.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

using namespace fly;

template <typename... T>
system::Supercell<system::TypeMap<>, Position, Frozen, T...> bcc_iron_motif(double a = 2.855300) {
  //
  system::TypeMap<> FeH(3);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "C");
  FeH.set(2, tp_, "V");

  Mat basis{
      {a, 0, 0},
      {0, a, 0},
      {0, 0, a},
  };

  system::Supercell motif
      = system::make_supercell<Position, Frozen, T...>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

int main() {
  //

  system::Supercell cell = motif_to_lattice(bcc_iron_motif(), {6, 6, 6});

  // Mat basis = cell.box().basis();

  // basis(2, 2) = 2 * basis(1, 1);

  // cell.set_box({basis, Arr<bool>{true, true, true}});

  // Vec r_C = {5.99564 + 2.86 / 4, 5.99576, 13.6908};

  Vec r_C = {8.85219 + 2.857 / 2, 11.7075, 8.85219};

  cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_C, false)});

  //   cell = add_atoms(
  //       cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857
  //       / 4 * 2 + 3.14}, false, 0)});

  //   cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {0.992116, 6.01736, 4.56979},
  //   false, 0)});

  //   cell(r_, 431) = Vec{4.5896, 4.5892, 5.91909};
  //   cell(r_, 0) = Vec{4.45211, 4.45172, 4.2526};

  fly::io::BinaryFile file("build/gsd/sim.gsd", fly::io::create);

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());

    file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));
    file.write("log/time", -1);
    file.write(id_, cell);
    file.write(r_, cell);
  });

  kinetic::SKMC runner = {
      {
          .debug = true,
          .fread = "build/gsd/cat.FeC.bin",
          .opt_cache = {
              .barrier_tol = 0.45,
              .debug = true,
              .opt_basin = {
                  .debug = true,
                  .temp = 300,
              },
              .opt_sb = {
                  .debug = true,
              },
          },
          .opt_master = {
              .num_threads = omp_get_max_threads(),
              .max_searches = 200,
              .max_failed_searches = 75,
              .debug = false,
          }
      },
      cell.box(),
      {
          {},
          cell.box(),
      },
      potential::Generic{
          potential::EAM{
              cell.map(),
              std::make_shared<potential::DataEAM>(
                potential::DataEAM{
                  {}, 
                  std::ifstream{"data/hepburn.eam.fs"}
                }
              ),
          },
      },
      {
          {},
          {},
          cell.box(),
      },
  };

  // minimise::LBFGS minimiser({}, cell.box());

  // potential::Generic pot{
  //     potential::EAM{
  //         cell.map(),
  //         std::make_shared<potential::DataEAM>(potential::DataEAM{
  //             {
  //                 .debug = false,
  //                 .symmetric = false,
  //             },
  //             std::ifstream{"data/hepburn.eam.fs"},
  //         }),
  //     },
  // };

  // system::SoA<Position &, PotentialGradient> mirror(cell.size());

  // mirror.rebind(r_, cell);

  // fmt::print("FoundMin?={}\n",
  //            !timeit("Min", [&] { return minimiser.minimise(mirror, cell, pot, omp_get_max_threads()); }));

  // saddle::Master mast{
  //     {
  //         .num_threads = omp_get_max_threads(),
  //         .max_searches = 200,
  //         .max_failed_searches = 75,
  //         .debug = true,
  //         .fout = &file,
  //     },
  //     cell.box(),
  //     pot,
  //     minimiser,
  //     {
  //         {.debug = true},
  //         {.debug = true},
  //         cell.box(),
  //     },
  // };

  // env::Catalogue cat({});

  // cat.rebuild(cell, omp_get_max_threads());

  // mast.find_mechs(mast.package({432}, cat), cell);

  // exit(0);

  double d_time = 0;
  int count = 0;

  runner.skmc(cell,
              omp_get_max_threads(),
              [&](double time,                        ///< Total time just after system at post
                  system::SoA<Position const &> pre,  ///< State just before mech applied
                  int atom,                           ///< Index of central atom of mechanism
                  env::Mechanism const &mech,         ///< Chosen mechanism
                  system::SoA<Position const &> post  ///< Final state of system after this iteration / mech
              ) {
                //
                d_time = time;

                fly::system::SoA<TypeID const &, Position const &> tmp(cell.size());

                tmp.rebind(r_, post);
                tmp.rebind(id_, cell);

                file.commit([&] {
                  file.write("log/time", -1);
                  file.write(r_, pre);
                });

                file.commit([&] {
                  file.write("log/time", d_time);
                  file.write(r_, post);
                });

                fmt::print("Just wrote frame index No. {}\n", file.n_frames() - 1);

                return false;
              });

  fmt::print("It took {:.3e}s for H do detrap\n", d_time);

  return 0;
}
