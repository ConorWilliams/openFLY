

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
  FeH.set(1, tp_, "H");
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

fly::system::SoA<TypeID, Position> explicit_V(std::vector<Vec> const &vac,
                                              system::viewSoA<TypeID, Position> cell) {
  //
  fly::system::SoA<TypeID, Position> special(cell.size() + fly::ssize(vac));

  special[id_].head(cell.size()) = cell[id_];  // Copy cells types
  special[id_].tail(fly::ssize(vac)) = 2;      // TypeID for vacacies == 2

  special[r_].head(cell[r_].size()) = cell[r_];

  Eigen::Index x = cell.size();

  for (auto const &v : vac) {
    special(r_, x++) = v;
  }

  return special;
}

struct Result {
  double v_v;  ///< The maximum of the V-V neighrest-neighbour distances. E.G for each vacancy compute the
               ///< distance to its closest neighbour then v_v is the maximum of these.
  double v_h;  ///< The minimum V-H distance.
};

/**
 * @brief
 */
Result distances(system::Box const &box, std::vector<Vec> const &vac, Vec hy) {
  //
  auto mi = box.slow_min_image_computer();

  Result r{
      0,
      std::numeric_limits<double>::max(),
  };

  for (auto const &v : vac) {
    r.v_h = std::min(r.v_h, mi(hy, v));

    for (auto const &n : vac) {
      r.v_v = std::max(r.v_v, mi(n, v));
    }
  }

  return r;
}

int main() {
  //

  // double a;

  // fmt::print("a = ");

  // std::cin >> a;

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(), {6, 6, 6});

  DetectVacancies detect(4, perfect.box(), perfect);

  system::Supercell cell = remove_atoms(perfect, {1});

  Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};

  cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_H, false)});

  //   cell = add_atoms(
  //       cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857
  //       / 4 * 2 + 3.14}, false, 0)});

  //   cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen, Hash>(1, {0.992116, 6.01736, 4.56979},
  //   false, 0)});

  //   cell(r_, 431) = Vec{4.5896, 4.5892, 5.91909};
  //   cell(r_, 0) = Vec{4.45211, 4.45172, 4.2526};

  fly::io::BinaryFile file("build/gsd/sim.gsd", fly::io::create);

  auto vac = detect.detect_vacancies(cell);

  fmt::print("Found {} vacancies @{:::.2f}\n", vac.size(), vac);

  auto const N = vac.size();

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());

    auto special = explicit_V(vac, cell);

    file.write("particles/N", fly::safe_cast<std::uint32_t>(special.size()));

    file.write(id_, special);
    file.write(r_, special);
  });

  kinetic::SKMC runner = {
      {
          .debug = true,
          .fread = "build/gsd/tmp.bin",
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
                  {
                    .debug = false, 
                    .symmetric = false,
                  }, 
                  std::ifstream{"data/wen.eam.fs"}
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

                auto v2 = detect.detect_vacancies(tmp);

                verify(v2.size() == N, "Num v changed");

                fmt::print("Found {} vacancies @{:::.2f}\n", v2.size(), v2);

                auto dist = distances(cell.box(), v2, post(r_, post.size() - 1));

                fmt::print("Max V-V = {:.3e}, min V-H = {:.3e}\n", dist.v_v, dist.v_h);

                auto vpost = explicit_V(v2, tmp);

                file.commit([&] {
                  auto copy = vpost;
                  copy[r_].head(pre[r_].size()) = pre[r_];
                  file.write(r_, copy);
                });

                file.commit([&] { file.write(r_, vpost); });

                fmt::print("Just wrote frame index No. {}\n", file.n_frames() - 1);

                if (dist.v_v > 6) {
                  ++count;
                } else {
                  count = 0;
                }

                return count >= 2;
              });

  fmt::print("It took {:.3e}s for H do detrap\n", d_time);

  return 0;
}
