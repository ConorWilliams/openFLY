

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
system::Supercell<system::TypeMap<>, Position, Frozen, T...> bcc_iron_motif() {
  //
  system::TypeMap<> FeH(3);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");
  FeH.set(2, tp_, "V");

  Mat basis{
      {2.855300, 0.000000, 0.000000},
      {0.000000, 2.855300, 0.000000},
      {0.000000, 0.000000, 2.855300},
  };

  system::Supercell motif
      = system::make_supercell<Position, Frozen, T...>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

template <typename... Ts>
void ignore(Ts const &...) {}

int main() {
  //

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(), {6, 6, 6});

  DetectVacancies detect(perfect.box(), perfect);

  system::Supercell cell = remove_atoms(perfect, {1, 3});

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

  fmt::print("Found {} vacancies @{}\n", vac.size(), vac);

  fly::system::SoA<TypeID, Position> special(cell.size() + fly::ssize(vac));

  special[id_].head(cell.size()) = cell[id_];
  special[id_].tail(fly::ssize(vac)) = 2;

  auto const N = vac.size();

  {
    special[r_].head(cell.size() * Position::size()) = cell[r_];

    Eigen::Index x = cell.size();

    for (auto const &v : vac) {
      special(r_, x++) = v;
    }
  }

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());
    file.write("particles/N", fly::safe_cast<std::uint32_t>(special.size()));

    file.write(id_, special);
    file.write(r_, special);
  });

  kinetic::SKMC runner = {
      {
          .debug = true,
          .fread = "build/gsd/cat.bin",
          .opt_cache = {
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
              .max_searches = 100,
              .max_failed_searches = 50,
              .debug = true,
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
              std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
          },
      },
      {
          {},
          {},
          cell.box(),
      },
  };

  runner.skmc(cell,
              omp_get_max_threads(),
              [&](double time,
                  system::SoA<Position const &> pre,
                  int atom,
                  env::Mechanism const &mech,
                  system::SoA<Position const &> post) {
                //
                ignore(time, pre, atom, mech, pre, post);

                fly::system::SoA<TypeID const &, Position const &> tmp(cell.size());

                tmp.rebind(r_, post);
                tmp.rebind(id_, cell);

                auto v2 = detect.detect_vacancies(tmp);

                verify(v2.size() == N, "Num v changed");

                {
                  special[r_].head(post.size() * Position::size()) = post[r_];

                  Eigen::Index x = post.size();

                  for (auto const &v : v2) {
                    special(r_, x++) = v;
                  }
                }

                fmt::print("Found {} vacancies @{}\n", vac.size(), vac);

                // file.commit([&] { file.write(r_, pre); });
                // file.commit([&] { file.write(r_, post); });
                file.commit([&] { file.write(r_, special); });
              });

  return 0;
}
