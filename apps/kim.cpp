

#include <fmt/core.h>

#include <cstddef>
#include <exception>
#include <nonstd/span.hpp>

//
#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/kinetic/basin.hpp"
#include "libfly/kinetic/skmc.hpp"
#include "libfly/kinetic/superbasin.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/EAM/data.hpp"
#include "libfly/potential/KIM/kim.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle//find.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/saddle/perturb.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/hessian.hpp"
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
  system::TypeMap<> FeH(2);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");

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

template <typename... T>
using cview = fly::system::SoA<T const &...>;

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> fd_hess(system::Box const &box,
                                                              potential::Generic &pot,
                                                              cview<Position, TypeID, Frozen> data) {
  system::SoA<PotentialGradient> glo(data.size());
  system::SoA<PotentialGradient> ghi(data.size());
  system::SoA<Position> xlo(data.size());
  system::SoA<Position> xhi(data.size());

  neigh::List nl(box, pot.r_cut());

  Eigen::Index const n = data.size() * Position::size();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H(n, n);

  double const dx = 1e-4;

  for (Eigen::Index i = 0; i < n; i++) {
    //
    xlo = data;
    xhi = data;

    if (!data(fzn_, i / 3)) {
      xlo(r_, i / 3)[i % 3] -= dx / 2;
      xhi(r_, i / 3)[i % 3] += dx / 2;
    }

    nl.rebuild(xlo);
    pot.gradient(glo, data, nl);

    nl.rebuild(xhi);
    pot.gradient(ghi, data, nl);

    H.col(i).array() = (ghi[g_] - glo[g_]) / dx;
  }

  return H;
}

int main(int, const char **) {
  //

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(), {2, 2, 2});

  system::Supercell cell = remove_atoms(perfect, {1, 3});

  Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};

  cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_H, false)});

  ///
  // MEAM_LAMMPS_LeeJang_2007_FeH__MO_095610951957_001
  // EAM_Dynamo_Wen_2021_FeH__MO_634187028437_000
  potential::KIM_API::Options opt{
      .model_name = "EAM_Dynamo_Wen_2021_FeH__MO_634187028437_000",
  };

  potential::KIM_API model(opt, cell.map());

  system::SoA<PotentialGradient> grad_kim(cell.size());
  system::SoA<PotentialGradient> grad_gen(cell.size());

  system::Hessian hess_kim;
  system::Hessian hess_gen;

  {
    neigh::List nl(cell.box(), model.r_cut());

    nl.rebuild(cell, omp_get_max_threads());

    for (int i = 0; i < 5; i++) {
      timeit("KIM    ", [&] { model.gradient(grad_kim, cell, nl, 1); });
    }

    fmt::print("energy={:.10f}\n", model.energy(cell, nl, 1));

    timeit("Hessian", [&] { model.hessian(hess_kim, cell, nl, 1); });
  }

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(potential::DataEAM::Options{},
                                               std::ifstream{"data/wen.eam.fs"}),
      },
  };

  {
    neigh::List nl(cell.box(), pot.r_cut());

    nl.rebuild(cell, omp_get_max_threads());

    for (int i = 0; i < 5; i++) {
      timeit("Generic", [&] { pot.gradient(grad_gen, cell, nl, 1); });
    }

    fmt::print("energy={:.10f}\n", pot.energy(cell, nl, 1));

    timeit("Hessian", [&] { pot.hessian(hess_gen, cell, nl, 1); });
  }

  // print the norm of the gradient differences
  fmt::print("Norm of the gradient difference: {}\n", gnorm(grad_gen[g_] - grad_kim[g_]));

  system::Hessian::Matrix tmp = hess_gen.get().triangularView<Eigen::Lower>();

  tmp -= hess_kim.get().triangularView<Eigen::Lower>();

  fmt::print("Norm of the hessian difference: {}\n", gnorm(tmp));

  fmt::print("r_cut = {}\n", model.r_cut());

  // std::cout << "hess1\n" << hess_gen.get() << std::endl;

  // std::cout << "hess2\n" << hess_kim.get() << std::endl;

  // std::cout << "hess1 - hess2\n" << tmp << std::endl;

  fmt::print("Working\n");

  return 0;
}