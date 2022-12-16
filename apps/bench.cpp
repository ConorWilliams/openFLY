

#include <fmt/core.h>
#include <omp.h>

#include <cstdlib>
#include <random>

#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/find.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/hessian.hpp"
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

// Perturb in-place positions around centre and write axis,
auto perturb(system::SoA<Position &, Axis &> out,
             system::SoA<Position const &, Frozen const &> in,
             Index::scalar_t centre,
             neigh::List const &nl,
             Xoshiro &prng,
             double ostddev = 0.6,
             double r_pert = 4.0) -> double {
  //
  std::normal_distribution normal(0., 1.);

  std::normal_distribution prime(ostddev, ostddev / 3.0);

  double stddev = -1;

  while (stddev < 0) {
    stddev = prime(prng);
  }

  fmt::print("Standard deviation={}\n", stddev);

  std::normal_distribution<double> gauss(0, stddev);

  out[r_] = in[r_];
  out[ax_] = 0;

  nl.for_neighbours(centre, r_pert, [&](auto n, double r, auto const &) {
    if (!in(Frozen{}, n)) {
      out(r_, n) += Vec::NullaryExpr([&] { return gauss(prng); }) * (1. - r / r_pert);
      out(ax_, n) += Vec::NullaryExpr([&] { return normal(prng); });
    }
  });

  ASSERT(!in(Frozen{}, centre), "perturbation centred on a frozen atom {}", centre);

  out(r_, centre) += Vec::NullaryExpr([&] { return gauss(prng); });
  out(ax_, centre) += Vec::NullaryExpr([&] { return normal(prng); });

  out[ax_] /= gnorm(out[ax_]);  // normalize

  return stddev;
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

int main() {
  /////////////////   Initialise cell   /////////////////

  int N = 5;

  system::Supercell cell = remove_atoms(motif_to_lattice(bcc_iron_motif(), {N, N, N}), {1});

  Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};

  cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_H, false)});

  //

  std::random_device rd{};

  Xoshiro prng(rd);

  std::uniform_real_distribution<double> d(-0.2, 0.2);

  /////////////////////   Relax   /////////////////////

  io::BinaryFile fout("build/gsd/bench.gsd", io::create);

  fout.commit([&fout, &cell] {
    fout.write(cell.box());
    fout.write(cell.map());
    fout.write("particles/N", safe_cast<std::uint32_t>(cell.size()));
    fout.write(id_, cell);
    fout.write(r_, cell);
  });

  minimise::LBFGS minimiser({.f2norm = 1e-7}, cell.box());

  potential::Generic pot{
      potential::EAM{
          cell.map(),
          std::make_shared<potential::DataEAM>(std::ifstream{"data/wen.eam.fs"}),
      },
  };

  system::SoA<Position &, PotentialGradient> mirror(cell.size());

  mirror.rebind(r_, cell);

  fmt::print("FoundMin?={}\n",
             !timeit("Min", [&] { return minimiser.minimise(mirror, cell, pot, omp_get_max_threads()); }));

  ///////////////////  Set up finder  ///////////////////

  neigh::List nl(cell.box(), pot.r_cut());

  nl.rebuild(cell, omp_get_max_threads());

  saddle::Dimer dimer{
      {.convex_max = 5, .debug = false, .fout = nullptr},
      {.iter_max_rot = 50, .theta_tol = 1 * M_PI / 360.},
      cell.box(),
  };

  //   //////////////////  Do benchmarking  //////////////////

  saddle::Master mast{
      {.num_threads = omp_get_max_threads(), .debug = true, .fout = &fout},
      cell.box(),
      pot,
      minimiser,
      dimer,
  };

  env::Catalogue cat({});

  cat.rebuild(cell, omp_get_max_threads());

  mast.find_mechs(mast.package({cell.size() - 1}, cat), cell);

  exit(0);

  // /////////////////////////

  double stddev = 0.4;

  for (size_t i = 0; i < 200; i++) {
    //
    system::SoA<Position, Axis, TypeID const &, Frozen const &> dim(cell.size());

    dim.rebind(id_, cell);
    dim.rebind(fzn_, cell);

    auto used = perturb(dim, cell, cell.size() - 1, nl, prng, stddev);

    fmt::print("Stddev = {}, used = {}\n", stddev, used);

    auto res = dimer.find_sp(dim, dim, cell, pot, {}, 0, omp_get_max_threads());

    if (res == saddle::Dimer::Exit::success) {
      stddev = std::min(1.0, used);

      fout.commit([&fout, &dim] { fout.write(r_, dim); });
    }
  }

  return 0;
}
