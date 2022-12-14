

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
             double r_pert = 4.0) -> void {
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

  int N = 2;

  std::cin >> N;

  system::Supercell cell = remove_atoms(motif_to_lattice(bcc_iron_motif(), {N, N, N}), {});

  cell(fzn_, 1) = true;
  cell(fzn_, 2) = true;

  //

  std::random_device rd{};

  Xoshiro prng(rd);

  std::uniform_real_distribution<double> d(-0.2, 0.2);

  /////////////////////   Relax   /////////////////////

  minimise::LBFGS minimiser({.f2norm = 1e-12}, cell.box());

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

  io::BinaryFile fout("build/gsd/bench.gsd", io::create);

  fout.commit([&fout, &cell] {
    fout.write(cell.box());
    fout.write(cell.map());
    fout.write("particles/N", safe_cast<std::uint32_t>(cell.size()));
    fout.write(id_, cell);
    fout.write(r_, cell);
  });

  ///////////////////  Set up finder  ///////////////////

  //   ///////////////////////////////

  neigh::List nl(cell.box(), pot.r_cut());

  nl.rebuild(cell, omp_get_max_threads());

  timeit("Bench gradient single threaded", [&] { pot.gradient(mirror, cell, nl, 1); });

  system::Hessian H;

  timeit("hess comp", [&] { pot.hessian(H, cell, nl, omp_get_max_threads()); });

  auto H2 = fd_hess(cell.box(), pot, cell);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dH = (H2 - H.get()).triangularView<Eigen::Lower>();

  fmt::print("hess ERROR = {}\n", gnorm(dH));

  Eigen::Index active = 0;

  for (auto const &elem : cell[fzn_]) {
    if (!elem) {
      ++active;
    }
  }

  using M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

  M r2a = M::Zero(active * Position::size(), cell.size() * Position::size());

  Eigen::Index j = 0;

  for (Eigen::Index i = 0; i < cell.size(); i++) {
    if (!cell(fzn_, i)) {
      r2a.block<3, 3>(j, i * 3) = Mat::Identity();
      j += 3;
    }
  }

  //   std::cout << "r2a\n" << r2a << std::endl;

  Position::array_t act;

  timeit("bench  mat * vec", [&] { act = r2a * cell[r_].matrix(); });

  M reduced = r2a * H.get() * r2a.transpose();

  Eigen::SelfAdjointEigenSolver<M> quick{reduced};

  fmt::print("quick eigen values:\n");

  for (auto &&v : quick.eigenvalues().head(10)) {
    fmt::print("{}\n", v);
  }

  auto [val, vec] = timeit("eigen", [&] { return H.eigen(); });

  Eigen::Index count = 0;

  while (count < val.size() && std::abs(val[count]) < 1e-5) {
    count++;
  }

  fmt::print("count={}\n", count);

  fmt::print("analytic eigen values:\n");

  for (auto &&v : val.head(13)) {
    fmt::print("{}\n", v);
  }

  //   ////

  system::Hessian::Matrix R = vec.rightCols(val.size() - count);

  for (Eigen::Index i = 0; i < R.cols(); i++) {
    R.col(i) *= 1 / std::sqrt(val[i + count]);
  }

  system::Hessian::Matrix Rt = R.transpose();

  exit(0);

  env::Catalogue cat({});

  cat.rebuild(cell, omp_get_max_threads());

  saddle::Dimer dimer{
      {.convex_max = 25, .use_history = false, .debug = true, .fout = &fout},
      {.relax_in_convex = false},
      cell.box(),
  };

  //   //////////////////  Do benchmarking  //////////////////

  for (size_t i = 0; i < 1; i++) {
    system::SoA<Position, Axis, TypeID const &, Frozen const &> dim(cell.size());

    dim.rebind(id_, cell);
    dim.rebind(fzn_, cell);

    perturb(dim, cell, 61, nl, prng);

    dimer.find_sp(dim, dim, cell, pot, {}, 0, omp_get_max_threads(), &R, &Rt);

    fout.commit([&fout, &dim] { fout.write(r_, dim); });
  }

  return 0;
}
