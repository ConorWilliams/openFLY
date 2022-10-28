

#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/saddle/find.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"

using namespace fly;

system::Supercell<system::TypeMap<>, Position, Frozen, Hash> bcc_motif() {
  //
  system::TypeMap<> FeH(2);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");

  Mat basis{
      {2.855300, 0.000000, 0.000000},
      {0.000000, 2.855300, 0.000000},
      {0.000000, 0.000000, 2.855300},
  };

  system::Supercell motif = system::make_supercell<Position, Frozen, Hash>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;
  motif[hash_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

void benchmark(saddle::Master &master,
               std::vector<env::Geometry<Index>> const &geos,
               system::SoA<Position const &, Frozen const &, TypeID const &> const &cell) {
  //
  master.find_mechs(geos, cell);
}

int main() {
  /////////////////   Initialise cell   /////////////////

  system::Supercell cell = remove_atoms(motif_to_lattice(bcc_motif(), {6, 6, 6}), {});

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

  bool done = timeit("Minimise", [&] { return minimiser.minimise(mirror, cell, pot, omp_get_max_threads()); });

  fmt::print("FoundMin?={}\n", !done);

  //   /////////////////////   Get geometries   /////////////////////

  env::Catalogue cat({});

  std::vector ix = timeit("cat.rebuild()", [&] { return cat.rebuild(cell, omp_get_max_threads()); });

  fmt::print("New envs @{}\n", ix);

  for (int i = 0; i < cell.size(); ++i) {
    cell(hash_, i) = cat.get_ref(i).cat_index();
  }

  std::vector<env::Geometry<Index>> geos;

  for (auto const &elem : ix) {
    geos.push_back(cat.get_geo(elem));
  }

  /////////////////////   WRITE    /////////////////////

  io::BinaryFile fout("build/gsd/bench.gsd", io::create);

  fout.commit([&fout, &cell] {
    fout.write(cell.box());
    fout.write(cell.map());
    fout.write("particles/N", safe_cast<std::uint32_t>(cell.size()));
    fout.write(id_, cell);
    fout.write(hash_, cell);
    fout.write(r_, cell);
  });

  ///////////////////  Set up finder  ///////////////////

  saddle::Dimer dimer{{}, {}, cell.box()};

  saddle::Master mast{
      {.num_threads = omp_get_max_threads(), .debug = true, .fout = &fout},
      cell.box(),
      pot,
      minimiser,
      dimer,
  };

  //////////////////  Do benchmarking  //////////////////

  benchmark(mast, {cat.get_geo(0)}, cell);

  return 0;
}
