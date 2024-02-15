

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <omp.h>

//

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
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
#include "libfly/utility/data.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

using namespace fly;

double sym_2_mass(std::string_view symbol) {
  using namespace fly::data;

  auto it = std::find_if(atoms.begin(), atoms.end(), [&](auto const &sym) { return sym.symbol == symbol; });

  if (it == atoms.end()) {
    throw error("Symbol \"{}\" not found in the periodic table!", symbol);
  }

  return it->mass;
}

template <typename... T>
system::Supercell<system::TypeMap<Mass>, Position, Frozen, T...> bcc_iron_motif(double a) {
  //
  system::TypeMap<Mass> FeH(3);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");
  FeH.set(2, tp_, "V");

  for (std::uint32_t i = 0; i < FeH.num_types(); i++) {
    FeH.set(i, m_, sym_2_mass(FeH.get(i, tp_)));
  }

  Mat basis{
      {a, 0, 0},
      {0, a, 0},
      {0, 0, a},
  };

  system::Supercell motif = system::make_supercell<Position, Frozen, T...>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

fly::system::SoA<TypeID, Position> explicit_V(std::vector<Vec> const &vac, system::viewSoA<TypeID, Position> cell) {
  //
  fly::system::SoA<TypeID, Position> special(cell.size() + fly::ssize(vac));

  special[id_].head(cell.size()) = cell[id_];  // Copy cells types
  special[id_].tail(fly::ssize(vac)) = 2;      // TypeID for vacancies == 2

  special[r_].head(cell[r_].size()) = cell[r_];

  Eigen::Index x = cell.size();

  for (auto const &v : vac) {
    special(r_, x++) = v;
  }

  return special;
}

enum class stop_criteria {
  h_escaped,
  v_dissociated,
  both,
};

struct result {
  stop_criteria crit;
  double msd;
  double msd_c_align;
  double time;
};

auto run_sim(std::string ofname, std::string ifname, double temp, int n_vac, bool hy) -> result {
  //
  constexpr auto a = 2.855;  // eam 55 or meam 67

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(a), {8, 8, 8});

  DetectVacancies detect(0.75, perfect.box(), perfect);

  std::vector<Eigen::Index> Vs;

  switch (n_vac) {
    default:
      throw error("n_vac={} oob", n_vac);
    case 5:
      Vs.push_back(129);
    case 4:
      Vs.push_back(130);
    case 3:
      Vs.push_back(146);
    case 2:
      Vs.push_back(3);
    case 1:
      Vs.push_back(1);
  }

  system::Supercell cell = remove_atoms(perfect, Vs);

  if (hy) {
    Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};
    cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>{1, r_H, false}});
  }

  fly::io::BinaryFile file(ofname, fly::io::create);

  // Write first frame to GSD

  file.commit([&] {
    file.write(cell.box());
    file.write(cell.map());

    file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));

    file.write(id_, cell);
    file.write(r_, cell);

    file.write("log/time", -1.);

    file.write("log/energy", -1.);
    file.write("log/barrier", -1.);
    file.write("log/kinetic", -1.);

    file.write("log/vv_max", -1.);
    file.write("log/vv_min", -1.);

    file.write("log/temp", temp);  // Meta-data
  });

  kinetic::SKMC runner = {
      {
          .debug = true,
          .fread = ifname,
          .opt_cache = {
              .barrier_tol = 0.45,
              .debug = true,
              .opt_basin = {
                  .debug = true,
                  .temp = temp,
              },
              .opt_sb = {
                  .debug = true,
              },
          },
          .opt_master = {
              .num_threads = omp_get_max_threads(),
              .max_searches = 200,
              .max_failed_searches = 75,
              .debug = true,
          }
      },
      cell.box(),
      cell.map(),
      {
          {},
          cell.box(),
      },
      potential::Generic{
            potential::EAM{
                system::TypeMap<>{cell.map()},
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

  auto const min_image = cell.box().slow_min_image_computer();

  double run_time = 0;

  bool v_dis = false;
  bool h_escaped = false;

  static constexpr double v_dis_crit = 5;
  static constexpr double h_esc_crit = 5;

  double msd = 0;
  double msd_c_align = 0;

  std::size_t n = 0;

  system::SoA<Position> init(cell.size());

  runner.skmc(cell,
              omp_get_max_threads(),
              [&](double time,                         ///< Total time just after system at post
                  system::SoA<Position const &> pre,   ///< State just before mech applied
                  double E0,                           ///< Energy of the system in state pre.
                  int,                                 ///< Index of central atom of mechanism
                  env::Mechanism const &mech,          ///< Chosen mechanism
                  system::SoA<Position const &> post,  ///< Final state of system after this iteration / mech
                  double Ef                            ///< Energy of system in state post.
              ) {
                //
                if (n == 0) {
                  init = pre;
                }

                //
                run_time = time;

                fly::system::SoA<TypeID const &, Position const &> tmp(cell.size());

                tmp.rebind(r_, post);
                tmp.rebind(id_, cell);

                std::vector vac = detect.detect_vacancies(tmp);

                // ----------- V-dis testing

                double const VV = kruskal_max(vac, min_image);

                v_dis = VV > v_dis_crit;

                double vh_min = std::numeric_limits<double>::max();

                //----------- H-escape testing

                if (auto last = post.size() - 1; hy) {
                  //
                  for (auto const &v : vac) {
                    vh_min = std::min(vh_min, min_image(post(r_, last), v));
                  }

                  h_escaped = vh_min > h_esc_crit;  // H-escape criterion set here.
                }

                bool end = h_escaped || v_dis;

                if (end) {
                  //
                  msd = 0;

                  for (std::size_t i = 0; i < pre.size(); ++i) {
                    if (cell(id_, i) == 0) {
                      msd += (pre(r_, i) - init(r_, i)).squaredNorm();
                    }
                  }

                  msd_c_align = 0;

                  centroid_align(init, pre);

                  for (std::size_t i = 0; i < pre.size(); ++i) {
                    if (cell(id_, i) == 0) {
                      msd_c_align += (pre(r_, i) - init(r_, i)).squaredNorm();
                    }
                  }
                }

                // Write to GSD

                if (n < 10 || end) {
                  timeit("IO", [&] {
                    //
                    file.commit([&] {
                      file.write("particles/N", fly::safe_cast<std::uint32_t>(cell.size()));
                      file.write(r_, pre);
                      file.write("log/energy", E0);
                    });

                    file.commit([&] {
                      //
                      auto vpost = explicit_V(vac, tmp);

                      file.write("particles/N", fly::safe_cast<std::uint32_t>(vpost.size()));

                      file.write(id_, vpost);
                      file.write(r_, vpost);

                      file.write("log/time", time);
                      file.write("log/energy", Ef);
                      file.write("log/barrier", mech.barrier);
                      file.write("log/kinetic", mech.kinetic_pre);
                      file.write("log/vv_max", VV);
                      file.write("log/vv_min", vh_min);
                    });
                  });
                }

                ++n;

                return end;
              });

  fmt::print("It took {:.3e}s to terminate\n", run_time);
  fmt::print("H escaped = {}, cluster dissociated = {}\n", h_escaped, v_dis);

  if (h_escaped && v_dis) {
    return {stop_criteria::both, msd, msd_c_align, run_time};
  }

  if (h_escaped) {
    return {stop_criteria::h_escaped, msd, msd_c_align, run_time};
  }

  if (v_dis) {
    return {stop_criteria::v_dissociated, msd, msd_c_align, run_time};
  }

  throw error("Unexpected termination condition");
}

void loop_main(int n_vac, bool hy) {
  //
  fmt::print("Requested {} vacancies\n", n_vac);

  verify(0 < n_vac && n_vac <= 5, "n_vac={} is oob", n_vac);

  std::string prefix = fmt::format("cluster/V{}", n_vac);

  fmt::print("PREFIX={}\n", prefix);

  fmt::print("threads={}\n", omp_get_max_threads());

  auto out_h_s = fmt::format("{}/escape{}.csv", prefix, hy ? "" : ".no_h");
  auto out_v_s = fmt::format("{}/dis{}.csv", prefix, hy ? "" : ".no_h");

  auto out_h = fmt::output_file(out_h_s, fmt::file::CREATE | fmt::file::APPEND | fmt::file::WRONLY);
  auto out_v = fmt::output_file(out_v_s, fmt::file::CREATE | fmt::file::APPEND | fmt::file::WRONLY);

  std::size_t n = 20;         /// Number of temperatures to sample.
  std::size_t rep_h = 10;     /// Number of repetitions at each temperature.
  std::size_t rep_v = 10;     //
  std::size_t max_reps = 50;  /// Max reps at each temperature.

  double lo = 300;
  double hi = 800;

  // Want to uniformly sample inverse temperature range between 300K and 800K

  double inv_lo = 1. / lo;
  double inv_hi = 1. / hi;

  double dif = (inv_lo - inv_hi) / static_cast<double>(n - 1);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<std::size_t> dis;

  for (double inv_T = inv_lo; inv_T > inv_hi - 0.5 * dif; inv_T -= dif) {
    //
    double temp = 1 / inv_T;

    fmt::print("T={}, inv={}\n", temp, inv_T);

    std::size_t n_h = 0;
    std::size_t n_v = 0;

    if (n_vac == 1) {
      n_v = rep_v;
    }

    if (!hy) {
      n_h = rep_h;
    }

    while (n_h < rep_h || n_v < rep_v) {
      //
      if (n_h >= max_reps || n_v >= max_reps) {
        fmt::print("Reached {} escapes or dissociations at {:.1f}K\n", max_reps, temp);
        break;
      }

      //
      std::size_t uid = dis(gen);

      auto gsd_file = fmt::format("{}/gsd/{:.1f}K.u{}.gsd", prefix, temp, uid);
      auto cat_file = fmt::format("{}/cat{}.bin", prefix, hy ? ".h" : "");

      auto [term, msd, msd_c, dt] = run_sim(gsd_file, cat_file, temp, n_vac, hy);

      if (term == stop_criteria::both || term == stop_criteria::h_escaped) {
        if (!hy) {
          throw error("no h cannot escape!");
        }
        out_h.print("{:%Y-%m-%d %H:%M:%S} {} {} {} {:.5f} {:e}\n", fmt::localtime(std::time(nullptr)), uid, msd, msd_c, temp, dt);
        out_h.flush();
        ++n_h;
      }
      if (term == stop_criteria::both || term == stop_criteria::v_dissociated) {
        out_v.print("{:%Y-%m-%d %H:%M:%S} {} {} {} {:.5f} {:e}\n", fmt::localtime(std::time(nullptr)), uid, msd, msd_c, temp, dt);
        out_v.flush();
        ++n_v;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  try {
    if (argc < 2) {
      throw error("Require hydrogen (y/n)");
    }

    bool hy = false;

    if (std::string hy_arg = argv[1]; hy_arg == "y") {
      hy = true;
    } else if (hy_arg != "n") {
      throw error("Invalid hydrogen argument");
    }

    for (int i = 2; i < argc; ++i) {
      loop_main(std::stoi(argv[i]), hy);
    }

  } catch (std::exception const &e) {
    fmt::print("Caught exception: {}\n", e.what());
    return 1;
  } catch (...) {
    fmt::print("Terminated with unknown exception\n");
    return 1;
  }

  return 0;
}
