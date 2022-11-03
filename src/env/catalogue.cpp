// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/env/catalogue.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "libfly/env/geometry.hpp"
#include "libfly/env/heuristics.hpp"
#include "libfly/utility/core.hpp"

namespace fly::env {

  std::vector<int> Catalogue::rebuild_impl(system::SoA<Position const &, TypeID const &, Frozen const &> const &info,
                                           Eigen::Index num_types,
                                           int num_threads) {
    // Prepare memory.
    m_real.resize(std::size_t(info.size()));

    // Prep neigh list.
    m_nl->rebuild(info, num_threads);

    Geometry<Index> scratch;

#pragma omp parallel for num_threads(num_threads) firstprivate(scratch) schedule(static)
    for (int i = 0; i < info.size(); i++) {
      //
      RelEnv &env = m_real[std::size_t(i)];

      // Build geometries
      rebuild_geo_from_nl(i, env.geo, info, *m_nl, m_opt.r_env);
      // Make fingerprint.
      env.f.rebuild(env.geo);
      //
      env.hash = canon_hash(env.geo, m_opt.r_edge, safe_cast<std::size_t>(2 * num_types), &scratch);
    }

    // Make sure every hash is in the map.
    for (auto const &elem : m_real) {
      m_cat.try_emplace(elem.hash);
    }

    bool flag = false;

    // Find in catalogue, optimising for found case
#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t i = 0; i < m_real.size(); i++) {
      if (!(m_real[i].ptr = canon_find(m_real[i]))) {
#pragma omp atomic write
        flag = true;
        if (m_opt.debug) {
          fmt::print("CAT: environment around {} is new\n", i);
        }
      }
    }

    bool missing;

#pragma omp atomic read
    missing = flag;

    if (!missing) {
      if (m_opt.debug) {
        fmt::print("CAT: No new environments\n");
      }
      return {};
    }

    // Now must operate single threaded as we modify buckets
    std::vector<int> new_idx;

    struct Pair {
      Pointer ptr;
      std::size_t hash;
    };

    std::vector<Pair> new_ptrs;

    // Insert missing
    for (std::size_t i = 0; i < m_real.size(); i++) {
      if (!m_real[i].ptr) {
        //
        auto it = std::find_if(new_ptrs.begin(), new_ptrs.end(), [&](Pair const &ref) {
          return m_real[i].hash == ref.hash && canon_equiv(m_real[i], *(ref.ptr));
          //
        });

        if (it == new_ptrs.end()) {
          // New environment.
          [[maybe_unused]] auto new_hash = canon_hash(m_real[i].geo, m_opt.r_edge, safe_cast<std::size_t>(2 * num_types), &scratch);
          ASSERT(new_hash == m_real[i].hash, "Hash has changed unexpectedly, {}!={}", new_hash, m_real[i].hash);

          m_real[i].ptr = insert(m_real[i]);
          new_idx.push_back(int(i));
          new_ptrs.push_back({*m_real[i].ptr, m_real[i].hash});

          if (m_opt.debug) {
            fmt::print("CAT: Unknown at {} is new\n", i);
          }

        } else {
          // New environment equivalent to other new environment.
          m_real[i].ptr = it->ptr;

          if (m_opt.debug) {
            fmt::print("CAT: Unknown at {} is a duplicate new\n", i);
          }
        }
        //
      }
    }

    if (m_opt.debug) {
      fmt::print("CAT: found {} new environments\n", new_idx.size());
    }

    // Update frequencies
    for (auto const &elem : m_real) {
      ASSERT(elem.ptr, "Null pointer in catalogue\n", 0);
      (**elem.ptr).m_freq++;
    }

    return new_idx;
  }

  std::optional<Catalogue::Pointer> Catalogue::canon_find(Catalogue::RelEnv &env) {
    // "it" always points to valid bucket (possibly empty)
    auto it = m_cat.find(env.hash);

    ASSERT(it != m_cat.end(), "Catalogue missing key!", 0);

    // Existing key, must search bucket for explicit match
    auto match = std::find_if(it->second.begin(), it->second.end(), [&](Env const &ref) { return canon_equiv(env, ref); });

    // If found a match, return it
    if (match != it->second.end()) {
      return Pointer(it, match - it->second.begin());
    } else {
      return {};
    }
  }

  Catalogue::Pointer Catalogue::insert(Catalogue::RelEnv &env) {
    //

    ASSERT(!canon_find(env), "Environment already in catalogue.", 0);

    auto it = m_cat.find(env.hash);

    ASSERT(it != m_cat.end(), "Catalogue missing key!", 0);

    it->second.emplace_back(Env{env.geo, env.f, m_size++, m_opt.delta_max});

    return {it, xise(it->second) - 1};
  }

  bool Catalogue::canon_equiv(Catalogue::RelEnv &mut, Env const &ref) const {
    //
    double r_min = std::min(mut.f.r_min(), ref.m_finger.r_min());

    double delta = std::min(0.4 * r_min, ref.m_delta_max);

    // Test if fuzzy keys match (fast)
    if (!ref.m_finger.equiv(mut.f, delta * m_opt.overfuzz)) {
      return false;
    }

    return static_cast<bool>(mut.geo.permute_onto(ref, delta));
  }

  auto Catalogue::set_mechs(int i, std::vector<Mechanism> const &m) -> void {
    //
    Env &env = **(m_real[std::size_t(i)].ptr);

    verify(env.m_mechs.empty(), "We already have {} mechanisms, set_mech() should only be called once.", env.m_mechs.size());

    verify(std::all_of(m.begin(),
                       m.end(),
                       [s = env.size()](Mechanism const &x) {
                         //
                         return x.delta_sp.size() == s && x.delta_fwd.size() == s;
                       }),
           "Wrong number of atoms.");

    env.m_mechs = m;
  }

  double Catalogue::refine_tol(int i, double min_delta) {
    //
    std::size_t si = safe_cast<std::size_t>(i);

    double r_min = std::min(m_real[si].f.r_min(), get_ref(i).m_finger.r_min());

    double delta = std::min(0.4 * r_min, get_ref(i).m_delta_max);

    std::optional res = m_real[si].geo.permute_onto(get_ref(i).ref_geo(), delta);

    verify(bool(res), "While tightening @{} perm failed with delta={}", delta, delta);

    double new_delta_max = std::max(min_delta, res->rmsd / 1.5);

    if (m_opt.debug) {
      fmt::print("Refining delta_max @{} from {} to {}\n", i, get_ref(i).m_delta_max, new_delta_max);
    }

    (**(m_real[si].ptr)).m_delta_max = new_delta_max;

    return new_delta_max;
  }

}  // namespace fly::env