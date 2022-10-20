#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.centroid>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

#include "libfly/env/geometry.hpp"
#include "libfly/env/heuristics.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle/find.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/mechanisms.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file catalogue.hpp
 *
 * @brief Classes responsible for mapping local environments to mechanisms.
 */

namespace fly::env {

  /**
   * @brief The Catalogue is a mapping from ``Geometry`` to local environments and mechanisms.
   *
   * The Catalogue is responsible for building geometries from a ``Supercell`` and associating them with equivalent, historically-seen
   * geometries via a three step matching process detailed in ``Catalogue::rebuild()``.
   */
  class Catalogue {
  public:
    /**
     * @brief Used to configure the catalogue.
     */
    struct Options {
      /** @brief Initial maximum difference in L2 norm between LEs for them to be considered the same. */
      double delta_max = 0.5;
      /** @brief Smaller to decrease false positives during fingerprint equivalence but, values < 1.0 introduce false negatives. */
      double overfuzz = 0.5;
      /** @brief Radius of a local environment. */
      double r_env = 5.2;
      /** @brief Maximum distance for atoms to be considered connected in the canonisation neighbour graph. */
      double r_edge = 3.0;
      /** @brief If true prints debugging info. */
      bool debug = false;
    };

    /**
     * @brief The representation of a local environment in the catalogue.
     */
    struct Env : private Geometry<> {
    public:
      /**
       * @brief Fetch the reference geometry this environment represents.
       */
      auto ref_geo() const -> Geometry<> const& { return *this; }

      /**
       * @brief Get a vector of the ``system::LocalMech``s centred on this environment.
       *
       * Pre-condition, set_mechs() must have been called on this ``Env``.
       */
      auto get_mechs() const -> std::vector<system::LocalMech> const& {
        ASSERT(!m_mechs.empty(), "Environment has not been completed", 0);
        return m_mechs;
      }

      /**
       * @brief Get the index (unique to this environment) in the catalogue.
       */
      auto cat_index() const noexcept -> int { return m_index; }

    private:
      friend class Catalogue;

      Fingerprint m_finger;                      ///< The fingerprint of this environment.
      std::vector<system::LocalMech> m_mechs{};  ///< The mechanisms accessible, centred on this environment.
      int m_freq = 0;                            ///< The number of times this environment has been discovered.
      int m_index;                               ///< The unique index in the catalogue.
      double m_delta_max;                        ///< The maximum value norm for environments to be considered equivilent.

      /**
       * @brief Construct a new ``Env`` object.
       */
      explicit Env(Geometry<Index> const& geo, Fingerprint const& f, int index, double del)
          : m_finger(f), m_index(index), m_delta_max(del) {
        for (auto const& elem : geo) {
          this->emplace_back(elem[r_], elem[col_]);
        }
      }

      /**
       * @brief Set the vector of ``system::LocalMech``s.
       *
       * Pre-condition, this may only be called once.
       */
      auto set_mech(std::vector<system::LocalMech>&& m) -> void {
        verify(m_mechs.empty(), "We already have {} mechanisms, set_mech() should only be called once.", m_mechs.size());
        verify(std::all_of(m.begin(),
                           m.end(),
                           [s = this->size()](system::LocalMech const& x) {
                             //
                             return x.delta_sp.size() == s && x.delta_fwd.size() == s;
                           }),
               "Wrong number of atoms.");
        m_mechs = std::move(m);
      }
    };

    /**
     * @brief Construct a new Catalogue object.
     */
    explicit Catalogue(Options const& opt) : m_opt(opt) {}

    /**
     * @brief Build the geometries and associate them with equivalent, historically-seen geometries.
     *
     * \rst
     *
     * This occurs via a three step matching process:
     *
     * 1. The new geometries are canonised (based on their graph representation) and a discrete key generated.
     * 2. Geometries with equal discrete keys are compared for equivalence of their fingerprints.
     * 3. If they are equivalent under (2) they are compared for equivalence via ``geometry::permute_onto()``.
     *
     * Any new geometries that have not been encountered before are inserted into the catalogue and the index of the atom atom which
     * the unknown geometry is centred on is returned.
     *
     * This process is deterministic, providing none of the reference geometries are modified (allowing new ones to be inserted) the
     * same match will always be returned.
     *
     * \endrst
     */
    template <typename Map, typename... T>
    auto rebuild(system::Supercell<Map, T...> const& cell, int num_threads = 1) -> std::vector<int> {
      if (!(m_box == cell.box())) {
        if (m_opt.debug) {
          fmt::print("CAT: reconstruct neighbour list\n");
        }
        m_box = cell.box();
        m_nl = neigh::List(*m_box, m_opt.r_env);
      }
      return rebuild_impl(cell, cell.map().num_types(), num_threads);
    }

    /**
     * @brief Convert the pathways into mechanisms and store them in the catalogue.
     *
     * This requires that this function is called following a call to ``rebuild()`` (with the same cell).
     */
    auto crunch_pathways(system::SoA<Position const&, Frozen const&> cell,
                         std::vector<typename saddle::MasterFinder::PathGroup> const& pathways) -> void;

    /**
     * @brief Get the number of local environments in the catalogue.
     */
    int size() const noexcept { return m_size; }

    /**
     * @brief Get the number of discrete keys in the catalogue.
     */
    std::size_t num_keys() const noexcept { return m_cat.size(); }

    /**
     * @brief Get a view of the (indexed) geometry, in canonical order, around atom ``i``.
     */
    auto get_geo(int i) const noexcept -> Geometry<Index> const& { return m_real[std::size_t(i)].geo; }

    /**
     * @brief Get the reference geometry stored in the catalogue that is equivalent to the geometry around atom ``i``.
     */
    auto get_ref(int i) noexcept -> Env& { return **(m_real[std::size_t(i)].ptr); }

  private:
    /**
     * @brief  A non-null pointer-like type to an ``Catalogue::Env``.
     *
     * Only invalidated if ...
     */
    struct Pointer {
    public:
      /**
       * @brief Pointer deference operation.
       */
      Env* operator->() const noexcept { return m_it->second.data() + m_offset; }

      /**
       * @brief Pointer deference operation.
       */
      Env& operator*() const noexcept { return *(m_it->second.data() + m_offset); }

    private:
      friend class Catalogue;

      Pointer() = default;

      // Maps iterators are stable.
      using iterator = std::map<std::size_t, std::vector<Env>, std::less<>>::iterator;

      Pointer(iterator it, std::ptrdiff_t offset) noexcept : m_it{it}, m_offset{offset} {}

      iterator m_it;
      std::ptrdiff_t m_offset;
    };

    struct RelEnv {
      Fingerprint f;
      Geometry<Index> geo;
      std::size_t hash;
      std::optional<Pointer> ptr;
    };

    Options m_opt;
    int m_size = 0;
    std::optional<system::Box> m_box;
    std::optional<neigh::List> m_nl;
    std::vector<RelEnv> m_real;
    std::map<std::size_t, std::vector<Env>, std::less<>> m_cat;

    void optimize() {
      for (auto&& [k, v] : m_cat) {
        std::sort(v.begin(), v.end(), [](Env const& e1, Env const& e2) {
          // Biggest first
          return e1.m_freq > e2.m_freq;
        });
      }
    }

    auto rebuild_impl(system::SoA<Position const&, TypeID const&, Frozen const&> const& info, int num_types, int num_threads)
        -> std::vector<int>;

    /**
     * @brief Simultaneously canonise the input local environment, ``env.geo``, and find its match in the catalogue.
     *
     * This is a deterministic search and will find the first match that is equivalent to ``env``.
     *
     * This function is effectively ``const`` but cannot be marked as so as it returns a mutable reference to the match.
     *
     * This function requires ``env.hash exists in the map``.
     *
     * @return Empty if no match found.
     */
    std::optional<Pointer> canon_find(RelEnv& env);

    /**
     * @brief Simultaneously canonise the input local environment and insert it into the Catalogue.
     *
     * Should only be called if ``env`` is not already in the catalogue.
     */
    Pointer insert(RelEnv& env);

    bool canon_equiv(RelEnv& mut, Env const& ref) const;
  };

}  // namespace fly::env