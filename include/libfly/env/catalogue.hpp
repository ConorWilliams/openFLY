#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.centroid>

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
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <vector>

//
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>

#include "libfly/env/geometry.hpp"
#include "libfly/env/heuristics.hpp"
#include "libfly/env/mechanisms.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
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
   * The Catalogue is responsible for building geometries from a ``Supercell`` and associating them with
   * equivalent, historically-seen geometries via a three step matching process detailed in
   * ``Catalogue::rebuild()``.
   */
  class Catalogue {
  public:
    /**
     * @brief Used to configure the catalogue.
     */
    struct Options {
      /** @brief Initial maximum difference in L2 norm between LEs for them to be considered the same. */
      double delta_max = std::numeric_limits<double>::max();
      /** @brief Smaller to decrease false positives during fingerprint equivalence but, values < 1.0
       * introduce false negatives. */
      double overfuzz = 0.5;
      /** @brief Radius of a local environment. */
      double r_env = 5.2;
      /** @brief Maximum distance for atoms to be considered connected in the canonisation neighbour graph. */
      double r_edge = 3.0;
      /** @brief Minimum value of delta_max during refinement (smaller is considered an error). */
      double min_delta_max = 1e-7;
      /** @brief If true prints debugging info. */
      bool debug = false;

      /**
       * @brief Lib cereal serialization support.
       */
      template <class Archive>
      void serialize(Archive& archive) {
        archive(delta_max, overfuzz, r_env, r_edge, debug);
      }
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
       * @brief Get a vector of the ``Mechanism``s centred on this environment.
       *
       * Pre-condition, set_mechs() must have been called on this ``Env``.
       */
      auto get_mechs() const -> std::vector<Mechanism> const& { return m_mechs; }

      /**
       * @brief Get an index (unique to this environment) in the catalogue.
       */
      auto cat_index() const noexcept -> int { return m_index; }

      /**
       * @brief Get the internal maximum norm between this environment and a geometry for them to be
       * considered equivalent.
       */
      auto delta_max() const noexcept -> double { return m_delta_max; }

      /**
       * @brief For lib cereal.
       */
      Env() = default;

      /**
       * @brief Lib cereal serialization support.
       */
      template <class Archive>
      void serialize(Archive& archive) {
        archive(
            static_cast<Geometry<>&>(*this), m_finger, m_mechs, m_freq, m_false_pos, m_index, m_delta_max);
      }

    private:
      friend class Catalogue;

      Fingerprint m_finger;              ///< The fingerprint of this environment.
      std::vector<Mechanism> m_mechs{};  ///< The mechanisms accessible, centred on this environment.
      int m_freq = 1;                    ///< The number of times this environment has been discovered.
      int m_false_pos = 0;               ///< The number of false positives encountered.
      int m_index;                       ///< The unique index in the catalogue.
      double m_delta_max;  ///< The maximum value norm for environments to be considered equivalent.

      /**
       * @brief Construct a new ``Env`` object.
       */
      Env(Geometry<Index> const& geo, Fingerprint const& f, int index, double del)
          : m_finger(f), m_index(index), m_delta_max(del) {
        for (auto const& elem : geo) {
          this->emplace_back(elem[r_], elem[col_]);
        }
      }
    };

    /**
     * @brief Describes a transformation + permutation that maps a geometry onto itself.
     *
     * Explicitly this a the transformation such that:
     */
    struct SelfSymetry {
      Mat O;                              ///< Transformation matrix.
      std::vector<Index::scalar_t> perm;  ///< Permutation.
    };

    /**
     * @brief Compute the permutation and transformation of the ``i``th geometry that maps it onto itself.

     * This is according to the ``delta_max()`` of the environment it is currently equivalent to.
     */
    auto calc_self_syms(int i) const -> std::vector<SelfSymetry>;

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
     * Any new geometries that have not been encountered before are inserted into the catalogue and the index
     * of the atom atom which the unknown geometry is centred on is returned.
     *
     * This process is deterministic, providing none of the reference geometries are modified (allowing new
     * ones to be inserted) the same match will always be returned.
     *
     * \endrst
     *
     * @param cell The supercell, must have the Position, TypeID and Frozen properties.
     * @param num_threads The number of openMP threads to use.
     * @return A ``std::vector`` of indices corresponding to the new geometries.
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
    auto get_geo(int i) const -> Geometry<Index> const& {
      //
      std::size_t si = safe_cast<std::size_t>(i);

      ASSERT(si < m_real.size(),
             "Accessing atom {} in catalogue, out of bounds as cat has {} active atoms",
             i,
             m_real.size());

      return m_real[si].geo;
    }

    /**
     * @brief Get the reference geometry stored in the catalogue that is equivalent to the geometry around
     * atom ``i``.
     */
    auto get_ref(int i) const -> Env const& {
      //
      std::size_t si = safe_cast<std::size_t>(i);

      ASSERT(si < m_real.size(),
             "Accessing atom {} in catalogue, out of bounds as cat has {} active atoms",
             i,
             m_real.size());

      return **(m_real[si].ptr);
    }

    /**
     * @brief Set the vector of ``Mechanism``s.
     *
     * Pre-condition, this may only be called once per LE.
     */
    auto set_mechs(int i, std::vector<Mechanism> const& m) -> void;

    /**
     * @brief Tighten the tolerance of the ``i``th environment.
     *
     * This is done such that the current geometry no longer matches the reference or ``delta_max =
     * min_delta``. Additionally, this function resets the frequency and false positive counters of the
     * reference environment.
     *
     * @return The new ``delta_max`` of the environment reference environment that atom ``i`` was equivalent
     * to.
     */
    auto refine_tol(int i, double min_delta = 0) -> double;

  private:
    Mat reconstruct_impl(Mechanism const& mech,
                         int i,
                         system::SoA<Position const&, TypeID const&, Frozen const&> in,
                         system::SoA<Position&> out,
                         bool in_ready_state,
                         Eigen::Index num_types,
                         int num_threads);

  public:
    /**
     * @brief Reconstruct a mechanisms ``mech`` onto the ``i`` atom of ``in``.
     */

    /**
     * @brief Reconstruct a mechanisms.
     *
     * @param out Write the reconstructed position here.
     * @param mech The mechanism to reconstruct.
     * @param i The index of the atom to reconstruct the mechanism onto.
     * @param in The initial state of the system before the reconstruction.
     * @param in_ready_state If ``true`` this function will assume the currently loaded geo/ref of the ``i``th
     * atom match the input
     * ``in`` otherwise, the geometry will be rebuilt.
     * @param num_threads Number of openMP threads to use.
     * @return Mat The matrix that transforms ``mech`` before reconstruction onto out.
     */
    template <typename Map, typename... T>
    auto reconstruct(system::SoA<Position&> out,
                     Mechanism const& mech,
                     int i,
                     system::Supercell<Map, T...> const& in,
                     bool in_ready_state,
                     int num_threads) -> Mat {
      return reconstruct_impl(mech, i, in, out, in_ready_state, in.map().num_types(), num_threads);
    }

    /**
     * @brief Lib cereal serialization support.
     */
    template <class Archive>
    void serialize(Archive& archive) {
      archive(m_opt, m_size, m_cat);
    }

    /**
     * @brief Load a catalogue from a binary archive.
     */
    Catalogue(Options const& opt, std::ifstream& fin) {
      //
      cereal::PortableBinaryInputArchive iarchive(fin);

      Catalogue loaded(opt);

      try {
        iarchive(loaded);
      } catch (std::exception const& e) {
        throw error("Catalogue is a bad binary: {}", e.what());
      }

      std::equal_to<> eq;

      verify(eq(opt.overfuzz, loaded.m_opt.overfuzz), "Loading a catalogue with contradicting .overfuzz");
      verify(eq(opt.r_env, loaded.m_opt.r_env), "Loading a catalogue with contradicting .r_env");
      verify(eq(opt.r_edge, loaded.m_opt.r_edge), "Loading a catalogue with contradicting .r_edge");

      *this = std::move(loaded);
    }

    /**
     * @brief Write this catalogue to a binary archive.
     */
    auto dump(std::ofstream& fout) const -> void {
      if (m_opt.debug) {
        fmt::print("Dump catalogue\n");
      }

      cereal::PortableBinaryOutputArchive iarchive(fout);
      iarchive(*this);
    }

    /**
     * @brief Re-order the catalogue into frequency order.
     */
    void optimize();

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

    // Members

    Options m_opt;
    int m_size = 0;
    std::optional<system::Box> m_box;
    std::optional<neigh::List> m_nl;
    std::vector<RelEnv> m_real;
    std::map<std::size_t, std::vector<Env>, std::less<>> m_cat;

    auto rebuild_impl(system::SoA<Position const&, TypeID const&, Frozen const&> const& info,
                      Eigen::Index num_types,
                      int num_threads) -> std::vector<int>;

    /**
     * @brief Simultaneously canonise the input local environment, ``env.geo``, and find its match in the
     * catalogue.
     *
     * This is a deterministic search and will find the first match that is equivalent to ``env``.
     *
     * This function is effectively ``const`` but cannot be marked as so as it returns a mutable reference to
     * the match.
     *
     * This function requires ``env.hash exists in the map``.
     *
     * @return Empty if no match found.
     */
    std::optional<Pointer> canon_find(RelEnv& env);

    /**
     * @brief Calculate the delta from the minimum fingerprint of f_mut and ref.m_finger and ref.delta_max
     */
    double calc_delta(Fingerprint const& f_mut, Env const& ref) const {
      return std::min(0.4 * std::min(f_mut.r_min(), ref.m_finger.r_min()), ref.m_delta_max);
    }

    /**
     * @brief Simultaneously canonise the input local environment and insert it into the Catalogue.
     *
     * Should only be called if ``env`` is not already in the catalogue.
     */
    Pointer insert(RelEnv& env);

    bool canon_equiv(RelEnv& mut, Env& ref) const;

    /**
     * @brief Rebuild the ith RelEnv.
     *
     * Neighbour list must be ready.
     */
    void rebuild_env(int i,
                     Eigen::Index num_types,
                     Geometry<Index>& scratch,
                     system::SoA<Position const&, TypeID const&, Frozen const&> const& info);
  };

}  // namespace fly::env