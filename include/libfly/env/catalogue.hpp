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

#include <cstddef>
#include <map>
#include <vector>

#include "libfly/env/geometry.hpp"
#include "libfly/env/local.hpp"

namespace fly::env {

  /**
   * @brief The Catalogue is a mapping from ``Local`` environments to mechanisms.
   */
  class Catalogue {
  public:
    /**
     * @brief Used to configure the catalogue.
     */
    struct Options {
      /** @brief Maximum difference in L2 norm between LEs for them to be considered the same. */
      double delta_max = 0.5;
      /** @brief Smaller to decrease false positives but, values < 1.0 introduce false negatives. */
      double overfuzz = 0.5;
    };

    /**
     * @brief The representation of a local environment in the catalogue.
     *
     * Very similar to ``Local`` however this environment does not store index in its Geometry.
     */
    struct Env : Geometry<> {
      //
      Fingerprint fingerprint;                  ///< The fingerprint of this environment.
      std::vector<system::VoS<Delta>> mechs{};  ///< The mechanisms accessible, centred on this environment.
      int freq = 1;                             ///< The number of times this environment has been discovered.
      double delta_mod = 1;                     ///< The tolerance modifier for this environment.

      //   Env() = default;

      /**
       * @brief Construct a new Env object.
       */
      explicit Env(Fingerprint const& f) : fingerprint(f) {}
    };

    /**
     * @brief  A pointer-like type to an ``Catalogue::Env``.
     *
     * Only invalidated if ...
     */
    struct Pointer {
    public:
      /**
       * @brief Pointer deference operation.
       */
      Env* operator->() const noexcept { return m_it->second.data() + m_offset; }

    private:
      friend class Catalogue;

      Pointer() = default;

      using iterator = std::map<Local::Key, std::vector<Env>, std::less<>>::iterator;

      Pointer(iterator it, std::ptrdiff_t offset) noexcept : m_it{it}, m_offset{offset} {}

      iterator m_it;
      std::ptrdiff_t m_offset;
    };

    /**
     * @brief Construct a new Catalogue object.
     */
    explicit Catalogue(Options const& opt) : m_opt(opt) {}

    /**
     * @brief Get the number of environments in the catalogue.
     */
    int size() const noexcept { return m_size; }

    /**
     * @brief Get the number of discrete keys in the catalogue.
     */
    std::size_t num_keys() const noexcept { return m_cat.size(); }

    /**
     * @brief Simultaneously canonise the input local environment, ``env``, and find its match in the catalogue.
     *
     * This is a deterministic search and will find the first match that is equivalent to ``env``.
     *
     * Each call to ``canon_find()`` increases the frequency count if a match is found.
     *
     * @return Empty if no match found.
     */
    std::optional<Pointer> canon_find(Local& env);

    /**
     * @brief Simultaneously canonise the input local environment and insert it into the Catalogue.
     *
     * Should only be called if ``env`` is not already in the catalogue.
     */
    Pointer insert(Local& env);

  private:
    Options m_opt;
    int m_size = 0;
    std::map<Local::Key, std::vector<Env>, std::less<>> m_cat;

    void optimize() {
      for (auto&& [k, v] : m_cat) {
        std::sort(v.begin(), v.end(), [](Env const& e1, Env const& e2) {
          // Biggest first
          return e1.freq > e2.freq;
        });
      }
    }
  };

}  // namespace fly::env