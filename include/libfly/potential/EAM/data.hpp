#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

#include <array>
#include <cstddef>
#include <fstream>

#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/spline.hpp"

/**
 * \file data.hpp
 *
 * @brief EAM parsing and storage.
 */

namespace fly::potential {

  /**
   * @brief Parses/stores tabulated EAM data and reconstructs f, phi, v smoothly via splines.
   */
  class DataEAM {
  public:
    /**
     * @brief Options for parsing the EAM file.
     */
    struct Options {
      bool debug = false;      ///< Controls debug printing.
      bool symmetric = false;  ///< Assume symmetric phi.
    };

    /**
     * @brief Parse the LAMMPS-style eam/fs  stream and construct an DataEAM object.
     *
     * For detail on file specification see: https://docs.lammps.org/pair_eam.html
     */
    DataEAM(Options const& opt, std::ifstream in);

    /**
     * @brief Fetch the cut-off radius.
     */
    auto r_cut() const noexcept -> double { return m_rcut; }

    /**
     * @brief Fetch the TypeMap.
     *
     * The map stores the mapping from TypeID's to Types that this potential has been tabulated as.
     */
    auto type_map() const noexcept -> system::TypeMap<Mass> const& { return tmap; }

    /**
     * @brief Fetch the embedding energy function.
     *
     * @param id The TypeID of the Type whose function you want to fetch.
     */
    auto f(TypeID::scalar_t id) const -> Spline const& {
      ASSERT(id < m_n, "id={} invalid with {} types in potential", id, m_n);
      return m_f[id];
    }

    /**
     * @brief Fetch the electron density function.
     *
     * Corresponding to the atom pair with TypeID's ``a`` and ``b``.
     *
     * @param a The TypeID of the Type of the first atom of the pair.
     * @param b The TypeID of the Type of the second atom of the pair.
     */
    auto phi(TypeID::scalar_t a, TypeID::scalar_t b) const -> Spline const& { return m_phi[index(a, b)]; }

    /**
     * @brief Fetch the symmetric pair-potential function.
     *
     * Corresponding to the atom pair with TypeID's ``a`` and ``b``.
     *
     * @param a The TypeID of the Type of the first atom of the pair.
     * @param b The TypeID of the Type of the second atom of the pair.
     */
    Spline const& v(TypeID::scalar_t a, TypeID::scalar_t b) const { return m_v[sym_index(a, b)]; }

  private:
    system::TypeMap<Mass> tmap{0};

    std::size_t m_n;

    double m_rcut;

    std::vector<Spline> m_f;
    std::vector<Spline> m_phi;
    std::vector<Spline> m_v;

    std::size_t index(std::size_t i, std::size_t j) const {
      ASSERT(i < m_n && j < m_n, "{}, {} is out of bounds with {} types", i, j, m_n);
      return i + m_n * j;
    }

    std::size_t sym_index(std::size_t i, std::size_t j) const {
      return index(std::max(i, j), std::min(i, j));
    }
  };

}  // namespace fly::potential
