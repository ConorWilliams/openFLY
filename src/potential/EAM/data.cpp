// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/potential/EAM/data.hpp"

#include <fmt/core.h>

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string_view>

#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

namespace fly::potential {

  /**
   * @brief Throwing version of std::getline
   */
  static std::istringstream getline(std::ifstream& in) {
    if (std::string line; !std::getline(in, line)) {
      throw std::runtime_error("File terminated too soon");
    } else {
      return std::istringstream{line};
    }
  }

  /**
   * @brief Read N elements in lines of length K into a vector
   */
  static std::vector<double> read_chunked(std::ifstream& file, std::size_t n, std::size_t k) {
    //
    std::vector<double> raw;

    std::istringstream line;

    for (std::size_t i = 0; i < n; ++i) {
      if (i % k == 0) {
        line = getline(file);
      }

      double tmp;
      line >> tmp;
      raw.push_back(tmp);
    }

    return raw;
  }

  DataEAM::DataEAM(std::ifstream in) {
    //
    verify(in.good(), "Could not open eam in");

    // // Skip to 4th line
    for (int i = 0; i < 3; ++i) {
      getline(in);
    }

    {  // Parse Type's
      int n;

      auto line = getline(in);

      line >> n;

      tmap = system::TypeMap<Mass>{n};

      std::string name;

      m_n = safe_cast<std::size_t>(tmap.num_types());

      for (std::uint32_t i = 0; i < m_n; i++) {
        line >> name;
        tmap.set(i, tp_, name);
        name.clear();
      }
    }

    // // Allocate space for splines
    m_f.resize(m_n);
    m_phi.resize(m_n * m_n);
    m_v.resize(m_n * m_n);

    // // Temporaries
    std::size_t numP, numR;
    double delP, delR;

    // // Parse tabulation info
    getline(in) >> numP >> delP >> numR >> delR >> m_rcut;

    for (std::size_t i = 0; i < m_n; ++i) {
      {  // Read species info
        std::size_t atomic;
        double mass;
        getline(in) >> atomic >> mass;

        tmap.set(0, m_, mass);
      }

      // Read F
      m_f[i] = Spline{read_chunked(in, numP, 5), delP};

      // Read phi
      for (std::size_t j = 0; j < m_n; ++j) {
        m_phi[index(i, j)] = Spline{read_chunked(in, numR, 5), delR};
      }
    }

    // Read v ***IMPORTANT**** tabulated as r*v and symmetric

    for (std::size_t i = 0; i < m_n; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        //
        std::vector raw = read_chunked(in, numP, 5);

        for (std::size_t k = 0; k < raw.size(); k++) {
          raw[k] /= delR * static_cast<double>(k);
        }

        verify(!raw.empty(), "no elements!");

        raw[0] = raw[1];  // Fix-up divide by zero for spline

        m_v[sym_index(i, j)] = Spline{std::move(raw), delR};
      }
    }
  }

}  // namespace fly::potential