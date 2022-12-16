
// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

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

#include "libfly/utility/lattice.hpp"

#include <cstddef>

namespace fly {

  DetectVacancies::DetectVacancies(double r_max,
                                   system::Box const &box,
                                   system::viewSoA<TypeID, Position> perfect_lat)
      : m_list(box, r_max), m_perfect(perfect_lat) {
    verify(perfect_lat.size() > 0, "Perfect lattice must have some atoms!");
    m_tp = perfect_lat(id_, 0);
    verify((perfect_lat[id_] == m_tp).all(), "Perfect lattice must only contain atoms of typeID={}", m_tp);

    for (Eigen::Index i = 0; i < m_perfect.size(); i++) {
      m_perfect(r_, i) = box.canon_image(m_perfect(r_, i));
    }
  };

  std::vector<Vec> DetectVacancies::detect_vacancies(system::viewSoA<TypeID, Position> lat) {
    //
    Eigen::Index count = count_type(lat);

    verify(count <= m_perfect.size(),
           "Input lattice (n={}) has more atoms than perfect lattice (n={})",
           count,
           m_perfect.size());

    centroid_align(m_perfect, lat);

    if (auto n = m_perfect.size() + count; n != m_combo.size()) {
      m_combo.destructive_resize(n);
    }

    m_combo[r_].head(m_perfect.size() * Position::size()) = m_perfect[r_];

    {  // Copy defective into combo

      Eigen::Index x = m_perfect.size();

      for (Eigen::Index i = 0; i < lat.size(); i++) {
        if (lat(id_, i) == m_tp) {
          m_combo(r_, x) = lat(r_, i);
          x += 1;
        }
      }
    }

    m_combo[i_] = -1;

    m_list.rebuild(m_combo, 1);

    for (Eigen::Index i = m_perfect.size(); i < m_combo.size(); i++) {
      //
      double r_min = std::numeric_limits<double>::max();
      Eigen::Index i_min = -1;

      m_list.for_neighbours(i, rcut, [&](auto n, auto r, auto const &) {
        if (n < m_perfect.size() && r < r_min) {
          // If lattice site and closer
          r_min = r;
          i_min = n;
        }
      });

      verify(i_min != -1, "Could not assign atom to perfect lattice site!");

      verify(m_combo(i_, i_min) == -1, "Many-2-one assignment!");

      m_combo(i_, i_min) = i;  // Assignment
    }

    std::vector<Vec> out;

    for (Eigen::Index i = 0; i < m_perfect.size(); i++) {
      if (m_combo(i_, i) == -1) {
        out.emplace_back(m_combo(r_, i));
      }
    }

    return out;
  }

}  // namespace fly