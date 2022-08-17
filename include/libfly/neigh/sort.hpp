#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <algorithm>

#include "libfly/system/SoA.hpp"
#include "libfly/system/VoS.hpp"
#include "libfly/system/box.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file sort.hpp
 *
 * @brief Utility to enhance cache locality for fly::neigh::List.
 */

namespace fly::system::detail {

  template <typename... T>
  auto SoA_to_VoS(SoA<T...> const& soa) -> VoS<remove_cref_t<T>...> {
    //
    VoS<remove_cref_t<T>...> vos;

    for (Eigen::Index i = 0; i < soa.size(); i++) {
      vos.emplace_back(soa(T{}, i)...);
    }

    return vos;
  }

  template <typename... T>
  auto VoS_to_SoA(VoS<T...> const& vos) -> SoA<remove_cref_t<T>...> {
    //
    SoA<remove_cref_t<T>...> soa(vos.size());

    for (Eigen::Index i = 0; i < vos.size(); i++) {
      ((soa(T{}, i) = vos[i][T{}]), ...);
    }

    return soa;
  }
}  // namespace fly::system::detail

namespace fly::neigh {

  /**
   * @brief Order the atoms in a SoA according to their grid index.
   *
   * This improves cache locality for neighbour reduction operations.
   *
   * @param box System box that will be used to build the grid.
   * @param r_cut Cut of radius for atomic interactions.
   * @param soa Input SoA to order.
   * @return system::SoA<remove_cref_t<T>...> A SoA ordered according to their grid index.
   */
  template <typename... T>
  auto sort(system::Box const& box, double r_cut, system::SoA<T...> const& soa) -> system::SoA<remove_cref_t<T>...> {
    //
    system::VoS vos = system::detail::SoA_to_VoS(soa);

    std::sort(vos.begin(), vos.end(), [&box, grid = box.make_grid(r_cut)](auto const& a, auto const& b) {
      auto i_a = visit(grid, [&box, &a](auto const& g) { return g.cell_idx(box.canon_image(a[r_])); });
      auto i_b = visit(grid, [&box, &b](auto const& g) { return g.cell_idx(box.canon_image(b[r_])); });
      return i_a < i_b;
    });

    return system::detail::VoS_to_SoA(vos);
  }

}  // namespace fly::neigh
