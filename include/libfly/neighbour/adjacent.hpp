#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <array>
#include <cstddef>
#include <nonstd/span.hpp>
#include <vector>

#include "libfly/utility/core.hpp"

/**
 * \file adjacent.hpp
 *
 * @brief Helper types for neighbour::List
 */

namespace fly::neighbour {

  /**
   * @brief Construct a list of adjacent cell lists.
   */
  class AdjacentCells {
  public:
    /**
     * @brief Construct a new AdjacentCells object from the ``shape`` of the Grid.
     */
    explicit AdjacentCells(Arr<int> const& shape) {
      //

      verify((shape > 0).all(), "Invalid shape: {}", shape);

      //
      m_adj_cells.resize(shape.prod());

      // ND -> 1D
      auto to_1D = [cumprod = product_scan(shape)](Arr<int> const& x) { return (x * cumprod).sum(); };

      template_for<int>(Arr<int>::Zero(), shape, [&](auto... centre) {
        //
        std::size_t slot = 0;

        Arr<int> c = {centre...};

        Arr<int> beg = (c - 1).cwiseMax(0);
        Arr<int> end = (c + 2).cwiseMin(shape);

        int n = to_1D({centre...});

        template_for(beg, end, [&](auto... offset) {
          if (int m = to_1D({offset...}); n != m) {
            m_adj_cells[n][slot++] = m;
          }
        });

        m_adj_cells[n].count = slot;
      });
    }

    /**
     * @brief Get the ``n``th adjacent cell list.
     *
     * @return nonstd::span<int const> A span containing the indexes of every cell adjacent to the ``n``th cell, does not include the
     * ``n``th cell.
     */
    nonstd::span<int const> operator[](Eigen::Index n) const { return {m_adj_cells[n].data(), m_adj_cells[n].count}; }

  private:
    struct adjbours : std::array<int, ipow<spatial_dims>(3) - 1> {
      std::size_t count = 0;
    };

    Vector<adjbours> m_adj_cells;
  };

}  // namespace fly::neighbour
