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
      ASSERT((shape > 0).all(), "Invalid shape!");
      //
      m_adj_cells.resize(static_cast<std::size_t>(shape.prod()));

      // Cumulative product of m_shape.
      Arr<int> cumprod = Arr<int>::Ones();

      for (int i = 1; i < spatial_dims; i++) {
        cumprod[i] = cumprod[i - 1] * shape[i - 1];
      }

      // ND -> 1D
      auto to_1D = [=](Arr<int> const& x) { return (x * cumprod).sum(); };

      template_for(Arr<int>::Zero(), shape, [&](auto... centre) {
        //
        std::size_t slot = 0;

        Arr<int> c = {centre...};

        Arr<int> beg = (c - 1).cwiseMax(0);
        Arr<int> end = (c + 2).cwiseMin(shape);

        int n = to_1D({centre...});

        template_for(beg, end, [&](auto... offset) {
          if (int m = to_1D({offset...}); n != m) {
            m_adj_cells[(std::size_t)n][slot++] = m;
          }
        });

        m_adj_cells[(std::size_t)n].count = slot;
      });
    }

    /**
     * @brief Get the ``n``th adjacent cell list.
     *
     * @return nonstd::span<int const> A span containing the indexes of every cell adjacent to the ``n``th cell, does not include the
     * ``n``th cell.
     */
    nonstd::span<int const> operator[](int n) const {
      auto sn = (std::size_t)n;
      ASSERT(n >= 0 && sn < m_adj_cells.size(), "Invalid cell index!");
      return {m_adj_cells[sn].data(), m_adj_cells[sn].count};
    }

  private:
    struct adjbours : std::array<int, ipow<spatial_dims>(3) - 1> {
      std::size_t count = 0;
    };

    std::vector<adjbours> m_adj_cells;

    /**
     * @brief Invoke ``f`` with every tuple of indexes between ``beg`` and ``end``.
     */
    template <int N = 0, typename F, typename... Args>
    void template_for(Arr<int> const& beg, Arr<int> const& end, F const& f, Args... args) {
      if constexpr (N == spatial_dims) {
        f(args...);
      } else {
        for (int i = beg[spatial_dims - 1 - N]; i < end[spatial_dims - 1 - N]; i++) {
          template_for<N + 1>(beg, end, f, i, args...);
        }
      }
    }
  };

}  // namespace fly::neighbour
