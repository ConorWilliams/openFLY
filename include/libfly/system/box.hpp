#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-2.0

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <cstddef>
#include <nonstd/span.hpp>
#include <vector>

#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file box.hpp
 *
 * @brief Static polymorphism for box.
 */

namespace fly::system {

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
     * A span containing the indexes of every cell adjacent to the nth cell, does not include the ``n``th cell.
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

}  // namespace fly::system

//    * The canonical grid cell is the parallelepiped of space spanned by the extents of the canonical cell, canon_grid_pos() maps a
//    * point into the canonical cell and then nudges it by one grid cell in the (1,1,1) direction. This leaves room for a layer of
//    * ghost atom cells around the canonical cell.