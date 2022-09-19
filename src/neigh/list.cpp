// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/neigh/list.hpp"

// #include <algorithm>
// #include <cassert>
// #include <cstddef>
// #include <limits>
#include <optional>
#include <utility>
#include <variant>

#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

// #include "libatom/asserts.hpp"
// #include "libatom/sim_cell.hpp"
// #include "libatom/utils.hpp"

namespace fly::neigh {

  auto List::update(system::SoA<Delta const&> x) -> void {
    //
    verify(x.size() == size(), "Input to update has {} atoms but list contains {}", x.size(), size());

    constexpr auto k = Position::size();

    auto kN = k * size();

    m_atoms[r_].head(kN) -= x[dr_];  // Update real atoms directly.

    // Update positions of ghosts.
    for (Eigen::Index i = size(); i < m_num_plus_ghosts; i++) {
      m_atoms(r_, i) -= x(dr_, image_to_real(i));
    }
  }

  auto List::rebuild(system::SoA<Position const&> positions, int num_threads) -> void {
    //
    verify(num_threads > 0, "{} is not a valid number of threads for rebuild()", num_threads);
    //
    if (size() != positions.size()) {
      // Allocate space if needed
      m_atoms.destructive_resize(positions.size() * (1 + MAX_GHOST_RATIO));
      m_neigh_lists.resize(positions.size());

      // Need to re-index if size changed
      for (Eigen::Index i = 0; i < positions.size(); ++i) {
        m_atoms(Index{}, i) = i;
      }
    }
    // Copy in atoms
    visit(m_box.get(), [this, &positions](auto const& box) {
      for (Eigen::Index i = 0; i < positions.size(); ++i) {
        m_atoms(r_, i) = box.canon_image(positions(r_, i));
      }
    });

    make_ghosts();

    int num_cells = visit(m_grid, [](auto const& grid) -> Arr<int> { return grid.shape(); }).prod();

    ASSERT(num_cells > 0, "{} is not a valid number of cells", num_cells);

    // Update head.
    m_head.assign(num_cells, NONE);

    // Build LCL.
    for (Eigen::Index i = 0; i < m_num_plus_ghosts; i++) {
      //
      int i_cell = visit(m_grid, [this, i](auto const& grid) { return grid.cell_idx(this->m_atoms(r_, i)); });

      m_atoms(Next{}, i) = std::exchange(m_head[i_cell], i);
    }

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (Eigen::Index i = 0; i < size(); i++) {
      build_neigh_list(i);
    }
  }

  /**
   *
   */
  void List::make_ghosts() {
    visit(m_grid, [this](auto const& grid) {
      //
      Eigen::Index next_slot = size();

      auto insert = [this](Eigen::Index slot, Vec im, Eigen::Index j) {
        ASSERT(slot < m_atoms.size(), "Not enough space for ghost number: {}", slot - size());
        m_atoms(r_, slot) = im;
        m_atoms(i_, slot) = m_atoms(i_, j);
      };

      for (int ax = 0; ax < spatial_dims; ++ax) {
        // Only make ghosts if axis is periodic
        if (m_box.periodic(ax)) {
          //
          Eigen::Index const end = next_slot;
          //
          for (Eigen::Index j = 0; j < end; ++j) {
            if (std::optional im = grid.template gen_image<Sign::plus>(m_atoms(r_, j), ax)) {
              insert(next_slot++, *im, j);
            }
          }

          for (Eigen::Index j = 0; j < end; ++j) {
            if (std::optional im = grid.template gen_image<Sign::minus>(m_atoms(r_, j), ax)) {
              insert(next_slot++, *im, j);
            }
          }
        }
      }
      m_num_plus_ghosts = next_slot;
    });
  }

  void List::build_neigh_list(Eigen::Index i) {
    //

    m_neigh_lists[i].clear();

    auto i_cell = visit(m_grid, [this, i](auto const& grid) { return grid.cell_idx(this->m_atoms(r_, i)); });

    {
      auto n = m_head[i_cell];

      while (n != NONE) {
        // In same cell must check not-self
        if (n != i && gnorm_sq(m_atoms(r_, i) - m_atoms(r_, n)) < m_r_cut * m_r_cut) {
          m_neigh_lists[i].push_back(n);
        }
        n = m_atoms(Next{}, n);
      }
    }

    for (auto&& n_cell : m_cells[i_cell]) {
      //
      auto n = m_head[n_cell];

      while (n != NONE) {
        // In adjacent cells -- don't need to check against self
        if (gnorm_sq(m_atoms(r_, i) - m_atoms(r_, n)) < m_r_cut * m_r_cut) {
          m_neigh_lists[i].push_back(n);
        }
        n = m_atoms(Next{}, n);
      }
    }
  }

}  // namespace fly::neigh