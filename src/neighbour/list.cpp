// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/neighbour/list.hpp"

// #include <algorithm>
// #include <cassert>
// #include <cstddef>
// #include <limits>
#include <optional>
#include <utility>
#include <variant>

#include "libfly/utility/core.hpp"

// #include "libatom/asserts.hpp"
// #include "libatom/sim_cell.hpp"
// #include "libatom/utils.hpp"

namespace fly::neighbour {

  auto List::rebuild(system::SoA<Position const&> positions, int num_threads) -> void {
    //
    verify(num_threads > 0, "{} is not a valid number of threads for rebuild()", num_threads);

    init_and_build_lcl(positions, num_threads);

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (Eigen::Index i = 0; i < m_neigh_lists.size(); i++) {
      build_neigh_list(i);
    }
  }

  void List::init_and_build_lcl(system::SoA<Position const&> atoms, int) {
    //
    if (m_neigh_lists.size() != atoms.size()) {
      // Allocate space if needed
      m_atoms.destructive_resize(atoms.size() * (1 + MAX_GHOST_RATIO));
      m_neigh_lists.resize(atoms.size());

      // Need to re-index if size changed
      for (Eigen::Index i = 0; i < atoms.size(); ++i) {
        m_atoms(Index{}, i) = i;
      }
    }
    // Copy in atoms
    std::visit(
        [this, &atoms](auto const& box) {
          for (Eigen::Index i = 0; i < atoms.size(); ++i) {
            m_atoms(Position{}, i) = box.canon_image(atoms(r_, i));
          }
        },
        m_box.get());

    make_ghosts();

    int num_cells = std::visit([](auto const& grid) -> Arr<int> { return grid.shape(); }, m_grid).prod();

    ASSERT(num_cells > 0, "{} is not a valid number of cells", num_cells);

    // Update head.
    m_head.assign(num_cells, NONE);

    // Build LCL.
    for (Eigen::Index i = 0; i < m_num_plus_ghosts; i++) {
      //
      int i_cell = std::visit([this, i](auto const& grid) { return grid.cell_idx(this->m_atoms(r_, i)); }, m_grid);

      m_atoms(Next{}, i) = std::exchange(m_head[i_cell], i);
    }
  }

  //   void List::update_positions(Position::matrix_type const& deltas) {
  //     // Copy in atoms
  //     for (std::size_t i = 0; i < m_neigh_lists.size(); ++i) {
  //       m_atoms(Position{}, i) -= deltas.col(i);
  //     }

  //     // Update ghosts
  //     for (std::size_t i = m_neigh_lists.size(); i < m_num_plus_ghosts; i++) {
  //       m_atoms(Position{}, i) -= deltas.col(image_to_real(i));
  //     }
  //   }

  void List::build_neigh_list(Eigen::Index i) {
    //

    m_neigh_lists[i].clear();

    auto i_cell = std::visit([this, i](auto const& grid) { return grid.cell_idx(this->m_atoms(r_, i)); }, m_grid);

    double r_cut = std::visit([](auto const& grid) { return grid.r_cut(); }, m_grid);

    {
      auto n = m_head[i_cell];

      while (n != NONE) {
        // In same cell must check not-self
        if (n != i && gnorm_sq(m_atoms(r_, i) - m_atoms(r_, n)) < r_cut * r_cut) {
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
        if (gnorm_sq(m_atoms(r_, i) - m_atoms(r_, n)) < r_cut * r_cut) {
          m_neigh_lists[i].push_back(n);
        }
        n = m_atoms(Next{}, n);
      }
    }
  }

  void List::make_ghosts() {
    std::visit(
        [this](auto const& grid) {
          //
          Eigen::Index next_slot = m_neigh_lists.size();

          for (int ax = 0; ax < spatial_dims; ++ax) {
            // Only make ghosts if axis is periodic
            if (m_box.periodic(ax)) {
              //
              Eigen::Index const end = next_slot;

              for (Eigen::Index j = 0; j < end; ++j) {
                //
                std::array const images = {
                    grid.template gen_image<Sign::plus>(m_atoms(r_, j), ax),
                    grid.template gen_image<Sign::minus>(m_atoms(r_, j), ax),
                };

                for (auto&& elem : images) {
                  if (elem) {
                    Eigen::Index slot = next_slot++;
                    ASSERT(slot < m_atoms.size(), "Not enough space for ghost number: {}", slot - m_atoms.size());
                    m_atoms(r_, slot) = *elem;
                    m_atoms(i_, slot) = m_atoms(i_, j);
                  }
                }
              }
            }
          }
          m_num_plus_ghosts = next_slot;
        },
        m_grid);
  }

}  // namespace fly::neighbour