#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <variant>
#include <vector>

#include "libfly/neighbour/adjacent.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"

/**
 * \file list.hpp
 *
 * @brief ...
 */

namespace fly::neighbour {

  /**
   * @brief The maximum number of ghosts neighbour cell supports is MAX_GHOST_RATIO * num_atoms
   */
  inline constexpr Eigen::Index MAX_GHOST_RATIO = 6;

  /**
   * @brief A class to contain, build and manage neighbour lists in shared memory.
   *
   * Designed with the intention of being reused List separates the building, updating and
   * using of neighbour lists (NLs). NLs can be constructed from a SimCell and then used to find all
   * atoms within some cut off of an atom efficiantly.
   *
   * List resolves periodicity using ghost atoms, these are stored and managed internally.
   *
   * An example of using a List to count the average number of atoms within r_cut of each
   * atom:
   */
  class List {
  public:
    /**
     * @brief Construct a new Neighbour List object. The cut off, r_cut, must be smaller than the
     * minimum OrthSimCell extent.
     */
    List(system::Box const& box, double r_cut)
        : m_box(box),
          m_grid(box.make_grid(r_cut)),
          m_cells(visit(m_grid, [](auto const& grid) -> Arr<int> { return grid.shape(); })) {}

    /**
     * @brief Build the internal neighbour lists in parallel with openMP.
     *
     * After a call to this function the for_neighbours methods can be used to iterate over all
     * atoms within r_cut of any atom.
     */
    auto rebuild(system::SoA<Position const&> positions, int num_threads = 1) -> void;

    // /**
    //  * @brief Update the positions of all atoms and ghosts to Position{} -= deltas
    //  *
    //  * Usefull if using a skin distance and wanting to avoid rebuilding the lists.
    //  */
    // void update_positions(Position::matrix_type const& deltas);

    // /**
    //  * @brief Call f(n, r, dr) for every neighbour n of atom i, within distance r_cut.
    //  *
    //  * n is the neighbour index which could be a ghost or real atom, to convert to the index of the
    //  * real atom regardless use .image_to_real(n)
    //  *
    //  * dr is the minimum image vector joining i to n and r is the norm of dr
    //  */
    // template <typename F>
    // void for_neighbours(std::size_t i, floating r_cut, F&& f) const {
    //   //
    //   ASSERT(r_cut <= m_r_cut, "Neighbour lists built with a larger r_cut.");

    //   for (auto&& n : m_neigh_lists[i]) {
    //     Vec3<floating> const dr = m_atoms(Position{}, n) - m_atoms(Position{}, i);
    //     floating const r = norm(dr);
    //     if (r < r_cut) {
    //       f(n, r, dr);
    //     }
    //   }
    // }

    // /**
    //  * @brief Call f(n, dr) for every neighbour n of atom i, within the cut off specified
    //  * during call to rebuild*.
    //  *
    //  * n is the neighbour index which could be a ghost or real atom, to convert to the index of the
    //  * real atom regardless use .image_to_real(n)
    //  *
    //  * dr is the minimum image vector joining i to n and r is the norm of dr
    //  */
    // template <typename F>
    // void for_neighbours(std::size_t i, F&& f) const {
    //   for (auto&& n : m_neigh_lists[i]) {
    //     f(n, m_atoms(Position{}, n) - m_atoms(Position{}, i));
    //   }
    // }

  private:
    struct Next : system::Property<Eigen::Index> {};  ///< Index of the next atom in the linked cell list.

    system::Box m_box;                           ///< Store the box.
    typename system::Box::Grid m_grid;           ///< Store the grid (made by the box).
    AdjacentCells m_cells;                       ///< Store the neighbour cell lists.
    system::SoA<Index, Next, Position> m_atoms;  ///< Store the canonical positions, Index = index-in-input for real+ghosts.
    Eigen::Index m_num_plus_ghosts = 0;          ///< Number of real atoms + number of ghost atoms.
    Vector<Eigen::Index> m_head;                 ///< Index of the start of each bucket.
    Vector<Vector<Eigen::Index>> m_neigh_lists;  ///< Neighbour list for each non-ghost atom.

    static constexpr auto NONE = std::numeric_limits<Eigen::Index>::max();

    /**
     * @brief Initialise memory, copy in atom atoms, build ghosts and build link cell lists.
     */
    void init_and_build_lcl(system::SoA<Position const&> positions, int num_threads);

    /**
     * @brief Set up all ghost indexes positions and offsets.
     *
     * A ghosts position can be calculated from the position of its image plus its offset.
     */
    void make_ghosts();

    /**
     * @brief Build the neighbour list of the ith atom.
     */
    void build_neigh_list(Eigen::Index i);

    // /**
    //  * @brief Convert the neighbour index of a real or ghost atom to the index of the real atom.
    //  *
    //  * This is separate, rather than the n being provided during for_neighbours so as not to pay for
    //  * a cache miss if the real index is not required.
    //  */
    // std::size_t image_to_real(std::size_t i) const { return m_atoms(Index{}, i); }
  };

}  // namespace fly::neighbour