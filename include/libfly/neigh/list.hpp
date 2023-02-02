#pragma once

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

#include <cmath>
#include <limits>

#include "libfly/neigh/adjacent.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

// Forward decl for friendship
namespace fly::potential {
  class KIM_API;
}

/**
 * \file list.hpp
 *
 * @brief Classes to build and update neighbour lists.
 */

namespace fly::neigh {

  /**
   * @brief A class to contain, build and manage neighbour lists in shared memory.
   *
   * Neighbour-lists are used to find the atoms within some cut-off of all atoms efficiently. The cache
   * efficiency of all List's major operations can be enhances by pre-sorting the atoms according to their
   * grid index, this can be done with the fly::neigh::sort() function.
   *
   * Designed with the intention of being reused; List separates the building, updating and using of the
   * neighbour-lists, performed by List::rebuild(), List::update() and List::for_neighbours() respectively.
   *
   * List resolves periodicity using *ghost* atoms, these are stored and managed internally.
   *
   * \rst
   *
   * .. admonition:: Implementation notes
   *
   *    Construction of a neighbour-list occurs as follows:
   *
   *    1. Real atoms are mapped into the canonical cell.
   *    2. Ghost atoms are generated.
   *    3. All atoms are linked into lists of atoms in the same grid-cell.
   *    4. Iterating over all adjacent grid-cells the neighbour-lists or real atoms are built.
   *
   *    The neighbour-lists contain the indexes of the real or ghost atoms. Alongside each atom we store the
   *    index of the real atom it may be an image of, this allows us to map ghost atoms to real atoms.
   *
   * \endrst
   *
   */
  class List {
  public:
    /**
     * @brief Control the  maximum number of ghosts.
     *
     * List supports up to MAX_GHOST_RATIO * List::size() ghost atoms.
     */
    static constexpr Eigen::Index MAX_GHOST_RATIO = 26;

    /**
     * @brief Construct a new List object.
     *
     * @param box The description of the simulation space.
     * @param r_cut The neighbour cut-off radius.
     */
    List(system::Box const& box, double r_cut)
        : m_box(box),
          m_grid(box.make_grid(r_cut)),
          m_cells(visit(m_grid, [](auto const& grid) { return grid.shape(); })),
          m_r_cut(r_cut) {}

    /**
     * @brief Builds the internal neighbour-lists.
     *
     * @param positions The positions of the atoms.
     * @param num_threads The number of (openMP) threads spawned to complete this operation.
     */
    auto rebuild(system::SoA<Position const&> positions, int num_threads = 1) -> void;

    /**
     * @brief Update the positions of the real + ghost atoms.
     *
     * \rst
     *
     * If the positions of the real atoms are stored in the :math:`3N \times 1` vector :math:`r` then after
     * calling this function:
     *
     * .. math::
     *
     *    r \gets r - x
     *
     * The ghost atoms are then updated accordingly.
     *
     * \endrst
     *
     * @param x The change in the positions of the real atoms.
     */
    auto update(system::SoA<Delta const&> x) -> void;

    /**
     * @brief Call ``f(n, r, dr)`` for every neighbour of atom ``i`` within cut-off ``r_cut``.
     *
     * \pre List::rebuild() must have be called.
     *
     * \rst
     *
     * An example to count the average number of neighbours:
     *
     * .. include:: ../../examples/neigh/list.cpp
     *    :code:
     *
     * \endrst
     *
     * @param i The index of the atom whose neighbours will be iterated over.
     * @param r_cut The neighbour cut-off.
     * @param f The functor with signature ``f`` : ``Eigen::Index``, ``double``, ``fly::Vec`` -> ``void``.
     */
    template <typename F>
    auto for_neighbours(Eigen::Index i, double r_cut, F&& f) const -> void {
      //
      if (r_cut > m_r_cut) {
        throw error("r_cut={} is bigger then than the internal r_cut={}", r_cut, m_r_cut);
      }

      ASSERT(i >= 0 && i < size(), "Atom {} is not in neighbour list length {}", i, size());

      for (auto&& n : m_neigh_lists[i]) {
        Vec dr = m_atoms(r_, n) - m_atoms(r_, i);
        double r_sq = gnorm_sq(dr);
        if (r_sq < r_cut * r_cut) {
          f(image_to_real(n), std::sqrt(r_sq), dr);
        }
      }
    }

    /**
     * @brief Get the number of real atoms in this neighbour list.
     *
     * @return The number of real atoms in this List.
     */
    auto size() const noexcept -> Eigen::Index { return m_neigh_lists.size(); }

  private:
    friend class ::fly::potential::KIM_API;

    struct Next : system::Property<Eigen::Index> {};  ///< Index of the next atom in the linked cell list.
    struct Contrib : system::Property<int> {};        /// To support the KIM-API

    using soa_type = system::SoA<Index, Next, Position, Contrib>;

    system::Box m_box;                           ///< Store the box.
    typename system::Box::Grid m_grid;           ///< Store the grid (made by the box).
    detail::AdjacentCells m_cells;               ///< Store the neighbour cell lists.
    double m_r_cut;                              ///< Store the cut of radius.
    soa_type m_atoms;                            ///< Canonical pos, Index = index-in-input for real+ghosts.
    Eigen::Index m_num_plus_ghosts = 0;          ///< Number of real atoms + number of ghost atoms.
    Vector<Eigen::Index> m_head;                 ///< Index of the start of each bucket.
    Vector<Vector<Eigen::Index>> m_neigh_lists;  ///< Neighbour list for each non-ghost atom.

    static constexpr auto NONE = std::numeric_limits<Eigen::Index>::max();

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

    /**
     * @brief Convert the neighbour index of a real or ghost atom to the index of the real atom.
     *
     * Very careful to avoid any side effects here (i.e. throwing) such that if the result is not needed this
     * function can be optimised away.
     */
    Eigen::Index image_to_real(Eigen::Index i) const noexcept {
      return *(m_atoms[Index{}].data() + Index::size() * i);
    }

    static int get_cluster_neigh(void* const dataObject,
                                 int const numberOfNeighborLists,
                                 double const* const cutoffs,
                                 int const neighborListIndex,
                                 int const particleNumber,
                                 int* const numberOfNeighbors,
                                 int const** const neighborsOfParticle);
  };

}  // namespace fly::neigh