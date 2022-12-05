#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "libfly/neigh/list.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file lattice.hpp
 *
 * @brief Utilities for building and manipulating basic lattices of atoms.
 */
namespace fly {

  /**
   * @brief Utility to build a supercell by replicating a ``motif`` along each axis ``extents`` times.
   *
   * @param motif Containing the ``TypeMap`` and positions of the atoms - in fractional coordinates - the box
   * is interpreted as the primitive cell.
   * @param extents Number of primitive-cells along each axis of the supercell.
   */
  template <typename Map, typename... T>
  auto motif_to_lattice(system::Supercell<Map, T...> motif, Arr<int> const& extents)
      -> system::Supercell<Map, T...> {
    //
    verify((extents > 0).all(), "Invalid extents={}", extents);

    Mat cell = motif.box().basis();

    // Transform motif's fractional coordinates to real coordinates:

    for (int i = 0; i < motif.size(); i++) {
      if ((motif(r_, i).array() < 0).any() || (motif(r_, i).array() >= 1).any()) {
        throw error("Atom #{} in motif is outside unit cell, r={}, are you using fractional coordinates?",
                    i,
                    motif(r_, i));
      }
      motif(r_, i) = cell * motif(r_, i);
    }

    // Construct empty supercell

    Mat super_basis = cell;
    //
    for (int i = 0; i < spatial_dims; i++) {
      super_basis.col(i) *= extents[i];
    }

    Arr<bool> periodic;
    //
    for (int i = 0; i < spatial_dims; i++) {
      if (!(periodic[i] = motif.box().periodic(i))) {
        verify(extents[i] == 1, "Non-periodic extent index {} must be 1 not {}", i, extents[i]);
      }
    }

    system::Supercell supercell
        = system::make_supercell<T...>({super_basis, periodic}, motif.map(), extents.prod() * motif.size());

    // Fill supercell

    int i = 0;

    template_for<int>(Arr<int>::Zero(), extents, [&](Arr<int> const& off) {
      //
      Vec super_off = cell * off.matrix().cast<double>();

      for (int j = 0; j < motif.size(); j++) {
        supercell(id_, i) = motif(id_, j);
        ((supercell(T{}, i) = motif(T{}, j)), ...);

        supercell(r_, i) += super_off;

        i++;
      }
    });

    supercell[r_] += M_PI;  // This irrational number removes atoms on the boundary

    return supercell;
  }

  /**
   * @brief Create a new supercell identical to ``cell`` but with the atoms in ``bad`` removed.
   *
   * @param cell Input ``SuperCell``.
   * @param bad Indexes of atoms to remove.
   */
  template <typename Map, typename... T>
  auto remove_atoms(system::Supercell<Map, T...> const& cell, std::vector<Eigen::Index> const& bad)
      -> system::Supercell<Map, T...> {
    //

    auto out_cell = system::make_supercell<T...>(
        cell.box(), cell.map(), cell.size() - safe_cast<Eigen::Index>(bad.size()));

    Eigen::Index x = 0;

    for (Eigen::Index i = 0; i < cell.size(); i++) {
      if (auto it = std::find(bad.begin(), bad.end(), i); it == bad.end()) {
        out_cell(id_, x) = cell(id_, i);
        ((out_cell(T{}, x) = cell(T{}, i)), ...);
        x++;
      }
    }

    return out_cell;
  }

  /**
   * @brief Add atoms to a supercell.
   *
   * @param cell Input ``SuperCell``.
   * @param atoms Atoms to add to ``cell``.
   */
  template <typename Map, typename... T>
  auto add_atoms(system::Supercell<Map, T...> const& cell,
                 std::vector<system::Atom<TypeID, T...>> const& atoms) {
    //
    auto out_cell = system::make_supercell<T...>(
        cell.box(), cell.map(), cell.size() + safe_cast<Eigen::Index>(atoms.size()));

    // Copy old atoms.
    for (Eigen::Index i = 0; i < cell.size(); i++) {
      out_cell(id_, i) = cell(id_, i);
      ((out_cell(T{}, i) = cell(T{}, i)), ...);
    }

    // Copy new atoms.
    for (std::size_t i = 0; i < atoms.size(); i++) {
      Eigen::Index x = cell.size() + safe_cast<Eigen::Index>(i);

      out_cell(id_, x) = atoms[i][id_];
      ((out_cell(T{}, x) = atoms[i][T{}]), ...);
    }

    return out_cell;
  }

  /**
   * @brief Create a new supercell identical to ``cell`` but with the atoms within ``r`` of the atom
   * ``centre`` removed.
   *
   * @param cell Input ``SuperCell``.
   * @param centre Index of central atom.
   * @param r Distance from ``centre`` to qualify for removal.
   */
  template <typename Map, typename... T>
  auto remove_sphere(system::Supercell<Map, T...> const& cell, Eigen::Index centre, double r)
      -> system::Supercell<Map, T...> {
    //
    neigh::List nl(cell.box(), r);

    nl.rebuild(cell);

    std::vector<Eigen::Index> bad = {centre};

    nl.for_neighbours(centre, r, [&bad](auto n, auto const&...) { bad.push_back(n); });

    return remove_atoms(cell, bad);
  }

  /**
   * @brief Compute the centroid of a vector of atoms, ``x``.
   *
   * \rst
   *
   * For a set of atoms with positions :math:`x_i` this is defined as:
   *
   * .. math::
   *    \frac{1}{N} \sum_{i=1}^{N}{x_i}
   *
   * \endrst
   */
  inline Vec centroid(system::SoA<Position const&> x) {
    Vec vsum = Vec::Zero();
    for (int i = 0; i < x.size(); i++) {
      vsum += x(r_, i);
    }
    return vsum / x.size();
  }

}  // namespace fly