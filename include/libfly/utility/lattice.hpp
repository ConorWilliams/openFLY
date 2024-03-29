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
#include <type_traits>

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
  auto motif_to_lattice(system::Supercell<Map, T...> motif, Arr<int> const &extents)
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

    template_for<int>(Arr<int>::Zero(), extents, [&](Arr<int> const &off) {
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
  auto remove_atoms(system::Supercell<Map, T...> const &cell, std::vector<Eigen::Index> const &bad)
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
  auto add_atoms(system::Supercell<Map, T...> const &cell,
                 std::vector<system::Atom<TypeID, T...>> const &atoms) {
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
  auto remove_sphere(system::Supercell<Map, T...> const &cell, Eigen::Index centre, double r)
      -> system::Supercell<Map, T...> {
    //
    neigh::List nl(cell.box(), r);

    nl.rebuild(cell);

    std::vector<Eigen::Index> bad = {centre};

    nl.for_neighbours(centre, r, [&bad](auto n, auto const &...) { bad.push_back(n); });

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
  inline Vec centroid(system::SoA<Position const &> x) {
    Vec vsum = Vec::Zero();
    for (int i = 0; i < x.size(); i++) {
      vsum += x(r_, i);
    }
    return vsum / x.size();
  }

  /**
   * @brief Translate the atoms in ``x`` such that the centroid of ``x`` is the centroid of ``with``.
   */
  inline void centroid_align(system::SoA<Position &> x, system::SoA<Position const &> with) {
    //
    Vec delta = centroid(with) - centroid(x);

    for (int i = 0; i < x.size(); i++) {
      x(r_, i) += delta;
    }

    ASSERT(gnorm(centroid(with) - centroid(x)) < 1e-10,
           "Norms should be aligned but {}!={}",
           centroid(with),
           centroid(x));
  }

  /**
   * @brief A class that can compare a defective lattice to a perfect lattice and detect the missing atoms.
   */
  class DetectVacancies {
  public:
    /**
     * @brief Construct a new Detect Vacancies object.
     *
     * @param box A description of the simulation space.
     * @param r_lat The "radius" of a lattice point i.e. the minimum distance to be considered at that point.
     * @param perfect_lat The perfect lattice, must consist only of atom of a single type.
     */
    DetectVacancies(double r_lat, system::Box const &box, system::viewSoA<TypeID, Position> perfect_lat);

    /**
     * @brief Find the canonical coordinate of the vacancies.
     *
     * @param lat The defective lattice.
     * @param num_threads The number of (openMP) threads to use.
     */
    std::vector<Vec> detect_vacancies(system::viewSoA<TypeID, Position> lat, int num_threads = 1);

  private:
    TypeID::scalar_t m_tp = 0;

    double m_r_lat;

    system::Box m_box;

    neigh::List m_list;

    system::SoA<Position> m_perfect;

    system::SoA<Position, Index> m_combo;

    Eigen::Index count_type(system::viewSoA<TypeID, Position> lat) const noexcept {
      int tmp = 0;

      for (Eigen::Index i = 0; i < lat.size(); i++) {
        if (lat(id_, i) == m_tp) {
          ++tmp;
        }
      }

      return tmp;
    }
  };

  /**
   * @brief Stores the integers 0... N as a collection of disjoint (non-overlapping) sets..
   */
  class DisjointSetForest {
  public:
    /**
     * @brief Initialise a forest of N sets each containing a single integer 0...N.
     */
    explicit DisjointSetForest(std::size_t N) : m_part(N) {
      for (std::size_t i = 0; i < N; i++) {
        m_part[i].parent = i;
      }
    }

    /**
     * @brief Return an index of the set that the integer ``x``, in 0... N, belongs to.
     *
     * This remains constant until merge() is called.
     */
    std::size_t find(std::size_t x) {
      //
      ASSERT(x < m_part.size(), "{} is not in the forest", x);

      if (m_part[x].parent != x) {
        return (m_part[x].parent = find(m_part[x].parent));  // Path compression
      } else {
        return x;
      }
    }

    /**
     * @brief Replace the set at index ``x`` and the set at index ``y`` with their union.
     */
    void merge(std::size_t x, std::size_t y) {
      //
      ASSERT(x < m_part.size(), "x={} is not in the forest", x);
      ASSERT(y < m_part.size(), "y={} is not in the forest", y);

      if (x == y) {
        return;
      }

      if (m_part[x].size < m_part[y].size) {
        std::swap(x, y);
      }

      m_part[y].parent = x;
      m_part[x].size = m_part[x].size + m_part[y].size;
    }

  private:
    struct Node {
      std::size_t parent;
      std::size_t size = 1;
    };

    std::vector<Node> m_part;
  };

  /**
   * @brief Get the longest edge in a minimum spanning tree of the nodes in ``v``.
   *
   * If there is less than two nodes then this function will return -1.
   *
   * @param v A vector of nodes.
   * @param norm A function which computes the distance between two of the nodes in ``v``.
   */
  template <typename T, typename F>
  auto kruskal_max(std::vector<T> const &v, F const &norm)
      -> std::enable_if_t<std::is_invocable_r_v<double, F const &, T const &, T const &>, double> {
    //

    struct Edge {
      std::size_t A;
      std::size_t B;

      double len;
    };

    std::vector<Edge> adj;

    // Initialise and sort edges

    for (std::size_t i = 0; i < v.size(); i++) {
      for (std::size_t j = 0; j < i; j++) {
        adj.push_back({i, j, std::invoke(norm, v[i], v[j])});
      }
    }

    std::sort(adj.begin(), adj.end(), [](Edge const &x, Edge const &y) { return x.len < y.len; });

    DisjointSetForest forest(v.size());

    double max = 0;
    std::size_t count = 0;

    for (auto const &edge : adj) {
      //
      auto set_A = forest.find(edge.A);
      auto set_B = forest.find(edge.B);

      if (set_A != set_B) {
        //
        forest.merge(set_A, set_B);

        max = std::max(max, edge.len);

        if (++count == v.size() - 1) {
          return max;
        }
      }
    }

    return -1;
  }

}  // namespace fly