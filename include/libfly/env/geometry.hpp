#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.centroid>

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

#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "libfly/system/VoS.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file geometry.hpp
 *
 * @brief Utilities for manipulating local environments of atoms.
 */

namespace fly::env {
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
   *
   * @param x Set of atoms with the Position property.
   */
  template <typename... M>
  auto centroid(system::VoS<M...> const &x) noexcept -> Vec {
    //
    Vec sum = Vec::Zero();

    for (auto const &elem : x) {
      sum += elem[r_];
    }

    return sum / x.size();
  }

  /**
   * @brief Compute the generalised RMSD between two sets of atoms.
   *
   * \rst
   *
   * The generalised RMSD (root mean squared distance) between two sets of atoms, ``x`` and ``y``, is the RMSD between them after
   * applying the matrix transformation ``M`` to each atom in ``x``. For two sets of atoms with positions :math:`x_i` and :math:`y_i`
   * this is:
   *
   * .. math::
   *    \Delta_\text{GRMSD} = \sqrt{\sum_{i=1}^N{\left\lVert Mx_i - y_i \right\rVert^2}}
   *
   * \endrst
   *
   * @param x Set of atoms with the Position property.
   * @param y Set of atoms with the Position property.
   * @param M Transformation matrix.
   */
  template <typename E, typename... M1, typename... M2>
  auto grmsd(E const &M, system::VoS<M1...> const &x, system::VoS<M2...> const &y) -> double {
    //
    ASSERT(x.size() == y.size(), "Sizes must match, {}!={}!", x.size(), y.size());

    double sum_sq = 0;

    for (int i = 0; i < x.size(); ++i) {
      sum_sq += gnorm_sq(y[i][r_] - M * x[i][r_]);
    }

    return std::sqrt(sum_sq);
  }

  /**
   * @brief Compute the RMSD between the positions of the atoms in x and y.
   *
   * \rst
   *
   * The RMSD (root mean squared distance) is the sqrt of the sum of the squared distances between
   * atoms, for two sets of atoms with positions :math:`x_i` and :math:`y_i` this is:
   *
   * .. math::
   *    \Delta_\text{RMSD} = \sqrt{\sum_{i=1}^N{\left\lVert x_i - y_i \right\rVert^2}}
   *
   * \endrst
   *
   * @param x Set of atoms with the Position property.
   * @param y Set of atoms with the Position property.
   */
  template <typename... M1, typename... M2>
  auto rmsd(system::VoS<M1...> const &x, system::VoS<M2...> const &y) -> double {
    return grmsd(Mat::Identity(), x, y);
  }

  /**
   * @brief Compute the orthogonal transformation matrix that best transforms ``x`` onto ``y``.
   *
   * This does not permute any of the atoms it just minimizes the RMSD between the two sets.
   *
   * \rst
   *
   * Effectively solves this optimisation problem:
   *
   * .. math::
   *    \min_{O}{\sum_{i=1}^N{\left\lVert O x_i - y_i \right\rVert^2}} \quad \text{with} \quad OO^\intercal = I
   *
   * With :math:`x_i` and :math:`y_i` the position of the ``i`` th atom in ``x`` onto ``y`` respectively.
   *
   * \endrst
   *
   * Uses the "Kabsch algorithm" see: https://en.wikipedia.org/wiki/Kabsch_algorithm
   *
   * @param x Set of atoms with the Position property.
   * @param y Set of atoms with the Position property.
   */
  template <typename... M1, typename... M2>
  auto ortho_onto(system::VoS<M1...> const &x, system::VoS<M2...> const &y) -> Mat {
    //
    Mat H = Mat::Zero();

    ASSERT(x.size() == y.size(), "Sizes must match, {}!={}!", x.size(), y.size());

    for (int i = 0; i < x.size(); i++) {
      H.noalias() += x[i][r_] * y[i][r_].transpose();
    }

    Eigen::JacobiSVD<Mat> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

    // Don't do sign correction as it's ok to reflect the coordinate system.

    return svd.matrixV() * svd.matrixU().transpose();
  }

  namespace detail {

    /** @brief The maximum number of atoms that can lie in the same plane. */
    inline constexpr int MAX_COPLANAR_ATOMS = 10;

    /**
     * @brief Helper function, verifies addition of n^th atom matches all previous atoms.
     */
    template <typename... M1, typename... M2>
    bool within_tol_up_to(system::VoS<M1...> const &mut, system::VoS<M2...> const &ref, double tol, int n) {
      //
      for (int i = 0; i < std::min(n, MAX_COPLANAR_ATOMS); ++i) {
        // Intra atomic distances
        double ref_ni = gnorm(ref[n][r_] - ref[i][r_]);
        double mut_ni = gnorm(mut[n][r_] - mut[i][r_]);

        if (std::abs(ref_ni - mut_ni) > tol) {
          return false;
        }
      }
      return true;
    }

    /**
     * @brief Recursive implementation of for_equiv_perms.
     *
     * Uses the GREEDY_PERM algorithm, recursively finds a new atom, such that the order that
     * the intra-atom distances in each system::VoS match (within tolerance SQRT_2 * delta). The
     * final permutation must be able to be transformed (via ortho_onto(mut, other)) such that the
     * RMSD between the two sets is less than delta.
     */
    template <typename... M1, typename... M2, typename F>
    bool for_equiv_perms_impl(system::VoS<M1...> &mut, system::VoS<M2...> const &ref, double delta, int n, F const &f) {
      // Termination criterion.
      if (n >= mut.size()) {
        //
        Mat O = ortho_onto(mut, ref);

        if (double dr = grmsd(O, mut, ref); dr < delta) {
          return std::invoke(f, O, dr);
        } else {
          return false;
        }
      }

      using std::swap;  // ADL

      // Attempt to find an atom to put in i^th position.
      for (int i = n; i < ref.size(); ++i) {
        if (mut[i][col_] == ref[n][col_]) {
          // Test next candidate
          swap(mut[n], mut[i]);

          // Verify all other distances then recurse
          if (within_tol_up_to(mut, ref, delta * M_SQRT2, n) && for_equiv_perms_impl(mut, ref, delta, n + 1, f)) {
            // Enable early exit.
            return true;
          }

          // Undo swap such that if we have to back-track order is the same.
          swap(mut[n], mut[i]);
        }
      }

      return false;
    }

  }  // namespace detail

  /**
   * @brief Explore equivalent permutations of atoms in two sets.
   *
   * This function is for greedily exploring the permutations of atoms in ``mut`` such that, after an orthogonal transformation
   * generated by ``ortho_onto()``, the ``grmsd()`` between the two sets of atoms is less than ``delta``.
   *
   *
   * \rst
   *
   * This effectively finds a matrix :math:`O` and permutation :math:`\pi` that satisfies:
   *
   * .. math::
   *    \sum_{i=1}^N{\left\lVert O x_{\pi(i)} - y_i \right\rVert^2} < \delta^2
   *
   * Subject to the constraints:
   *
   * .. math::
   *    OO^\intercal = I \quad \text{and} \quad \pi(i) = i \,\,\,\, \forall \,\,\,\,  i < n
   *
   * With :math:`x_i` and :math:`y_i` the position of the ``i`` th atom in ``x`` onto ``y`` respectively.
   *
   * Example:
   *
   * .. include:: ../../examples/env/geometry.cpp
   *    :code:
   *
   * \endrst
   *
   * ``for_equiv_perms()`` will not permute the first ``n`` (== 1 in the above example) atoms in ``mut``. This
   * is useful if there exist a number of atoms who's permutation is known.
   *
   *
   * @param mut Set of atoms with the Position & Colour properties that will be permuted to match ``ref``.
   * @param ref The reference set of atoms with the Position & Colour properties.
   * @param delta The maximum ``grmsd()`` between the permuted ``mut`` and ``ref``
   * @param n The first ``n`` atoms in ``mut`` will not be permuted.
   * @param f An invokable invoked with signature ``(fly::Mat const& O, double rmsd) -> bool`` called for each permutation of ``mut``
   * that matches ``ref``, see the example above for details.
   */
  template <typename... M1, typename... M2, typename F>
  auto for_equiv_perms(system::VoS<M1...> &mut, system::VoS<M2...> const &ref, double delta, int n, F const &f) -> void {
    //
    ASSERT(ref.size() == mut.size(), "Sizes must match, {}!={}!", ref.size(), mut.size());

    if constexpr (std::is_same_v<system::VoS<M1...>, system::VoS<M2...>>) {
      ASSERT(&ref != &mut, "Cannot perm onto self.", 0);
    }

    detail::for_equiv_perms_impl(mut, ref, delta, n, f);
  }

  /**
   * @brief A Geometry models a local distribution of periodically resolved Atoms.
   *
   * The distribution of atoms is centred on the first atom, i.e. the first atom cannot move and is always fixed during permutation
   * operations.
   *
   * A Geometry is a ``fly::system::VoS``  with the default properties Position and Colour.
   */
  template <typename... Pr>
  class Geometry : public system::VoS<Position, Colour, Pr...> {
  public:
    /**
     * @brief Returned by permutation methods, contains extra info about the permutation.
     */
    struct GeoInfo {
      Mat O;        ///< Orthogonal transformation required to map this onto ``x``.
      double rmsd;  ///<  Equal to ``grmsd(O, *this, x)``.
    };

    /**
     * @brief Find the first permutation of the atoms in this Geometry onto ``other`` within tolerance ``delta``.
     *
     * @param other The reference set of atoms with the Position & Colour properties.
     * @param delta The maximum ``grmsd()`` between the permuted ``this`` and ``ref`` if this operation succeeds.
     * @return std::optional<GeoInfo> An engaged ``GeoInfo`` if a permutation was found.
     */
    template <typename... Ts>
    std::optional<GeoInfo> permute_onto(system::VoS<Ts...> const &other, double delta) {
      std::optional<GeoInfo> res = std::nullopt;

      for_equiv_perms(*this, other, delta, 1, [&res](Mat const &O, double rmsd) {
        res = GeoInfo{O, rmsd};
        return true;  // First match accepted
      });

      return res;
    }

    /**
     * @brief Find the best permutation of the atoms in this Geometry onto ``other``.
     *
     * @param other The reference set of atoms with the Position & Colour properties.
     * @param delta The maximum ``grmsd()`` between the permuted ``this`` and ``ref`` for this operation to succeed.
     * @param scratch_best An optional pointer to scratch space that the algorithm can use (otherwise this function will allocate).
     * @return std::optional<GeoInfo> An engaged ``GeoInfo`` if a permutation was found.
     */
    template <typename... Ts>
    std::optional<GeoInfo> best_perm_onto(system::VoS<Ts...> const &other, double delta, Geometry *scratch_best = nullptr) {
      //
      std::optional<GeoInfo> res = std::nullopt;

      Geometry tmp;

      if (!scratch_best) {
        scratch_best = &tmp;
      }

      for_equiv_perms(*this, other, delta, 1, [&](Mat const &O, double rmsd) {
        if (!res || res->rmsd < rmsd) {
          res = GeoInfo{O, rmsd};  // Cache best info.
          scratch_best = *this;    // Copy best so far into scratch.
        }

        return false;  // Explore all.
      });

      if (res) {
        using std::swap;
        swap(*this, scratch_best);
      }

      return res;
    }
  };

}  // namespace fly::env