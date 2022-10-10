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
#include "libfly/utility/core.hpp"

/**
 * \file geometry.hpp
 *
 * @brief ...
 *
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
   *    \Delta_\text{RMSD} = \sqrt{\sum_{i=1}^N{\left\lVert Mx_i - y_i \right\rVert^2}}
   *
   * \endrst
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
   *    \Delta_\text{GRMSD} = \sqrt{\sum_{i=1}^N{\left\lVert x_i - y_i \right\rVert^2}}
   *
   * \endrst
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
   * \endrst
   *
   * Uses the "Kabsch algorithm" see: https://en.wikipedia.org/wiki/Kabsch_algorithm
   */
  template <typename... M1, typename... M2>
  auto ortho_onto(system::VoS<M1...> const &x, system::VoS<M2...> const &y) -> Mat {
    //
    Mat H = Mat::Zero();

    ASSERT(x.size() == y.size(), "Sizes must match, {}!={}!", x.size(), y.size());

    // ASSERT(near(gnorm(centroid(x)), 0.0), "Centroid of x={}", centroid(x));
    // ASSERT(near(gnorm(centroid(y)), 0.0), "Centroid of y={}", centroid(y));

    for (int i = 0; i < x.size(); i++) {
      H.noalias() += x[i][r_] * y[i][r_].transpose();
    }

    Eigen::JacobiSVD<Mat> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

    // Don't do sign correction as it's ok to reflect the coordinate system.

    return svd.matrixV() * svd.matrixU().transpose();
  }

  //   namespace detail {

  //     /** @brief The maximum number of atoms that can lie in the same plane. */
  //     inline constexpr std::size_t MAX_COPLANAR_ATOMS = 10;

  //     /**
  //      * @brief Helper function, verifies addition of n^th atom matches all previous atoms.
  //      */
  //     template <typename... M1, typename... M2>
  //     bool within_tol_up_to(system::VoS<M1...> const &mut, system::VoS<M2...> const &ref, double tol, std::size_t n) {
  //       //
  //       for (std::size_t i = 0; i < std::min(n, MAX_COPLANAR_ATOMS); ++i) {
  //         // Intra atomic distances
  //         double ref_ni = norm(ref[n][r_] - ref[i][r_]);
  //         double mut_ni = norm(mut[n][r_] - mut[i][r_]);

  //         if (std::abs(ref_ni - mut_ni) > tol) {
  //           return false;
  //         }
  //       }
  //       return true;
  //     }

  //     /**
  //      * @brief Recursive implementation of for_equiv_perms.
  //      *
  //      * Uses the GREEDY_PERM algorithm, recursively finds a new atom, such that the order that
  //      * the intra-atom distances in each system::VoS match (within tolerance SQRT_2 * delta). The
  //      * final permutation must be able to be transformed (via ortho_onto(mut, other)) such that the
  //      * RMSD between the two sets is less than delta.
  //      */
  //     template <typename... M1, typename... M2, typename F>
  //     bool for_equiv_perms_impl(system::VoS<M1...> &mut, system::VoS<M2...> const &ref, double delta, std::size_t n, F const &f) {
  //       // Termination criterion.
  //       if (n >= mut.size()) {
  //         //
  //         Mat O = ortho_onto(mut, ref);

  //         if (double dr = rmsd(O, mut, ref); dr < delta) {
  //           return std::invoke(f, O, dr);
  //         } else {
  //           return false;
  //         }
  //       }

  //       using std::swap;  // ADL

  //       // Attempt to find an atom to put in i^th position.
  //       for (std::size_t i = n; i < ref.size(); ++i) {
  //         if (mut[i](Colour{}) == ref[n](Colour{})) {
  //           // Test next candidate
  //           swap(mut[n], mut[i]);

  //           // Verify all other distances then recurse
  //           if (within_tol_up_to(mut, ref, delta * M_SQRT2, n)) {
  //             if (for_equiv_perms_impl(mut, ref, delta, n + 1, f)) {
  //               // Enable early exit.
  //               return true;
  //             }
  //           }

  //           // Undo swap such that if we have to back-track order is the same.
  //           swap(mut[n], mut[i]);
  //         }
  //       }

  //       return false;
  //     }

  //   }  // namespace detail

  //   /**
  //    * @brief Create a function object to explore equivalent permutations of another system::VoS.
  //    *
  //    * This function is for greedily exploring the permutations of atoms in another system::VoS such
  //    * that, after an orthogonal transformation generated by ortho_onto(), the RMSD between the two
  //    * sets of atoms is less than delta.
  //    *
  //    * Example:
  //    *
  //    * @code{.cpp}
  //    *
  //    * for_equiv_perms(ref, 0.2)(1, mut, [](Mat const & O, double RMSD){
  //    *     // If this function is called then "mut" has been permuted into an equivalent permutation.
  //    *     // "O" is the matrix that maps "mut" to "ref" such that that:
  //    *     //     "RMSD" == RMSD(O, mut, ref) and "RMSD" < delta.
  //    *
  //    *     // We could now do something with "RMSD" and "O". We must not mutate "ref" or "mut".
  //    *
  //    *     // If we "return true" then exploration/function terminates.
  //    *     // If we "return false" then exploration/function continues.
  //    * });
  //    *
  //    * @endcode
  //    *
  //    * for_equiv_perms will not permute the first "n" (== 1 in the above example) atoms in "mut". This
  //    * is useful if there exist a number of atoms who's permutation is known.
  //    *
  //    */
  //   template <typename... M>
  //   auto for_equiv_perms(system::VoS<M...> const &ref, double delta) {
  //     return [&ref, delta](std::size_t n, auto &mut, auto &&callback) {
  //       ASSERT(ref.size() == mut.size(), "Sizes must match!");

  //       if constexpr (std::is_same_v<std::remove_cv_t<decltype(mut)>, system::VoS<M...>>) {
  //         ASSERT(&ref != &mut, "Cannot perm onto self.");
  //       }

  //       detail::for_equiv_perms_impl(mut, ref, delta, n, std::forward<decltype(callback)>(callback));
  //     };
  //   }

  //   /**
  //    * @brief A geometry models a local distribution of periodically resolved Atoms.
  //    *
  //    * The distribution of atoms is centred on the first atom.
  //    */
  //   template <typename... Mems>
  //   class Geometry : public system::VoS<Mems...> {
  //   public:
  //     /**
  //      * @brief Returned by permute_onto(delta, x), contains extra info about the permutation.
  //      */
  //     struct PermResult {
  //       /** @brief  Orthogonal transformation required to map this onto x. */
  //       Mat O;
  //       /** @brief Equal to rmsd(O, *this, x). */
  //       double rmsd;
  //     };

  //     /**
  //      * @brief Attempt to permute the atoms in this Geometry into the same order as the atoms in
  //      * other.
  //      *
  //      * @return std::optional<PermResult> An engaged PermResult if a permutation was found.
  //      */
  //     template <typename... Ts>
  //     std::optional<PermResult> permute_onto(system::VoS<Ts...> const &other, double delta) {
  //       std::optional<PermResult> res = std::nullopt;

  //       for_equiv_perms(other, delta)(1, *this, [&](Mat const &O, double rmsd) {
  //         res = PermResult{O, rmsd};
  //         // First match accepted
  //         return true;
  //       });

  //       return res;
  //     }

  //     /**
  //      * @brief Find the best permutation of the atoms in this geomety onto other
  //      *
  //      * @return std::optional<PermResult> An engaged PermResult if a permutation was found.
  //      */
  //     template <typename... Ts>
  //     std::optional<PermResult> best_perm_onto(system::VoS<Ts...> const &other, double delta, Geometry<Mems...> &scratch) {
  //       std::optional<PermResult> res = std::nullopt;

  //       for_equiv_perms(other, delta)(1, *this, [&](Mat const &O, double rmsd) {
  //         if (!res || res->rmsd < rmsd) {
  //           res = PermResult{O, rmsd};
  //           scratch = *this;
  //         }

  //         return false;
  //       });

  //       if (res) {
  //         *this = scratch;
  //       }

  //       return res;
  //     }
  //   };

}  // namespace fly::env