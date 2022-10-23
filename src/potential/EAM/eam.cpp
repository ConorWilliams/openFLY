// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/potential/EAM/eam.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>

#include "fmt/core.h"
#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::potential {

  double EAM::energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int num_threads) {
    double v_sum = 0;
    double f_sum = 0;

#pragma omp parallel for reduction(+ : v_sum, f_sum) num_threads(num_threads) schedule(static)
    for (Eigen::Index a = 0; a < in.size(); a++) {
      // Skip contribution from frozen atoms.
      if (!in(fzn_, a)) {
        double rho = 0;

        nl.for_neighbours(a, r_cut(), [&](auto b, double r, Vec const&) {
          //
          v_sum += m_data->v(in(id_, a), in(id_, b)).f(r);
          rho += m_data->phi(in(id_, b), in(id_, a)).f(r);
        });

        f_sum += m_data->f(in(id_, a)).f(rho);
      }
    }

    return (0.5 * v_sum) + f_sum;
  }

  void EAM::gradient(system::SoA<PotentialGradient&> out,
                     system::SoA<TypeID const&, Frozen const&> in,
                     neigh::List const& nl,
                     int num_threads) {
    //
    verify(in.size() == out.size(), "EAM gradient size mismatch in={} out={}", in.size(), out.size());

    // Usually a noop
    m_aux.destructive_resize(in.size());

// First sum computes density  at each real atom, runs over frozen+active atoms.
#pragma omp parallel for num_threads(num_threads) schedule(static)
    // Compute rho at all
    for (Eigen::Index b = 0; b < in.size(); b++) {
      double rho = 0;

      // Computes rho at atom
      nl.for_neighbours(b, r_cut(), [&](auto a, double r, Vec const&) {
        //
        rho += m_data->phi(in(id_, a), in(id_, b)).f(r);
      });

      // Compute F'(rho) at atom
      m_aux(Fprime{}, b) = m_data->f(in(id_, b)).fp(rho);
    }

// Second sum computes gradient, only runs over active atoms
#pragma omp parallel for num_threads(num_threads) schedule(static)
    //
    for (Eigen::Index g = 0; g < in.size(); ++g) {
      //
      Vec grad = Vec::Zero();

      if (!in(fzn_, g)) {
        nl.for_neighbours(g, r_cut(), [&](auto a, double r, Vec const& dr) {
          //

          double mag = m_data->v(in(id_, a), in(id_, g)).fp(r) + m_aux(Fprime{}, g) * m_data->phi(in(id_, a), in(id_, g)).fp(r)
                       + m_aux(Fprime{}, a) * m_data->phi(in(id_, g), in(id_, a)).fp(r);

          grad -= (mag / r) * dr;
        });
      }

      // Write grad to atom
      out(g_, g) = grad;
    }
  }

  void EAM::hessian(system::Hessian& out, system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int num_threads) {
    // Usually a noop, make space in aux
    m_aux.destructive_resize(in.size());

    // Compute hessian indexes and count active
    Eigen::Index num_active = 0;

    for (Eigen::Index i = 0; i < in.size(); i++) {
      if (!in(fzn_, i)) {
        m_aux(Hidx{}, i) = num_active++;
      }
    }

    out.zero_for(num_active);

// First sum computes  rho & mu  at each atom, runs over active + frozen atoms
#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (Eigen::Index b = 0; b < in.size(); ++b) {
      double rho = 0;
      Vec mu = Vec::Zero();
      // Compute rho and mu at atom via local sum
      nl.for_neighbours(b, r_cut(), [&](auto a, double r, Vec const& dr) {
        rho += m_data->phi(in(id_, a), in(id_, b)).f(r);
        mu -= m_data->phi(in(id_, a), in(id_, b)).fp(r) * dr / r;  // A.13
      });

      // Write
      m_aux(Rho{}, b) = rho;
      m_aux(Mu{}, b) = mu;
    }

    // Second sums computes hessian, running over all atoms, only writing to z^th column block. As responsible for LOWER z^th column
    // use dynamic scheduling.
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (Eigen::Index z = 0; z < in.size(); ++z) {
      // During this section we only write to the z^th column block

      if (!in(fzn_, z)) {
        //
        Vec v = m_aux(Mu{}, z);
        // A.15 pre sum term, including off-diagonal elements along block diagonal make code simpler.
        out(m_aux(Hidx{}, z), m_aux(Hidx{}, z)).noalias() = m_data->f(in(id_, z)).fpp(m_aux(Rho{}, z)) * v * v.transpose();

        // Note: dr = r^{za}
        nl.for_neighbours(z, r_cut(), [&](auto a, double r, Vec const& dr) {
          // First compute sum over neigh for block-diagonal hessian-elements, neighbours can be
          // frozen or not

          // A.14
          double A = (m_data->v(in(id_, a), in(id_, z)).fp(r) +

                      m_data->f(in(id_, z)).fp(m_aux(Rho{}, z)) * m_data->phi(in(id_, a), in(id_, z)).fp(r) +

                      m_data->f(in(id_, a)).fp(m_aux(Rho{}, a)) * m_data->phi(in(id_, z), in(id_, a)).fp(r))
                     / r;

          // A.14
          double B = m_data->v(in(id_, a), in(id_, z)).fpp(r) +

                     m_data->f(in(id_, z)).fp(m_aux(Rho{}, z)) * m_data->phi(in(id_, a), in(id_, z)).fpp(r) +

                     m_data->f(in(id_, a)).fp(m_aux(Rho{}, a)) * m_data->phi(in(id_, z), in(id_, a)).fpp(r);

          double ABFpp
              = (A - B - m_data->f(in(id_, a)).fpp(m_aux(Rho{}, a)) * ipow<2>(m_data->phi(in(id_, z), in(id_, a)).fp(r))) / (r * r);

          out(m_aux(Hidx{}, z), m_aux(Hidx{}, z)) += A * Mat::Identity() - ABFpp * dr * dr.transpose();

          // Now we will compute the off diagonal element in this column block of H, the a's must now not be frozen, we do not
          // compute overlap here. Only need to compute lower triangular part of hessian hence, drop writes to out(a, z) with a < z.
          ASSERT(z != a, "Atoms {} is interacting with itself", z);

          if (a >= z && !in(fzn_, a)) {
            double BArr = (B - A) / (r * r);

            double ur = m_data->f(in(id_, z)).fpp(m_aux(Rho{}, z)) * m_data->phi(in(id_, a), in(id_, z)).fp(r) / r;

            double ru = m_data->f(in(id_, a)).fpp(m_aux(Rho{}, a)) * m_data->phi(in(id_, z), in(id_, a)).fp(r) / r;

            out(m_aux(Hidx{}, a), m_aux(Hidx{}, z)).noalias() += -BArr * dr * dr.transpose() - A * Mat::Identity()
                                                                 + ur * m_aux(Mu{}, z) * dr.transpose()
                                                                 - ru * m_aux(Mu{}, a) * dr.transpose();
          }

          // Finally compute overlap terms of z coupling to active atom g via a. Each thread only writes to e^th column blocks.

          double ddFg = m_data->f(in(id_, a)).fpp(m_aux(Rho{}, a)) * m_data->phi(in(id_, z), in(id_, a)).fp(r);

          nl.for_neighbours(a, r_cut(), [&](auto g, double r_g, Vec const& dr_ag) {
            if (g > z && !in(fzn_, g)) {
              // Now iterating over all pair of unfrozen neighbours of a
              double mag = ddFg / (r * r_g) * m_data->phi(in(id_, g), in(id_, a)).fp(r_g);

              out(m_aux(Hidx{}, g), m_aux(Hidx{}, z)).noalias() -= mag * dr * dr_ag.transpose();
            }
          });
        });
      }
    }

// Transform into mass-weighted hessian
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (Eigen::Index i = 0; i < in.size(); ++i) {
      if (!in(fzn_, i)) {
        //
        double mass_i = m_data->type_map().get(in(id_, i), m_);

        for (Eigen::Index j = i; j < in.size(); ++j) {
          if (!in(fzn_, j)) {
            //
            double mass_j = m_data->type_map().get(in(id_, j), m_);

            out(m_aux(Hidx{}, j), m_aux(Hidx{}, i)) *= 1 / std::sqrt(mass_i * mass_j);
          }
        }
      }
    }
  }

}  // namespace fly::potential