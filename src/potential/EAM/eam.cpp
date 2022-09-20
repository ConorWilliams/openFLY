

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

  auto EAM::gradient(system::SoA<TypeID const&, Frozen const&, PotentialGradient&> inout, neigh::List const& nl, int threads) -> void {
  }

  //   std::optional<double> EAM::gradient(SimCell& x, neighbour::List& nl, std::size_t num_threads) {
  //     // Usually a noop
  //     m_aux.destructive_resize(x.size());

  // // First sum computes density  at each atom, runs over active+boundary atoms
  // #pragma omp parallel for num_threads(num_threads) schedule(static)
  //     for (std::size_t b = 0; b < x.size(); b++) {
  //       double rho = 0;

  //       // Computes rho at atom
  //       nl.for_neighbours(b, r_cut(), [&](std::size_t a, double r, Vec3<double> const&) {
  //         rho += m_data->phi(x(id_, nl.image_to_real(a)), x(id_, b)).f(r);
  //       });

  //       // Compute F'(rho) at atom
  //       m_aux(Fprime{}, b) = m_data->f(x(id_, b)).fp(rho);
  //     }

  // // Second sum computes gradient, only runs over active atoms
  // #pragma omp parallel for num_threads(num_threads) schedule(static)
  //     for (std::size_t g = 0; g < x.size(); ++g) {
  //       //
  //       Vec3<double> grad = Vec3<double>::Zero();

  //       if (!x(fzn_, g)) {
  //         nl.for_neighbours(g, r_cut(), [&](std::size_t ap, double r, Vec3<double> const& dr) {
  //           //
  //           std::size_t a = nl.image_to_real(ap);

  //           double mag = m_data->v(x(id_, a), x(id_, g)).fp(r)
  //                        + m_aux(Fprime{}, g) * m_data->phi(x(id_, a), x(id_, g)).fp(r)
  //                        + m_aux(Fprime{}, a) * m_data->phi(x(id_, g), x(id_, a)).fp(r);

  //           grad -= (mag / r) * dr;
  //         });
  //       }

  //       // Write grad to atom
  //       x(Gradient{}, g) = grad;
  //     }

  //     return std::nullopt;
  //   }

  //   void EAM::hessian(SimCell& x, neighbour::List& nl, std::size_t) {
  //     // Usually a noop, make space in aux
  //     m_aux.destructive_resize(x.size());

  //     // Compute hessian indexes
  //     std::size_t count = 0;
  //     for (std::size_t i = 0; i < x.size(); i++) {
  //       if (!x(fzn_, i)) {
  //         m_aux(Hidx{}, i) = count++;
  //       }
  //     }

  //     // First sum computes  rho & mu  at each atom, runs over active + frozen atoms
  //     for (std::size_t b = 0; b < x.size(); ++b) {
  //       double rho = 0;
  //       Vec3<double> mu = Vec3<double>::Zero();
  //       // Compute rho and mu at atom via local sum
  //       nl.for_neighbours(b, r_cut(), [&](std::size_t ap, double r, Vec3<double> const& dr) {
  //         //
  //         std::size_t a = nl.image_to_real(ap);

  //         rho += m_data->phi(x(id_, a), x(id_, b)).f(r);
  //         mu -= m_data->phi(x(id_, a), x(id_, b)).fp(r) * dr / r;  // A.13
  //       });

  //       // Write
  //       m_aux(Rho{}, b) = rho;
  //       m_aux(Mu{}, b) = mu;
  //     }

  //     x.zero_hess();

  //     // Second sums computes hessian, running over all atoms
  //     for (std::size_t z = 0; z < x.size(); ++z) {
  //       // During this section we only write to the z^th column triplet

  //       if (!x(fzn_, z)) {
  //         auto v = m_aux(Mu{}, z).matrix();
  //         // A.15 pre sum term
  //         x.hess(m_aux(Hidx{}, z), m_aux(Hidx{}, z)).matrix() = m_data->f(x(id_, z)).fpp(m_aux(Rho{}, z)) * v *
  //         v.transpose();

  //         // Note: dr = r^{za}
  //         nl.for_neighbours(z, r_cut(), [&](std::size_t ap, double r, Vec3<double> const& dr) {
  //           // First compute sum over neigh for block-diagonal hessian-elements, neighbours can be
  //           // frozen or not

  //           std::size_t a = nl.image_to_real(ap);

  //           // A.14
  //           double A = (m_data->v(x(id_, a), x(id_, z)).fp(r) +

  //                       m_data->f(x(id_, z)).fp(m_aux(Rho{}, z)) * m_data->phi(x(id_, a), x(id_, z)).fp(r) +

  //                       m_data->f(x(id_, a)).fp(m_aux(Rho{}, a)) * m_data->phi(x(id_, z), x(id_, a)).fp(r))
  //                      / r;

  //           // A.14
  //           double B = m_data->v(x(id_, a), x(id_, z)).fpp(r) +

  //                      m_data->f(x(id_, z)).fp(m_aux(Rho{}, z)) * m_data->phi(x(id_, a), x(id_, z)).fpp(r) +

  //                      m_data->f(x(id_, a)).fp(m_aux(Rho{}, a)) * m_data->phi(x(id_, z), x(id_, a)).fpp(r);

  //           double ABFpp = (A - B
  //                           - m_data->f(x(id_, a)).fpp(m_aux(Rho{}, a))
  //                                 * ipow<2>(m_data->phi(x(id_, z), x(id_, a)).fp(r)))
  //                          / (r * r);

  //           x.hess(m_aux(Hidx{}, z), m_aux(Hidx{}, z)).matrix()
  //               += A * Eigen::Matrix<double, spatial_dims, spatial_dims>::Identity() - ABFpp * dr.matrix() *
  //               dr.matrix().transpose();

  //           // Now we will compute the off diagonal element in this column block of H, the a's must
  //           // now not be frozen, we do not commpute overlap here
  //           if (!x(fzn_, a)) {
  //             double BArr = (B - A) / (r * r);

  //             double ur
  //                 = m_data->f(x(id_, z)).fpp(m_aux(Rho{}, z)) * m_data->phi(x(id_, a), x(id_, z)).fp(r) / r;

  //             double ru
  //                 = m_data->f(x(id_, a)).fpp(m_aux(Rho{}, a)) * m_data->phi(x(id_, z), x(id_, a)).fp(r) / r;

  //             // TODO : .noaliase()

  //             x.hess(m_aux(Hidx{}, z), m_aux(Hidx{}, a)).matrix()
  //                 = -BArr * dr.matrix() * dr.matrix().transpose() - A * Eigen::Matrix<double, spatial_dims,
  //                 spatial_dims>::Identity()
  //                   + ur * m_aux(Mu{}, z).matrix() * dr.matrix().transpose() - ru * m_aux(Mu{}, a).matrix() *
  //                   dr.matrix().transpose();
  //           }
  //         });
  //       }
  //     }

  //     // Finally we must compute the overlap terms
  //     for (std::size_t a = 0; a < x.size(); ++a) {
  //       //
  //       double ddFg = m_data->f(x(id_, a)).fpp(m_aux(Rho{}, a));

  //       nl.for_neighbours(a, r_cut(), [&](std::size_t ep, double r_e, auto const& dr_ae) {
  //         nl.for_neighbours(a, r_cut(), [&](std::size_t gp, double r_g, auto const& dr_ag) {
  //           std::size_t e = nl.image_to_real(ep);
  //           std::size_t g = nl.image_to_real(gp);
  //           if (e != g && !x(fzn_, e) && !x(fzn_, g)) {
  //             // Now iterating over all pair of unfrozen neighbours of a
  //             double mag = ddFg / (r_e * r_g) * m_data->phi(x(id_, g), x(id_, a)).fp(r_g)
  //                          * m_data->phi(x(id_, e), x(id_, a)).fp(r_e);

  //             x.hess(m_aux(Hidx{}, e), m_aux(Hidx{}, g)).matrix() += mag * dr_ae.matrix() * dr_ag.matrix().transpose();
  //           }
  //         });
  //       });
  //     }
  //   }
}  // namespace fly::potential