

#include "libfly/potential/ROT/dimer.hpp"

#include <cmath>
#include <cstddef>
#include <memory>

#include "libfly/minimise/LBFGS/lbfgs_core.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::potential {

  Dimer::Dimer(Options const& opt, potential::Generic const& to_wrap)
      : m_opt(opt), m_core(opt.n), m_wrapped(std::make_unique<potential::Generic>(to_wrap)) {}

  Dimer::Dimer(Dimer const& other)
      : m_opt{other.m_opt}, m_core(m_opt.n), m_wrapped(std::make_unique<potential::Generic>(*other.m_wrapped)) {}

  Dimer::Dimer(Dimer&& other) noexcept : m_opt{other.m_opt}, m_core(std::move(other.m_core)), m_wrapped(std::move(other.m_wrapped)) {}

  double Dimer::r_cut() const noexcept { return m_wrapped->r_cut() + m_opt.delta_r * 2; }

  double Dimer::eff_gradient(system::SoA<PotentialGradient&, Axis&> out,
                             system::SoA<TypeID const&, Frozen const&, Axis const&> in,
                             neigh::List& nl,
                             int num_threads) {
    //
    verify(in.size() == out.size(), "Effective gradient size mismatch in={} out={}", in.size(), out.size());

    m_delta.destructive_resize(in.size());       // Store displacement for updates
    m_delta_prev.destructive_resize(in.size());  // Store previous displacement for updates

    m_axisp.destructive_resize(in.size());  // Temporary axis

    m_g0.destructive_resize(in.size());   // Central grad
    m_g1.destructive_resize(in.size());   // Temp end grad
    m_g1p.destructive_resize(in.size());  // Temp primed end grad

    m_delta_g.destructive_resize(in.size());  // Gradient difference

    out[ax_] = in[ax_];

    //
    m_core.clear();

    m_wrapped->gradient(m_g0, in, nl, num_threads);  // Gradient at centre (g0)

    m_delta[del_] = -m_opt.delta_r * out[ax_];

    nl.update(m_delta);

    using std::swap;  // ADL

    swap(m_delta, m_delta_prev);

    m_wrapped->gradient(m_g1, in, nl, num_threads);  // Gradient at end (g1)

    double curv = [&] {
      for (int i = 0;; i++) {
        m_delta_g[del_] = m_g1[g_] - m_g0[g_];
        m_delta_g[del_] -= gdot(m_delta_g[del_], out[ax_]) * out[ax_];  // Torque

        // Use lbfgs to find rotation plane
        auto& theta = m_core.newton_step<Axis, Delta>(out, m_delta_g);

        theta[del_] -= gdot(theta[del_], out[ax_]) * out[ax_];  // Ortho
        theta[del_] *= 1 / gnorm(theta[del_]);                  //      normalization

        double b_1 = gdot(m_g1[g_] - m_g0[g_], theta[del_]) / m_opt.delta_r;
        double c_x0 = gdot(m_g1[g_] - m_g0[g_], out[ax_]) / m_opt.delta_r;
        double theta_1 = -0.5 * std::atan(b_1 / std::abs(c_x0));  // Trial rotation angle

        if (std::abs(theta_1) < m_opt.theta_tol || i == m_opt.iter_max_rot) {
          return c_x0;
        } else {
          // Trial rotation
          m_axisp[ax_] = out[ax_] * std::cos(theta_1) + theta[del_] * std::sin(theta_1);

          m_delta[del_] = -m_opt.delta_r * m_axisp[ax_];  // Temporarily store next m_delta_prev into m_delta
          swap(m_delta, m_delta_prev);                    // Now put it in the correct place
          m_delta[del_] = m_delta_prev[del_] - m_delta[del_];
          nl.update(m_delta);

          m_wrapped->gradient(m_g1p, in, nl, num_threads);  // Gradient at primed end (g1p)

          double c_x1 = gdot(m_g1p[g_] - m_g0[g_], m_axisp[ax_]) / m_opt.delta_r;
          double a_1 = (c_x0 - c_x1 + b_1 * sin(2 * theta_1)) / (1 - std::cos(2 * theta_1));
          double theta_min = 0.5 * std::atan(b_1 / a_1);  // Optimal rotation

          // Flip if extrema is maxima
          if (a_1 * std::cos(2 * theta_min) - a_1 + b_1 * std::sin(2 * theta_min) > 0) {
            theta_min += M_PI / 2;
          }

          out[ax_] = out[ax_] * std::cos(theta_min) + theta[del_] * std::sin(theta_min);

          // Interpolate force at new rotation
          m_g1[g_] = (std::sin(theta_1 - theta_min) / std::sin(theta_1)) * m_g1[g_]
                     + (std::sin(theta_min) / std::sin(theta_1)) * m_g1p[g_]
                     + (1 - std::cos(theta_min) - std::sin(theta_min) * std::tan(0.5 * theta_1)) * m_g0[g_];

          if (m_opt.debug) {
            fmt::print("\tDimer: i={:<4} theta={:f} curv={:f}\n",
                       i,
                       std::abs(theta_min),
                       gdot(m_g1[g_] - m_g0[g_], out[ax_]) / m_opt.delta_r);
          }

          if (std::abs(theta_min) < m_opt.theta_tol) {
            return gdot(m_g1[g_] - m_g0[g_], out[ax_]) / m_opt.delta_r;
          }
        }
      }
    }();

    m_delta[del_] = -m_delta_prev[del_];  // Flip to reverse last.
    nl.update(m_delta);                   // Reset neighbour lists to input state.

    if (!m_opt.relax_in_convex && curv > 0) {
      out[g_] = -gdot(m_g0[g_], out[ax_]) * out[ax_];
    } else {
      out[g_] = m_g0[g_] - 2 * gdot(m_g0[g_], out[ax_]) * out[ax_];
    }

    return curv;
  }

}  // namespace fly::potential