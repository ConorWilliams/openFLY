#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <cstddef>
#include <exception>
#include <memory>
#include <type_traits>
#include <utility>
#include <variant>

#include "libfly/neigh/list.hpp"
#include "libfly/potential/EAM/eam.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file generic.hpp
 *
 * @brief Generic potential energy class.
 *
 * In libFLY potentials have concrete implementations to cut-down compile times. Potentials get info about atoms via a ``SoA``. This
 * allows potentials to slice generic supercells while remaining concrete. However, a traditional virtual interface for potentials
 * would overly constrain the input to potentials (e.g what about a potential that needs a different per-atom property e.g. charge).
 * Hence, the generic potential interface accepts generic ``SoA``s, if the ``SoA`` can be statically sliced to the potential's desired
 * input then this computation occurs otherwise, an error is thrown at run-time. This delaying of errors until run-time enables a
 * single variant-like class to store a set of potentials that do not share a common interface.
 */

namespace fly::potential {

  /**
   * @brief Generalised interatomic potential.
   *
   * Essentially a ``std::variant`` of all supported potentials.
   */
  class Generic {
  private:
    // clang-format off

    template <typename Potential, typename SOA>
    using Energy = decltype(
        std::declval<Potential>().energy(
            std::declval<SOA>(), std::declval<neigh::List const&>(), std::declval<int>()
        )
    );

    template <typename Potential, typename SOA>
    using Gradient = decltype(
        std::declval<Potential>().gradient(
            std::declval<SOA>(), std::declval<neigh::List const&>(), std::declval<int>()
        )
    );

    template <typename Potential, typename SOA>
    using Hessian = decltype(
        std::declval<Potential>().hessian(
            std::declval<SOA>(), std::declval<system::Hessian &>(), std::declval<neigh::List const&>(), std::declval<int>()
        )
    );

    // clang-format on
  public:
    /**
     * @brief Construct a new Generic object
     *
     * @param args Forwarded to ``std::variant``'s constructor.
     */
    template <typename... Args>
    explicit Generic(Args&&... args) : m_pot(std::forward<Args>(args)...) {}

    /**
     * @brief Get this potentials cut-off radius.
     *
     * This is the maximum distance two atoms can interact. The neighbour::List passed to the other functions should be configured with
     * a cut-off equal or greater than this.
     */
    auto r_cut() const noexcept -> double {
      return visit(m_pot, [](auto const& pot) -> double { return pot.r_cut(); });
    }

    /**
     * @brief Compute the potential energy.
     *
     * Can assume the neighbour list are ready.
     *
     * @param in Per-atom data used by potential for computation.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     * @return double The potential energy of the system of atoms.
     */
    template <typename... Ts>
    auto energy(system::SoA<Ts...> const& in, neigh::List const& nl, int threads = 1) -> double {
      return visit(m_pot, [&, threads](auto& pot) -> double {
        if constexpr (is_detected_v<Energy, decltype(pot), system::SoA<Ts...> const&>) {
          return pot.energy(in, nl, threads);
        } else {
          throw error("Generic potential {}, does not support .energy(...) of this system.", m_pot.index());
        }
      });
    }

    /**
     * @brief Compute potential energy gradient.
     *
     * Can assume the neighbour list are ready, force on frozen atoms must be zero.
     *
     * @param inout Per-atom data used by potential for computation. Result is written to the PotentialGradient property.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     */
    template <typename... Ts>
    auto gradient(system::SoA<Ts...>& inout, neigh::List const& nl, int threads = 1) -> void {
      return visit(m_pot, [&](auto& pot) {
        if constexpr (is_detected_v<Energy, decltype(pot), system::SoA<Ts...>&>) {
          pot.gradient(inout, nl, threads);
        } else {
          throw error("Generic potential {}, does not support .gradient(...) of this system.", m_pot.index());
        }
      });
    }

    /**
     * @brief Compute hessian matrix..
     *
     * Can assume the neighbour list are ready. The resulting hessian must be m by m and only include contributions from the m active
     * atoms. As hessian matrices are always symmetric this function is only required to compute the lower diagonal portion.
     *
     * @param in er-atom data used by potential for computation.
     * @param out Hessian matrix to write output to.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     */
    template <typename... Ts>
    auto hessian(system::SoA<Ts...> const& in, system::Hessian& out, neigh::List const& nl, int threads = 1) -> void {
      return visit(m_pot, [&](auto& pot) {
        if constexpr (is_detected_v<Hessian, decltype(pot), system::SoA<Ts...>&>) {
          pot.hessian(in, out, nl, threads);
        } else {
          throw error("Generic potential {}, does not support .hessian(...) of this system.", m_pot.index());
        }
      });
    }

  private:
    std::variant<EAM> m_pot;
  };

}  // namespace fly::potential