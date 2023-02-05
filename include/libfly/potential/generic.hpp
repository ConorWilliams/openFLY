#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlooK.com>

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
 * The potential completely describes the interatomic interactions present in a system.
 */

namespace fly::potential {

  /**
   * @brief Generalised interatomic potential.
   *
   * Essentially a ``std::variant`` of all supported potentials.
   */
  class Generic {
  private:
    using variant = std::variant<EAM>;

  public:
    /**
     * @brief Construct a new Generic object
     *
     * @param args Forwarded to ``std::variant``'s constructor.
     */
    template <typename... Args, typename = std::enable_if_t<std::is_constructible_v<variant, Args&&...>>>
    explicit Generic(Args&&... args) : m_pot(std::forward<Args>(args)...) {}

    /**
     * @brief Construct a new Generic object.
     */
    Generic(Generic const&) = default;

    /**
     * @brief Construct a new Generic object.
     */
    Generic(Generic&&) = default;

    /**
     * @brief Get this potentials cut-off radius.
     *
     * This is the maximum distance two atoms can interact. The neighbour::List passed to the other functions
     * should be configured with a cut-off equal or greater than this.
     */
    auto r_cut() const noexcept -> double {
      return ::fly::visit(m_pot, [](auto const& pot) -> double { return pot.r_cut(); });
    }

    /**
     * @brief Compute the potential energy.
     *
     * Can assume the neighbour list are ready.
     *
     * @param in Per-atom data used by potential for computation.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called)
     * configured with a cut-off at least ``r_cut()``.
     * @param threads Number of openMP threads to use.
     * @return The potential energy of the system of atoms.
     */
    auto energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int threads = 1)
        -> double {
      return ::fly::visit(m_pot, [&, threads](auto& pot) -> double { return pot.energy(in, nl, threads); });
    }

    /**
     * @brief Compute potential energy gradient.
     *
     * Can assume the neighbour list are ready, force on frozen atoms must be zero.
     *
     * @param out Result is written here.
     * @param in Per-atom data used by potential for computation.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called)
     * configured with a cut-off at least ``r_cut()``.
     * @param threads Number of openMP threads to use.
     */
    auto gradient(system::SoA<PotentialGradient&> out,
                  system::SoA<TypeID const&, Frozen const&> in,
                  neigh::List const& nl,
                  int threads = 1) -> void {
      return ::fly::visit(m_pot, [&, threads](auto& pot) { pot.gradient(out, in, nl, threads); });
    }

    /**
     * @brief Compute hessian matrix.
     *
     * Can assume the neighbour list are ready. The resulting hessian must be n by n (n = number of atoms *
     * spatial dims) and only include contributions from the m active atoms i.e. have zeros for frozen atoms.
     * As hessian matrices are always symmetric this function is only required to compute the lower diagonal
     * portion.
     *
     * @param in Per-atom data used by hessian for computation.
     * @param out Hessian matrix to write output to.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called)
     * configured with a cut-off at least ``r_cut()``.
     * @param threads Number of openMP threads to use.
     */
    auto hessian(system::Hessian& out,
                 system::SoA<TypeID const&, Frozen const&> in,
                 neigh::List const& nl,
                 int threads = 1) -> void {
      return ::fly::visit(m_pot, [&](auto& pot) { pot.hessian(out, in, nl, threads); });
    }

  private:
    variant m_pot;
  };

}  // namespace fly::potential