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
#include <memory>

#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file kim.hpp
 *
 * @brief This file contains a potential class for interfacing with KIM models.
 */

// Forward declarations of KIM API classes.
namespace KIM {

  class Model;

  class ComputeArguments;

}  // namespace KIM

namespace fly::potential {

  /**
   * @brief This class provides an interface to KIM potential models (https://openkim.org/).
   */
  class KIM_API {
  public:
    struct Options {
      std::string model_name = "No default";  ///< The exact name of the KIM model.
      double epsilon = 1e-2;                  ///< Hessian finite-difference step size.
      bool debug = false;                     ///< Print debug information.
    };

    explicit KIM_API(Options const& opt, system::TypeMap<> const& map);

    ~KIM_API() noexcept;

    /**
     * @brief Get this potentials cut-off radius.
     *
     * This is the maximum distance two atom can interact. The neighbour::List passed to the other functions
     * should be configured with a cut-off equal or greater than this.
     */
    auto r_cut() const noexcept -> double;

    /**
     * @brief Compute the potential energy.
     *
     * Assumes the neighbour list are ready.
     *
     * @param in Per-atom data used by potential for computation.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called)
     * configured with a cut-off at least ``r_cut()``.
     * @param threads Dummy parameter, only here to match the interface of other potentials.
     * @return The potential energy of the system of atoms.
     */
    auto energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int threads = 1)
        -> double;

    /**
     * @brief Compute potential energy gradient.
     *
     * Assumes the neighbour list are ready, force on frozen atoms will be zero.
     *
     * @param in Per-atom TypeID's and Frozen properties
     * @param out Potential gradient written to this.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild()
     * called).
     * @param threads Dummy parameter, only here to match the interface of other potentials.
     */
    auto gradient(system::SoA<PotentialGradient&> out,
                  system::SoA<TypeID const&, Frozen const&> in,
                  neigh::List const& nl,
                  int threads = -1) -> void;

  private:
    Options m_opt;

    struct KIM_code : system::Property<int> {};
    struct KIM_contrib : system::Property<int> {};
    struct KIM_energy : system::Property<double> {};
    struct KIM_force : system::Property<double, 3> {};

    system::TypeMap<KIM_code> m_map;  ///< Maps TypeID to KIM code.

    system::SoA<KIM_code, KIM_contrib, KIM_energy, KIM_force> m_kim_io;

    KIM::Model* m_model = nullptr;
    KIM::ComputeArguments* m_args = nullptr;
  };

}  // namespace fly::potential