// // Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// // SPDX-License-Identifier: GPL-3.0-or-later

// // This file is part of openFLY.

// // OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// // as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// // OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// // warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// // You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

// #include "libatom/minimise/LBFGS/lbfgs.hpp"

// #include <fmt/core.h>
// #include <fmt/os.h>

// #include <cmath>
// #include <cstddef>
// #include <optional>
// #include <utility>

// #include "libatom/asserts.hpp"
// #include "libatom/atom.hpp"
// #include "libatom/io/xyz.hpp"
// #include "libatom/minimise/LBFGS/core.hpp"
// #include "libatom/neighbour/list.hpp"
// #include "libatom/potentials/base.hpp"
// #include "libatom/sim_cell.hpp"
// #include "libatom/utils.hpp"

// namespace otf::minimise {

//   bool LBFGS::minimise(SimCell &atoms, potentials::Base &pot, std::size_t num_threads) {
//     //
//     // Clear history from previous runs;

//     m_core.clear();

//     floating skin = std::max(std::pow(m_opt.skin_frac, 1. / 3.) - 1, 0.0) * pot.rcut();

//     if (m_opt.debug) {
//       fmt::print("Skin = {}\n", skin);
//     }

//     m_nl = neighbour::List(atoms, pot.rcut() + skin);

//     m_nl->rebuild(atoms, num_threads);

//     pot.gradient(atoms, *m_nl, num_threads);

//     floating trust = m_opt.min_trust;

//     auto file = [&]() -> std::optional<fmt::ostream> {
//       if (m_opt.debug) {
//         return fmt::output_file("lbfgs_debug.xyz");
//       } else {
//         return std::nullopt;
//       }
//     }();

//     floating acc = 0;
//     floating convex_count = 0;

//     for (std::size_t i = 0; i < m_opt.iter_max; ++i) {
//       //
//       floating mag_g = norm_sq(atoms(Gradient{}));

//       if (m_opt.debug) {
//         constexpr auto str = "LBFGS: i={:<4} trust={:f} acc={:f} norm(g)={:e}\n";
//         fmt::print(str, i, trust, acc, std::sqrt(mag_g));
//         io::dump_xyz(*file, atoms, "Debug");
//       }

//       if (mag_g < m_opt.f2norm * m_opt.f2norm) {
//         return true;
//       } else if (convex_count >= m_opt.convex_max) {
//         return false;
//       }

//       Position::matrix_type &Hg = m_core.newton_step(atoms(Position{}), atoms(Gradient{}));

//       ASSERT(gdot(atoms(Gradient{}), Hg) > 0, "Ascent direction");

//       // Limit step size.
//       Hg *= std::min(1.0, trust / norm(Hg));

//       // Add distance of most displaced atom
//       acc += std::sqrt((Hg * Hg).colwise().sum().maxCoeff());

//       // Update positions in real space
//       atoms(Position{}) -= Hg;

//       if (acc > 0.5 * skin) {
//         m_nl->rebuild(atoms, num_threads);
//         acc = 0;
//       } else {
//         m_nl->update_positions(Hg);
//       }

//       if (std::optional curv = pot.gradient(atoms, *m_nl, num_threads); curv > 0) {
//         convex_count += 1;
//       } else {
//         convex_count = 0;
//       }

//       floating proj = gdot(atoms(Gradient{}), Hg);

//       if (proj < -m_opt.proj_tol) {
//         trust = std::max(m_opt.min_trust, m_opt.shrink_trust * trust);
//       } else if (proj > m_opt.proj_tol) {
//         trust = std::min(m_opt.max_trust, m_opt.grow_trust * trust);
//       }
//     }

//     return false;
//   }

// }  // namespace otf::minimise
