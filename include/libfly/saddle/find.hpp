// #pragma once

// // Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// // SPDX-License-Identifier: GPL-3.0-or-later

// // This file is part of openFLY.

// // OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// // as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// // OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// // warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// // You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

// #include <vector>

// #include "libfly/minimise/LBFGS/lbfgs.hpp"
// #include "libfly/potential/ROT/dimer.hpp"
// #include "libfly/potential/generic.hpp"
// #include "libfly/system/SoA.hpp"
// #include "libfly/system/box.hpp"
// #include "libfly/system/property.hpp"
// #include "libfly/utility/core.hpp"
// #include "libfly/utility/random.hpp"

// /**
//  * \file find.hpp
//  *
//  * @brief Utility for perturbing a supercell.
//  */

// namespace fly::saddle {

//   /**
//    * @brief a
//    *
//    */
//   class MasterFinder {
//   public:
//     /**
//      * @brief Construct a new MasterFinder object.
//      *
//      * @param min
//      * @param pot
//      * @param opt
//      * @param num_threads
//      */
//     MasterFinder(minimise::LBFGS const& min, potential::Generic const& pot, potential::Dimer::Options const& opt, int num_threads)
//         : m_data(safe_cast<std::size_t>(num_threads),
//                  ThreadData{
//                      min,
//                      pot,
//                      potential::Generic{
//                          // Construct a dimer wrapping a copy of the input potential.
//                          potential::Dimer{
//                              opt,
//                              pot,
//                          },
//                      },
//                  }) {}

//     // template <typename... T>
//     // void find_all(std::vector<int> const& unknown, system::SoA<T...> const& in) {

//     // }

//   private:
//     struct ThreadData {
//       minimise::LBFGS min;
//       potential::Generic pot;
//       potential::Generic dimer;
//     };

//     std::vector<ThreadData> m_data;
//   };

// }  // namespace fly::saddle