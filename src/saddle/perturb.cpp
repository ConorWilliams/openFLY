// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <algorithm>
#include <random>
#include <vector>

#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/random.hpp"

namespace fly::saddle {

  /**
   * @brief Canonicalize centre and return it and its 3^spatial_dims - 1 canonical images.
   *
   * @param centre
   * @param box
   * @return std::vector<Position::matrix_t>
   */
  std::vector<Position::matrix_t> gen_images(Position::matrix_t const& centre, system::Box const& box) {
    //
    std::vector<Position::matrix_t> images = {box.canon_image(centre)};

    Mat basis = box.basis();

    template_for<int>(Arr<int>::Constant(-1), Arr<int>::Constant(2), [&](auto... args) {
      //
      Arr<int> signs{args...};

      Vec offset = Vec::Zero();

      for (int i = 0; i < spatial_dims; i++) {
        if (box.periodic(i)) {  // Only in periodic axis
          offset += basis.col(i) * signs(i);
        }
      }

      if ((signs == 0).all() || offset != Vec::Zero()) {
        images.push_back(images[0] + offset);
      }
    });

    return images;
  }

  void perturb(system::SoA<Position&, Axis&> out,
               Xoshiro& urbg,
               system::Box const& box,
               Position::matrix_t const& centre,
               system::SoA<Position const&, Frozen const&> cell,
               double rcut,
               double stddev) {
    // PRNG

    std::normal_distribution<double> normal(0, 1);
    std::normal_distribution<double> gauss(0, stddev);

    out[r_] = cell[r_];
    out[ax_] = 0;  // Zero the rotation axis.

    std::vector c_image = gen_images(centre, box);

    for (int i = 0; i < cell.size(); i++) {
      //
      if (!cell(Frozen{}, i)) {
        //
        Position::matrix_t tmp = box.canon_image(cell(r_, i));

        auto it = std::min_element(c_image.begin(), c_image.end(), [&](auto const& a, auto const& b) {
          //
          return gnorm_sq(tmp - a) < gnorm_sq(tmp - b);
        });

        //
        double dr2 = gnorm_sq(tmp - *it);

        if (dr2 < rcut * rcut) {
          //
          double env = 1 - std::sqrt(dr2) / rcut;

          out(r_, i) += env * Vec{gauss(urbg), gauss(urbg), gauss(urbg)};

          out(ax_, i) += Vec{
              normal(urbg),
              normal(urbg),
              normal(urbg),
          };
        }
      }
    }

    out[ax_] /= gnorm(out[ax_]);  // normalize
  }

}  // namespace fly::saddle