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

#include "KIM/meam_c.hpp"
#include "KIM/meam_spline.hpp"
#include "KIM/meam_sw_spline.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/hessian.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file meam.hpp
 *
 * @brief MEAM potential implementation.
 */

namespace fly::potential {

  /**
   * @brief MEAM potential child class.
   */
  class MEAM {
  public:
    /**
     * @brief Construct a new MEAM object.
     *
     * @param map TypeMap of expected types to find in data.
     */
    MEAM(system::TypeMap<> const& map, std::string const&) {
      {
        std::FILE* parameter_file_pointers[kMaxNumberParameterFiles];

        int number_parameter_files(0);

        // Getting the number of parameter files, open the files and process them
        model_driver_create->GetNumberOfParameterFiles(&number_parameter_files);

        if (number_parameter_files > kMaxNumberParameterFiles) {
          LOG_ERROR("Too many input parameter files!\n");
          *ier = true;
          return;
        }

        if (!number_parameter_files) {
          LOG_ERROR("There is no parameter file!\n");
          *ier = true;
          return;
        }

        *ier = OpenParameterFiles(model_driver_create, number_parameter_files, parameter_file_pointers);
        if (*ier) {
          return;
        }

        *ier = ProcessParameterFiles(model_driver_create, number_parameter_files, parameter_file_pointers);

        CloseParameterFiles(number_parameter_files, parameter_file_pointers);

        if (*ier) {
          return;
        }
      }
    }

  private:
    /*!
     * \brief Flag indicating we do not need to request neighbors from
     *        non-contributing particles for this model driver
     */
    static constexpr int model_will_not_request_neighbors_of_non_contributing_particles_{1};

    /*!
     * \brief Number of particles
     *
     * \note
     * This is a Mutable value that can change with each
     * call to \b Refresh() and \b Compute().
     */
    int cached_number_of_particles_{0};

    /*! # of unique elements */
    int number_of_elements_{0};

    /*! meam/c style flag */
    int is_meam_c_{0};

    /*! meam/spline style flag */
    int is_meam_spline_{0};

    /*! meam/sw/spline style flag */
    int is_meam_sw_spline_{0};

    /*! max cutoff for all elements */
    double max_cutoff_{0.0};

    /*! max cutoff squared for all elements */
    double max_cutoff_squared_{0.0};

    /*!
     * \brief Cutoff value in %KIM API object
     *
     * \note
     * This is a Mutable value that can change with each
     * call to \b Refresh() and \b Compute().
     */
    double influence_distance_{0.0};

    /*! names of unique elements */
    std::vector<std::string> element_name_;

    /*! MEAMC object */
    std::unique_ptr<openKIM::MEAMC> meam_c_;

    /*! MEAMSpline object */
    std::unique_ptr<openKIM::MEAMSpline> meam_spline_;

    /*! MEAMSWSpline object */
    std::unique_ptr<openKIM::MEAMSWSpline> meam_sw_spline_;
  };

}  // namespace fly::potential