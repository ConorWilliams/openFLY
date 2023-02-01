//
// meam_sw_spline.hpp
//
// LGPL Version 2.1 HEADER START
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
//
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301  USA
//
// LGPL Version 2.1 HEADER END
//

//
// Copyright (c) 2020--2021, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//
// Brief: This file is adapted from the LAMMPS software package
//        `lammps/src/USER-MISC/pair_meam_sw_spline.h`
//        `lammps/src/USER-MISC/pair_meam_sw_spline.cpp`
//        and it is rewritten and updated for the KIM-API by
//        Yaser Afshar
//

#ifndef MEAM_SW_SPLINE_HPP
#define MEAM_SW_SPLINE_HPP

#pragma GCC system_header

#include <cstddef>
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "spline.hpp"
#include "utility.hpp"

#ifdef MEAM_SW_SPLINE_LOG_ERROR
#undef MEAM_SW_SPLINE_LOG_ERROR
#endif

namespace openKIM {

/*!
 * \brief Helper macro for printing error message
 * formats messages, filename, line number and function
 * name into an std::ostringstream object
 *
 */
#define MEAM_SW_SPLINE_LOG_ERROR(msg)                                         \
  {                                                                           \
    std::ostringstream ss;                                                    \
    ss << "\nError :" << __FILE__ << ":" << __LINE__ << ":@(" << __FUNCTION__ \
       << ")\n"                                                               \
       << msg << "\n\n";                                                      \
    std::cerr << ss.str();                                                    \
  }

/*! \class MEAMSWSpline
 * \brief Modified Embedded Atom Method with an additional
 *        Stillinger-Weber (SW) term.
 *
 * The MEAMSWSpline class is a base class to compute pairwise interactions for
 * metals using a variant of modified embedded-atom method (MEAM) potentials
 * with an additional Stillinger-Weber (SW) term (Stillinger) in the energy.
 *
 */
class MEAMSWSpline {
 public:
  /*!
   * \brief Grow local array if necessary
   *
   * \param nall Total number of atoms contributing and non-contributing
   */
  inline void ResizeLocalArray(std::size_t const nall);

  /*!
   * \brief Parse the input potential file and process the read parameters
   *
   * \param potential_file_pointer FILE pointer to the opened potential file
   * \param max_line_size Maximum size of the line
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  inline int ProcessPotentialFile(std::FILE *const potential_file_pointer,
                                  int const max_line_size);

  /*!
   * \brief Convert units of the parameters
   *
   * \param convert_length_factor length unit conversion factor
   * \param convert_energy_factor energy unit conversion factor
   */
  inline void ConvertUnit(double const convert_length_factor,
                          double const convert_energy_factor);

  /*!
   * \brief Complete the set up
   *
   * \param max_cutoff
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  inline int CompleteSetup(double *max_cutoff);

 public:
  /*! Flag indicating whether we should use a regular grid. */
  int use_regular_grid_;

  /*! Shift embedding energy by this value to make
   *  it zero for a single atom in vacuum. */
  double zero_atom_energy_;

  /*! Used for temporary storage of U'(rho) values */
  std::vector<double> e_u_prime_;

  /*! Phi_i(r_ij) */
  Spline e_phi_;
  /*! U_i(rho) */
  Spline e_u_;

  /*! Rho_ij(r_ij) */
  Spline rho_r_;
  /*! f_i(r_ij) */
  Spline rho_f_;
  /*! g_ij(cos_theta) */
  Spline rho_g_;

  /*! F(r_ij) */
  Spline esw_f_;
  /*! G(cos_theta) */
  Spline esw_g_;
};

inline void MEAMSWSpline::ResizeLocalArray(std::size_t const nall) {
  // grow local arrays if necessary
  if (nall > e_u_prime_.size()) {
    // 2**14
    constexpr std::size_t kDelta = 16384;

    // Maximum number of atoms
    std::size_t const nmax = nall / kDelta * kDelta + kDelta;

    e_u_prime_.resize(nmax);
  }
}

inline int MEAMSWSpline::ProcessPotentialFile(
    std::FILE *const potential_file_pointer, int const max_line_size) {
  std::rewind(potential_file_pointer);

  std::unique_ptr<char> next_line_(new char[max_line_size]);
  char *next_line = next_line_.get();

  // Create a utility object
  Utility ut;

  // Skip the first line. The first line is always a comment
  if (ut.GetFirstLine(potential_file_pointer, next_line, max_line_size)) {
    std::string msg = "End of file while reading the first line ";
    msg += "from the `meam/sw/spline` potential file.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  // Read and parse the input file
  if (e_phi_.ProcessParameterFile(potential_file_pointer, max_line_size)) {
    std::string msg = "Unable to read and parse the `meam/sw/spline` ";
    msg += "potential file for e_phi_.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  use_regular_grid_ = e_phi_.AreKnotsOnRegularGrid();

  if (esw_f_.ProcessParameterFile(potential_file_pointer, max_line_size)) {
    std::string msg = "Unable to read and parse the `meam/sw/spline` ";
    msg += "potential file for esw_f_.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (use_regular_grid_ != esw_f_.AreKnotsOnRegularGrid()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (esw_g_.ProcessParameterFile(potential_file_pointer, max_line_size)) {
    std::string msg = "Unable to read and parse the `meam/sw/spline` ";
    msg += "potential file for esw_g_.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (use_regular_grid_ != esw_g_.AreKnotsOnRegularGrid()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (rho_r_.ProcessParameterFile(potential_file_pointer, max_line_size)) {
    std::string msg = "Unable to read and parse the `meam/sw/spline` ";
    msg += "potential file for rho_r_.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (use_regular_grid_ != rho_r_.AreKnotsOnRegularGrid()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (e_u_.ProcessParameterFile(potential_file_pointer, max_line_size)) {
    std::string msg = "Unable to read and parse the `meam/sw/spline` ";
    msg += "potential file for e_u_.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (use_regular_grid_ != e_u_.AreKnotsOnRegularGrid()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (rho_f_.ProcessParameterFile(potential_file_pointer, max_line_size)) {
    std::string msg = "Unable to read and parse the `meam/sw/spline` ";
    msg += "potential file for rho_f_.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (use_regular_grid_ != rho_f_.AreKnotsOnRegularGrid()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (rho_g_.ProcessParameterFile(potential_file_pointer, max_line_size)) {
    std::string msg = "Unable to read and parse the `meam/sw/spline` ";
    msg += "potential file for rho_g_.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  if (use_regular_grid_ != rho_g_.AreKnotsOnRegularGrid()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  // Everything is good
  return false;
}

inline void MEAMSWSpline::ConvertUnit(double const convert_length_factor,
                                      double const convert_energy_factor) {
  e_phi_.ConvertUnit(convert_length_factor, convert_energy_factor);
  e_u_.ConvertUnit(1.0, convert_energy_factor);

  rho_r_.ConvertUnit(convert_length_factor, 1.0);
  rho_f_.ConvertUnit(convert_length_factor, 1.0);

  esw_f_.ConvertUnit(convert_length_factor, 1.0);
}

inline int MEAMSWSpline::CompleteSetup(double *max_cutoff) {
  if (e_phi_.UpdateSpline()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }
  if (e_u_.UpdateSpline()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }
  if (rho_r_.UpdateSpline()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }
  if (rho_f_.UpdateSpline()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }
  if (rho_g_.UpdateSpline()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }
  if (esw_f_.UpdateSpline()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }
  if (esw_g_.UpdateSpline()) {
    std::string msg = "This driver does not support both uniform & ";
    msg += "non-uniform cubic splines at the same time.\n";
    MEAM_SW_SPLINE_LOG_ERROR(msg);
    return true;
  }

  // Calculate 'zero-point energy' of a single atom in vacuum.
  zero_atom_energy_ =
      use_regular_grid_ ? e_u_.Eval<true>(0.0) : e_u_.Eval<false>(0.0);

  // Determine maximum cutoff radius of all relevant spline functions.
  *max_cutoff = 0.0;

  if (e_phi_.GetCutoff() > *max_cutoff) {
    *max_cutoff = e_phi_.GetCutoff();
  }

  if (rho_r_.GetCutoff() > *max_cutoff) {
    *max_cutoff = rho_r_.GetCutoff();
  }

  if (rho_f_.GetCutoff() > *max_cutoff) {
    *max_cutoff = rho_f_.GetCutoff();
  }

  if (esw_f_.GetCutoff() > *max_cutoff) {
    *max_cutoff = esw_f_.GetCutoff();
  }

  // Everything is good
  return false;
}

}

#undef MEAM_SW_SPLINE_LOG_ERROR
#endif  // MEAM_SW_SPLINE_HPP