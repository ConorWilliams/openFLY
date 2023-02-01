//
// meam_spline.hpp
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
//        `lammps/src/USER-MISC/pair_meam_spline.h`
//        `lammps/src/USER-MISC/pair_meam_spline.cpp`
//        and it is rewritten and updated for the KIM-API by
//        Yaser Afshar
//

#ifndef MEAM_SPLINE_HPP
#define MEAM_SPLINE_HPP

#pragma GCC system_header

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "special.hpp"
#include "spline.hpp"
#include "utility.hpp"

namespace openKIM {

#ifdef MEAM_SPLINE_LOG_ERROR
#undef MEAM_SPLINE_LOG_ERROR
#endif

/*!
 * \brief Helper macro for printing error message
 * formats messages, filename, line number and function
 * name into an std::ostringstream object
 *
 */
#define MEAM_SPLINE_LOG_ERROR(msg)                                            \
  {                                                                           \
    std::ostringstream ss;                                                    \
    ss << "\nError :" << __FILE__ << ":" << __LINE__ << ":@(" << __FUNCTION__ \
       << ")\n"                                                               \
       << msg << "\n\n";                                                      \
    std::cerr << ss.str();                                                    \
  }

/*! \class MEAMSpline
 * \brief Modified Embedded Atom Method
 *
 * The MEAMSpline class is a base class to cmpute interactions for metals using
 * a variant of modified embedded-atom method (MEAM) potentials. It acounts for
 * a single species ('old-style') MEAM, as well as a new style multicomponent
 * modified embedded-atom method (MEAM) potential.
 *
 */
class MEAMSpline {
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
   * \param element_name Names of unique elements
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  inline int ProcessPotentialFile(std::FILE *const potential_file_pointer,
                                  int const max_line_size,
                                  std::vector<std::string> const &element_name);

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
  std::vector<double> zero_atom_energies_;

  /*! Used for temporary storage of U'(rho) values */
  std::vector<double> e_u_prime_;

  /*! Phi_i(r_ij)*/
  std::vector<Spline> e_phi_;
  /*! U_i(rho) */
  std::vector<Spline> e_u_;

  /*! Rho_ij(r_ij) */
  std::vector<Spline> rho_r_;
  /*! f_i(r_ij) */
  std::vector<Spline> rho_f_;
  /*! g_ij(cos_theta) */
  std::vector<Spline> rho_g_;
};

inline void MEAMSpline::ResizeLocalArray(std::size_t const nall) {
  // grow local arrays if necessary
  if (nall > e_u_prime_.size()) {
    // 2**14
    constexpr std::size_t kDelta = 16384;

    // Maximum number of atoms
    std::size_t const nmax = nall / kDelta * kDelta + kDelta;

    e_u_prime_.resize(nmax);
  }
}

inline int MEAMSpline::ProcessPotentialFile(
    std::FILE *const potential_file_pointer, int const max_line_size,
    std::vector<std::string> const &element_name) {
  std::rewind(potential_file_pointer);

  std::unique_ptr<char> next_line_(new char[max_line_size]);
  char *next_line = next_line_.get();

  // Create a utility object
  Utility ut;

  // Skip the first line. The first line is always a comment
  if (ut.GetFirstLine(potential_file_pointer, next_line, max_line_size)) {
    std::string msg = "End of file while reading the first line ";
    msg += "from the `meam/spline` potential file.\n";
    MEAM_SPLINE_LOG_ERROR(msg);
    return true;
  }

  int const number_of_elements = static_cast<int>(element_name.size());

  std::vector<int> species_codes(number_of_elements, 0);

  std::map<std::string, int> species_map;
  for (int ielem = 0; ielem < number_of_elements; ++ielem) {
    species_map[element_name[ielem]] = ielem;
  }

  bool is_new_format;

  // Second line holds potential type ("meam/spline")
  // in the new potential format.
  if (ut.GetNextLine(potential_file_pointer, next_line, max_line_size)) {
    std::string msg = "End of file while reading a line ";
    msg += "from the `meam/spline` potential file.\n";
    MEAM_SPLINE_LOG_ERROR(msg);
    return true;
  }

  char *word = std::strtok(next_line, " \t\n\r\f");

  if (std::strcmp(word, "meam/spline") == 0) {
    is_new_format = true;

    // parse the rest of the line!
    word = std::strtok(nullptr, " \t\n\r\f");
    if (word == nullptr) {
      std::string msg = "Missing to include the number of atomic species";
      msg += "in the multi-element `meam/spline` potential file.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }

    if (number_of_elements != std::atoi(word)) {
      std::string msg = "Invalid number of atomic species provided ";
      msg += "in the `meam/spline` potential file.\n";
      msg += "The number of atomic species in `meam/spline` ";
      msg += "potential file =" + std::to_string(std::atoi(word));
      msg += " != " + std::to_string(number_of_elements);
      msg += " which is provided in the MEAM element file.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }

    for (int ielem = 0; ielem < number_of_elements; ++ielem) {
      word = std::strtok(nullptr, " \t\n\r\f");
      if (word == nullptr) {
        std::string msg = "Not enough atomic species ";
        msg += "in the `meam/spline` potential file.\n";
        MEAM_SPLINE_LOG_ERROR(msg);
        return true;
      }

      std::string species_name(word);

      // Check for the species
      auto iter = species_map.find(species_name);
      if (iter == species_map.end()) {
        std::string msg = "Incorrect format in the `meam/spline` ";
        msg += "potential file.\nThe Species '" + species_name;
        msg += "' is not defined in the MEAM element file.\n";
        MEAM_SPLINE_LOG_ERROR(msg);
        return true;
      }

      species_codes[ielem] = species_map[species_name];
    }  // Loop over elements names
  } else {
    is_new_format = false;

    if (number_of_elements > 1) {
      std::string msg = "Incorrect format in the `meam/spline` ";
      msg += "potential file.\nThe `meam/spline` potential file ";
      msg += "is in the old format and the old format only handles ";
      msg += "one species, but there are more than one species ";
      msg += "are defined in the MEAM element file.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }

    std::rewind(potential_file_pointer);

    // Skip the first line. The first line is always a comment
    ut.GetFirstLine(potential_file_pointer, next_line, max_line_size);
  }

  int const nmultichoose2 = number_of_elements * (number_of_elements + 1) / 2;

  // Change the functional form
  // f_ij->f_i
  // g_i(cos\theta_ijk)->g_jk(cos\theta_ijk)
  zero_atom_energies_.resize(number_of_elements);

  e_phi_.resize(nmultichoose2);
  e_u_.resize(number_of_elements);

  rho_r_.resize(number_of_elements);
  rho_f_.resize(number_of_elements);
  rho_g_.resize(nmultichoose2);

  {
    int const species_i = species_codes[0];
    int const species_j = species_codes[0];

    int const ind = species_j + species_i * number_of_elements -
                    species_i * (species_i + 1) / 2;

    // Read and parse the input data
    if (e_phi_[ind].ProcessParameterFile(potential_file_pointer, max_line_size,
                                         is_new_format)) {
      std::string msg = "Unable to read and parse the `meam/spline` ";
      msg += "potential file for e_phi_.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }

    use_regular_grid_ = e_phi_[ind].AreKnotsOnRegularGrid();
  }

  for (int i = 0; i < number_of_elements; ++i) {
    int const species_i = species_codes[i];
    for (int j = i; j < number_of_elements; ++j) {
      if (i == 0 && j == 0) {
        continue;
      }
      int const species_j = species_codes[j];
      int const ind = species_j + species_i * number_of_elements -
                      species_i * (species_i + 1) / 2;
      if (e_phi_[ind].ProcessParameterFile(potential_file_pointer,
                                           max_line_size, is_new_format)) {
        std::string msg = "Unable to read and parse the `meam/spline` ";
        msg += "potential file for e_phi_.\n";
        MEAM_SPLINE_LOG_ERROR(msg);
        return true;
      }
      if (use_regular_grid_ != e_phi_[ind].AreKnotsOnRegularGrid()) {
        std::string msg = "This driver does not support both uniform & ";
        msg += "non-uniform cubic splines at the same time.\n";
        MEAM_SPLINE_LOG_ERROR(msg);
        return true;
      }
    }
  }

  for (int i = 0; i < number_of_elements; ++i) {
    int const ind = species_codes[i];
    if (rho_r_[ind].ProcessParameterFile(potential_file_pointer, max_line_size,
                                         is_new_format)) {
      std::string msg = "Unable to read and parse the `meam/spline` ";
      msg += "potential file for rho_r_.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
    if (use_regular_grid_ != rho_r_[ind].AreKnotsOnRegularGrid()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  for (int i = 0; i < number_of_elements; ++i) {
    int const ind = species_codes[i];
    if (e_u_[ind].ProcessParameterFile(potential_file_pointer, max_line_size,
                                       is_new_format)) {
      std::string msg = "Unable to read and parse the `meam/spline` ";
      msg += "potential file for e_u_.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
    if (use_regular_grid_ != e_u_[ind].AreKnotsOnRegularGrid()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  for (int i = 0; i < number_of_elements; ++i) {
    int const ind = species_codes[i];
    if (rho_f_[ind].ProcessParameterFile(potential_file_pointer, max_line_size,
                                         is_new_format)) {
      std::string msg = "Unable to read and parse the `meam/spline` ";
      msg += "potential file for rho_f_.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
    if (use_regular_grid_ != rho_f_[ind].AreKnotsOnRegularGrid()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  for (int i = 0; i < number_of_elements; ++i) {
    int const species_i = species_codes[i];
    for (int j = i; j < number_of_elements; ++j) {
      int const species_j = species_codes[j];
      int const ind = species_j + species_i * number_of_elements -
                      species_i * (species_i + 1) / 2;
      if (rho_g_[ind].ProcessParameterFile(potential_file_pointer,
                                           max_line_size, is_new_format)) {
        std::string msg = "Unable to read and parse the `meam/spline` ";
        msg += "potential file for rho_g_.\n";
        MEAM_SPLINE_LOG_ERROR(msg);
        return true;
      }
      if (use_regular_grid_ != rho_g_[ind].AreKnotsOnRegularGrid()) {
        std::string msg = "This driver does not support both uniform & ";
        msg += "non-uniform cubic splines at the same time.\n";
        MEAM_SPLINE_LOG_ERROR(msg);
        return true;
      }
    }
  }

  // Everything is good
  return false;
}

inline void MEAMSpline::ConvertUnit(double const convert_length_factor,
                                    double const convert_energy_factor) {
  for (std::size_t i = 0; i < e_phi_.size(); ++i) {
    e_phi_[i].ConvertUnit(convert_length_factor, convert_energy_factor);
  }

  for (std::size_t i = 0; i < e_u_.size(); ++i) {
    e_u_[i].ConvertUnit(1.0, convert_energy_factor);
  }

  for (std::size_t i = 0; i < rho_r_.size(); ++i) {
    rho_r_[i].ConvertUnit(convert_length_factor, 1.0);
  }

  for (std::size_t i = 0; i < rho_f_.size(); ++i) {
    rho_f_[i].ConvertUnit(convert_length_factor, 1.0);
  }
}

inline int MEAMSpline::CompleteSetup(double *max_cutoff) {
  for (std::size_t i = 0; i < e_phi_.size(); ++i) {
    if (e_phi_[i].UpdateSpline()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  for (std::size_t i = 0; i < e_u_.size(); ++i) {
    if (e_u_[i].UpdateSpline()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  for (std::size_t i = 0; i < rho_r_.size(); ++i) {
    if (rho_r_[i].UpdateSpline()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  for (std::size_t i = 0; i < rho_f_.size(); ++i) {
    if (rho_f_[i].UpdateSpline()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  for (std::size_t i = 0; i < rho_g_.size(); ++i) {
    if (rho_g_[i].UpdateSpline()) {
      std::string msg = "This driver does not support both uniform & ";
      msg += "non-uniform cubic splines at the same time.\n";
      MEAM_SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  // Calculate 'zero-point energy' of a single atom in vacuum.
  if (use_regular_grid_) {
    for (std::size_t i = 0; i < zero_atom_energies_.size(); ++i) {
      zero_atom_energies_[i] = e_u_[i].Eval<true>(0.0);
    }

  } else {
    for (std::size_t i = 0; i < zero_atom_energies_.size(); ++i) {
      zero_atom_energies_[i] = e_u_[i].Eval<false>(0.0);
    }
  }

  // Determine the maximum cutoff radius of all relevant spline functions.

  *max_cutoff = 0.0;

  for (std::size_t i = 0; i < e_phi_.size(); ++i) {
    if (e_phi_[i].GetCutoff() > *max_cutoff) {
      *max_cutoff = e_phi_[i].GetCutoff();
    }
  }

  for (std::size_t i = 0; i < rho_r_.size(); ++i) {
    if (rho_r_[i].GetCutoff() > *max_cutoff) {
      *max_cutoff = rho_r_[i].GetCutoff();
    }
  }

  for (std::size_t i = 0; i < rho_f_.size(); ++i) {
    if (rho_f_[i].GetCutoff() > *max_cutoff) {
      *max_cutoff = rho_f_[i].GetCutoff();
    }
  }

  // Everything is good
  return false;
}

}

#undef MEAM_SPLINE_LOG_ERROR
#endif  // MEAM_SPLINE_HPP