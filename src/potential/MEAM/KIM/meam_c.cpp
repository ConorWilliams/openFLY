//
// meam_c.cpp
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
//        `lammps/src/USER-MEAMC/meam_dens_final.cpp`
//        `lammps/src/USER-MEAMC/meam_dens_init.cpp`
//        `lammps/src/USER-MEAMC/meam_force.cpp`
//        `lammps/src/USER-MEAMC/meam_funcs.cpp`
//        `lammps/src/USER-MEAMC/meam_impl.cpp`
//        `lammps/src/USER-MEAMC/meam_setup_done.cpp`
//        `lammps/src/USER-MEAMC/meam_setup_global.cpp`
//        `lammps/src/USER-MEAMC/meam_setup_param.cpp`
//        `lammps/src/USER-MEAMC/meam.h`
//        `lammps/src/USER-MEAMC/pair_meamc.cpp`
//        `lammps/src/USER-MEAMC/pair_meamc.h`
//        and it is rewritten and updated for the KIM-API by
//        Yaser Afshar
//

#include "libfly/potential/MEAM/KIM/meam_c.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <memory>

#include "libfly/potential/MEAM/KIM/special.hpp"
#include "libfly/potential/MEAM/KIM/utility.hpp"

namespace {
  static constexpr std::size_t kDelta = 16384;                  // 2**14
  static constexpr std::size_t kDeltaN = 1024;                  // 2**10
  static constexpr double kThird = 0.33333333333333333333;      // 1/3
  static constexpr double kTwoThird = 0.66666666666666666666;   // 2/3
  static constexpr double kSqrtTwo = 1.41421356237309504880;    // sqrt(2)
  static constexpr double kSqrtThree = 1.73205080756887729352;  // sqrt(3)

  /*!
   * \brief Array of factors to apply for Voight notation.
   *        multiplicity of Voigt index (i.e. [1] -> xy+yx = 2)
   */
  static constexpr int const kVoight2dFactor[6] = {1, 2, 2, 1, 2, 1};

  /*!
   * \brief Array of factors to apply for Voight notation.
   *        (multiplicity of Voigt index)
   */
  static constexpr int const kVoight3dFactor[10] = {1, 3, 3, 3, 6, 3, 1, 3, 3, 1};

  static constexpr char const *const keywords[]
      = {"Ec",    "alpha",    "rho0",          "delta",      "lattce", "attrac",      "repuls",
         "nn2",   "Cmin",     "Cmax",          "rc",         "delr",   "augt1",       "gsmooth_factor",
         "re",    "ialloy",   "mixture_ref_t", "erose_form", "zbl",    "emb_lin_neg", "bkgd_dyn",
         "theta", "fcut_form"};

  static constexpr int nkeywords = 23;
}  // namespace

namespace openKIM {

  void MEAMC::GetShapeFactors(Lattice const &lat,
                              double const stheta,
                              double const ctheta,
                              VectorOfSizeDIM &shape_factors) {
    switch (lat) {
      case Lattice::FCC:
      case Lattice::BCC:
      case Lattice::B1:
      case Lattice::B2: {
        shape_factors[0] = 0.0;
        shape_factors[1] = 0.0;
        shape_factors[2] = 0.0;
        return;
      }
      case Lattice::HCP: {
        shape_factors[0] = 0.0;
        shape_factors[1] = 0.0;
        shape_factors[2] = kThird;
        return;
      }
      // CH4 actually needs shape factor for diamond for C, dimer for H
      case Lattice::CH4:
      case Lattice::DIA:
      case Lattice::DIA3: {
        shape_factors[0] = 0.0;
        shape_factors[1] = 0.0;
        shape_factors[2] = 32.0 / 9.0;
        return;
      }
      case Lattice::DIM: {
        shape_factors[0] = 1.0;
        shape_factors[1] = kTwoThird;
        shape_factors[2] = 0.40;
        return;
      }
      // linear, theta being 180
      case Lattice::LIN: {
        shape_factors[0] = 0.0;
        shape_factors[1] = 8.0 / 3.0;
        shape_factors[2] = 0.0;
        return;
      }
      // zig-zag
      case Lattice::ZIG:
      // trimer e.g. H2O
      case Lattice::TRI: {
        shape_factors[0] = 4.0 * special::Square(ctheta);
        shape_factors[1] = 4.0 * (special::PowInt(ctheta, 4) + special::PowInt(stheta, 4) - kThird);
        shape_factors[2]
            = 4.0 * (special::Square(ctheta) * (3 * special::PowInt(stheta, 4) + special::PowInt(ctheta, 4)));
        // legend in dyn, 0.6 is the default value.
        shape_factors[2] -= 0.6 * shape_factors[0];
        return;
      }
      // C11, L12,
      default: {
        shape_factors[0] = 0.0;
      }
    }
  }

  std::string MEAMC::LatticeToString(Lattice const &lat) {
    switch (lat) {
      case Lattice::FCC: {
        return "fcc";
      }
      case Lattice::BCC: {
        return "bcc";
      }
      case Lattice::HCP: {
        return "hcp";
      }
      case Lattice::DIM: {
        return "dim";
      }
      case Lattice::DIA: {
        return "dia";
      }
      case Lattice::DIA3: {
        return "dia3";
      }
      case Lattice::LIN: {
        return "lin";
      }
      case Lattice::ZIG: {
        return "zig";
      }
      case Lattice::TRI: {
        return "tri";
      }
      case Lattice::B1: {
        return "b1";
      }
      case Lattice::C11: {
        return "c11";
      }
      case Lattice::L12: {
        return "l12";
      }
      case Lattice::B2: {
        return "b2";
      }
      case Lattice::CH4: {
        return "ch4";
      }
      default: {
        // unknown lattic flag
        return "";
      }
    }
  }

  int MEAMC::ProcessLibraryFile(std::FILE *const library_file_pointer,
                                int const max_line_size,
                                std::vector<std::string> const &element_name) {
    std::rewind(library_file_pointer);

    std::unique_ptr<char> next_line_(new char[max_line_size]);
    char *next_line = next_line_.get();

    // Create a utility object
    Utility ut;

    // # of unique element types
    number_of_element_types_ = static_cast<int>(element_name.size());

    // Grow element arrays and initialize them to some default
    ResizeElementArrays();

    // Temporray arrays
    std::vector<bool> found(number_of_element_types_, false);
    std::vector<double> mass(number_of_element_types_);
    std::vector<double> z(number_of_element_types_);
    std::vector<double> element_t0(number_of_element_types_);
    std::vector<double> element_lattice_constant(number_of_element_types_);

    // Read each set of parameters from the MEAM file.
    // One set of params can span multiple lines.
    // Store the parameters, if element name is in the element list.
    // If the element name appears multiple times, only store the 1st entry.
    constexpr int params_per_line = 19;
    char *words[params_per_line + 1];

    int nset = 0;
    while (true) {
      // read a line, strip comment, skip line if blank
      if (ut.GetNextLine(library_file_pointer, next_line, max_line_size)) {
        // We reached the end of the file and it should exit the loop
        break;
      }

      int nwords = ut.GetNumberOfWordsInLine(next_line);

      // concatenate additional lines until have params_per_line words
      while (nwords < params_per_line) {
        int const n = std::strlen(next_line);

        if (ut.GetNextLine(library_file_pointer, &next_line[n], max_line_size - n)) {
          std::string msg = "End of file while reading a line from the ";
          msg += "`meam/c` library file.\n";
          HELPER_LOG_ERROR(msg);
          return true;
        }

        nwords = ut.GetNumberOfWordsInLine(next_line);
      }

      if (nwords != params_per_line) {
        HELPER_LOG_ERROR("Incorrect format in the `meam/c` library file.\n");
        return true;
      }

      // Because the library file is used by Fortran MD codes, these
      // strings may be enclosed in single quotes. Thus strip single
      // quotes from words
      nwords = 0;
      words[nwords++] = std::strtok(next_line, "' \t\n\r\f");

      std::string const species_name(words[0]);

      int index;
      for (index = 0; index < number_of_element_types_; ++index) {
        if (species_name == element_name[index]) {
          break;
        }
      }

      // Skip if the species name isn't in the element list
      if (index == number_of_element_types_) {
        continue;
      }

      // Element name is already in the element list
      if (found[index]) {
        // skip if the element already appeared (technically error in
        // the library file, but always ignored)
        continue;
      }

      while ((words[nwords] = std::strtok(nullptr, "' \t\n\r\f"))) {
        ++nwords;
      }

      found[index] = true;

      // map element_lattice string to an integer
      if (StringToLattice(words[1], true, element_lattice_[index][index])) {
        std::string msg = "Unrecognized lattice type in the `meam/c` ";
        msg += "library file: '" + std::string(words[1]) + "'.\n";
        HELPER_LOG_ERROR(msg);
        return true;
      }

      // Store parameters
      z[index] = std::atof(words[2]);
      // if z given is mismatched -> fatal error
      if (special::IsNotZero(z[index]
                             - NumNearestNeighborsInReferenceStructure(element_lattice_(index, index)))) {
        std::string msg = "Mismatched parameter in the `meam/c` library file: ";
        msg += "'z != number of nearest neighbors in the reference structure ";
        msg += "lattice'.\n";
        HELPER_LOG_ERROR(msg);
        return true;
      }

      element_atomic_number_[index] = std::atof(words[3]);
      mass[index] = std::atof(words[4]);
      element_alpha_(index, index) = std::atof(words[5]);
      element_beta0_[index] = std::atof(words[6]);
      element_beta1_[index] = std::atof(words[7]);
      element_beta2_[index] = std::atof(words[8]);
      element_beta3_[index] = std::atof(words[9]);
      element_lattice_constant[index] = std::atof(words[10]);
      element_Ec_(index, index) = std::atof(words[11]);
      element_A_[index] = std::atof(words[12]);

      element_t0[index] = std::atof(words[13]);
      if (special::IsNotOne(element_t0[index])) {
        std::string msg = "Unsupported parameter in the `meam/c` ";
        msg += "library file: 't0 != 1'.\n";
        HELPER_LOG_ERROR(msg);
        return true;
      }

      element_t1_[index] = std::atof(words[14]);
      element_t2_[index] = std::atof(words[15]);
      element_t3_[index] = std::atof(words[16]);
      element_rho0_[index] = std::atof(words[17]);

      element_ibar_[index] = std::atoi(words[18]);
      if (!(element_ibar_[index] == 0 || element_ibar_[index] == 1 || element_ibar_[index] == 3
            || element_ibar_[index] == 4 || element_ibar_[index] == -5)) {
        std::string msg = "Unsupported form of the function `G(Gamma)` ";
        msg += "in the `meam/c` library file is requested.\n";
        HELPER_LOG_ERROR(msg);
        return true;
      }

      ++nset;
    }  // while (true)

    // error if didn't find all elements in file
    if (nset != number_of_element_types_) {
      std::string msg = "Did not find all the requested elements in ";
      msg += "the `meam/c` library file.\nElements: '";
      msg += element_name[0];
      for (int i = 1; i < number_of_element_types_; ++i) {
        if (!found[i]) {
          msg += ", " + element_name[i];
        }
      }
      msg += " ' are missing in the `meam/c` library file.\n";
      HELPER_LOG_ERROR(msg);
      return true;
    }

    for (int i = 0; i < number_of_element_types_; ++i) {
      switch (element_lattice_(i, i)) {
        case Lattice::FCC: {
          element_re_(i, i) = element_lattice_constant[i] / kSqrtTwo;
          break;
        }
        case Lattice::BCC: {
          element_re_(i, i) = element_lattice_constant[i] * kSqrtThree * 0.5;
          break;
        }
        case Lattice::HCP:
        case Lattice::DIM:
        case Lattice::CH4:
        case Lattice::LIN:
        case Lattice::ZIG:
        case Lattice::TRI: {
          element_re_(i, i) = element_lattice_constant[i];
          break;
        }
        case Lattice::DIA:
        case Lattice::DIA3: {
          element_re_(i, i) = element_lattice_constant[i] * kSqrtThree * 0.25;
          break;
        }
        case Lattice::B1:
        case Lattice::B2:
        case Lattice::C11:
        case Lattice::L12:
          // do nothing
          break;
        default:;
          //  error
      }
    }

    // Everything is good
    return false;
  }

  int MEAMC::ProcessParameterFile(std::FILE *const parameter_file_pointer, int const max_line_size) {
    std::unique_ptr<char> next_line_(new char[max_line_size]);
    char *next_line = next_line_.get();

    // Create a utility object
    Utility ut;

    // read settings and pass them one at a time to MEAM package
    // match strings to list of corresponding ints

    constexpr int kMaxParams = 6;
    char *params[kMaxParams];

    int nset = 0;
    while (true) {
      // strip comment, skip line if blank
      if (ut.GetNextLine(parameter_file_pointer, next_line, max_line_size)) {
        break;
      }

      // params = ptrs to all fields in line

      int nparams = 0;
      params[nparams++] = std::strtok(next_line, "=(), '\t\n\r\f");

      while (nparams < kMaxParams && (params[nparams] = std::strtok(nullptr, "=(), '\t\n\r\f"))) {
        ++nparams;
      }

      int which;
      for (which = 0; which < nkeywords; ++which) {
        if (std::strcmp(params[0], keywords[which]) == 0) {
          break;
        }
      }

      if (which == nkeywords) {
        std::string msg = "The Keyword '" + std::string(params[0]);
        msg += "' in the `meam/c` parameter file is not recognized.\n";
        HELPER_LOG_ERROR(msg);
        return true;
      }

      int const nindex = nparams - 2;
      int index[3];
      for (int i = 0; i < nindex; ++i) {
        index[i] = std::atoi(params[i + 1]) - 1;
      }

      // map lattce_meam value to an integer

      double value;

      if (which == 4) {
        Lattice latt;
        if (StringToLattice(params[nparams - 1], false, latt)) {
          std::string msg = "Unrecognized lattice type in the `meam/c` ";
          msg += "parameter file: " + std::string(params[nparams - 1]);
          HELPER_LOG_ERROR(msg);
          return true;
        }
        value = static_cast<double>(latt);
      } else {
        value = std::atof(params[nparams - 1]);
      }

      // pass single setting to MEAM package
      int error_flag = 0;
      SetParameter(which, value, nindex, index, &error_flag);

      if (error_flag) {
        if ((error_flag < 0) || (error_flag > 5)) {
          error_flag = 0;
        }

        std::vector<std::string> descr = {"has an unknown error",                              // 0
                                          "is out of range (please report a bug)",             // 1
                                          "expected more indices",                             // 2
                                          "has out of range element index",                    // 3
                                          "has an unknown form of the MEAM cutoff function",   // 4
                                          "has an unknown form of the rose energy function"};  // 5

        std::string msg = "Error in the `meam/c` parameter file: keyword '";
        msg += std::string(params[0]) + "' " + descr[error_flag] + "\n";
        HELPER_LOG_ERROR(msg);
        return true;
      }

      ++nset;
    }

    // error if the `meam/c` parameter file is empty
    if (nset == 0) {
      HELPER_LOG_ERROR("The `meam/c` parameter file is empty.\n");
      return true;
    }

    // Everything is good
    return false;
  }

  void MEAMC::ConvertUnit(double const convert_length_factor, double const convert_energy_factor) {
    // cutoff_radius_, delr, re are in distance units
    // (Angstroms in the case of metal units).
    if (special::IsNotOne(convert_length_factor)) {
      cutoff_radius_ *= convert_length_factor;

      delr_ *= convert_length_factor;

      for (int i = 0; i < number_of_element_types_; ++i) {
        for (int j = 0; j < number_of_element_types_; ++j) {
          element_re_(i, j) *= convert_length_factor;
        }
      }
    }

    if (special::IsNotOne(convert_energy_factor)) {
      //  Ec and element_delta_ are in energy units
      //  (eV in the case of metal units).
      for (int i = 0; i < number_of_element_types_; ++i) {
        for (int j = 0; j < number_of_element_types_; ++j) {
          element_Ec_(i, j) *= convert_energy_factor;
        }
      }

      for (int i = 0; i < number_of_element_types_; ++i) {
        for (int j = 0; j < number_of_element_types_; ++j) {
          element_delta_(i, j) *= convert_energy_factor;
        }
      }
    }

    zbl_->ConvertUnit(convert_length_factor, convert_energy_factor);
  }

  void MEAMC::CompleteSetup(double *max_cutoff) {
    // Pass cutoff back to calling program
    *max_cutoff = cutoff_radius_;

    // squared force cutoff
    cutoff_radius_squared_ = cutoff_radius_ * cutoff_radius_;

    // Augment t1 term
    for (int i = 0; i < number_of_element_types_; ++i) {
      element_t1_[i] += augment_t1_ * 3.0 / 5.0 * element_t3_[i];
    }

    // we don't use theta, instead stheta and ctheta
    // Loop over pairs
    for (int i = 0; i < number_of_element_types_; ++i) {
      for (int j = i; j < number_of_element_types_; ++j) {
        if (special::IsZero(element_theta_(i, j) - 180.0)) {
          // stheta = sin(theta/2*pi/180) where theta is 180, so 1.0
          element_stheta_(i, j) = 1.0;
          // ctheta = cos(theta/2*pi/180) where theta is 180, so 0
          element_ctheta_(i, j) = 0.0;
        } else {
          double const value = element_theta_(i, j) * 0.5 * special::MY_PI / 180.0;
          element_stheta_(i, j) = std::sin(value);
          element_ctheta_(i, j) = std::cos(value);
        }
      }
    }

    // Compute off-diagonal alloy parameters
    FillOffDiagonalAlloyParameters();

    // indices number of pairs
    for (int m = 0, nv2 = 0; m < number_of_element_types_; ++m) {
      for (int n = m; n < number_of_element_types_; ++n, ++nv2) {
        element_pair_index_(m, n) = element_pair_index_(n, m) = nv2;
      }
    }

    // Compute the composition dependent electron density scaling Eq.4.5
    ComputeCompositionDependentDensityScaling();

    // Pair function discretization parameters
    dr_ = 1.1 * cutoff_radius_ / static_cast<double>(nr_);

    // Compute pair potentials and setup arrays for interpolation
    ComputePairPotential();
  }

  void MEAMC::ResizeElementArrays() {
    std::size_t const n = static_cast<std::size_t>(number_of_element_types_);

    element_ibar_.resize(n);
    element_atomic_number_.resize(n);
    element_beta0_.resize(n);
    element_beta1_.resize(n);
    element_beta2_.resize(n);
    element_beta3_.resize(n);
    element_A_.resize(n);
    element_t1_.resize(n);
    element_t2_.resize(n);
    element_t3_.resize(n);
    element_rho0_.resize(n);
    element_ref_rho_.resize(n, 0.0);
    element_lattice_.resize(n, n, Lattice::FCC);
    element_nn2_.resize(n, n, 0);
    element_zbl_.resize(n, n, 1);
    element_pair_index_.resize(n, n, 0);
    element_alpha_.resize(n, n, 0.0);
    element_re_.resize(n, n, 0.0);
    element_Ec_.resize(n, n, 0.0);
    element_delta_.resize(n, n, 0.0);
    element_attrac_.resize(n, n, 0.0);
    element_repuls_.resize(n, n, 0.0);
    // for trimer, zigzag, line refernece structure, sungkwang
    element_theta_.resize(n, n, 180.);
    element_stheta_.resize(n, n);
    element_ctheta_.resize(n, n);
    element_ebound_.resize(n, n, (2.8 * 2.8) / (4.0 * (2.8 - 1.0)));
    element_Cmin_.resize(n, n, n, 2.0);
    element_Cmax_.resize(n, n, n, 2.8);

    // Construct the ZBL object
    zbl_.reset(new ZBL(n));
  }

  void MEAMC::ResizeDenistyArrays(std::size_t const nall) {
    // grow local arrays if necessary

    if (nall > rho_.size()) {
      // Maximum number of atoms
      std::size_t const nmax = nall / kDelta * kDelta + kDelta;

      rho_.resize(nmax);
      frhop_.resize(nmax);

      rho0_.resize(nmax);
      rho1_.resize(nmax);
      rho2_.resize(nmax);
      rho3_.resize(nmax);

      gamma_.resize(nmax);
      dgamma1_.resize(nmax);
      dgamma2_.resize(nmax);
      dgamma3_.resize(nmax);

      arho1_.resize(nmax, 3);
      arho2_.resize(nmax, 6);
      arho2b_.resize(nmax);
      arho3_.resize(nmax, 10);
      arho3b_.resize(nmax, 3);

      t_ave_.resize(nmax, 3);
      tsq_ave_.resize(nmax, 3);
    }

    // zero out local arrays
    std::fill_n(rho0_.data(), nall, 0.0);

    std::fill_n(arho1_.data(), nall * 3, 0.0);
    std::fill_n(arho2_.data(), nall * 6, 0.0);
    std::fill_n(arho2b_.data(), nall, 0.0);
    std::fill_n(arho3_.data(), nall * 10, 0.0);
    std::fill_n(arho3b_.data(), nall * 3, 0.0);

    std::fill_n(t_ave_.data(), nall * 3, 0.0);
    std::fill_n(tsq_ave_.data(), nall * 3, 0.0);
  }

  void MEAMC::ResizeScreeningArrays(std::size_t const n_neigh) {
    // grow local arrays if necessary
    if (n_neigh > scrfcn_.size()) {
      // Maximum number of atoms neighbors
      std::size_t const nmax = n_neigh / kDeltaN * kDeltaN + kDeltaN;

      scrfcn_.resize(nmax);
      dscrfcn_.resize(nmax);
    }
  }

  void MEAMC::InitializeDensityCalculation(int const i,
                                           int const number_of_neighbors,
                                           int const *const neighbors_of_particle,
                                           int &offset,
                                           const VectorOfSizeDIM *const coordinates,
                                           int const *const particle_species_codes,
                                           int const *const particle_contributing) {
    // Compute screening function and derivative
    ComputeScreeningAndDerivative(i,
                                  number_of_neighbors,
                                  neighbors_of_particle,
                                  offset,
                                  coordinates,
                                  particle_species_codes,
                                  particle_contributing);

    // Calculate intermediate density terms
    ComputeIntermediateDensityTerms(i,
                                    number_of_neighbors,
                                    neighbors_of_particle,
                                    offset,
                                    coordinates,
                                    particle_species_codes,
                                    particle_contributing);
  }

  void MEAMC::FinalizeDensityCalculation(int const i, int const species_i, double &embedding) {
    double ts;

    double const *const arho1_1d_i = &arho1_(i, 0);

    // Eq.4.7b
    ts = arho1_1d_i[0] * arho1_1d_i[0];
    ts += arho1_1d_i[1] * arho1_1d_i[1];
    ts += arho1_1d_i[2] * arho1_1d_i[2];
    rho1_[i] = ts;

    double const *const arho2_1d_i = &arho2_(i, 0);

    // Eq.4.7c
    ts = -kThird * arho2b_[i] * arho2b_[i];
    ts += kVoight2dFactor[0] * arho2_1d_i[0] * arho2_1d_i[0];
    ts += kVoight2dFactor[1] * arho2_1d_i[1] * arho2_1d_i[1];
    ts += kVoight2dFactor[2] * arho2_1d_i[2] * arho2_1d_i[2];
    ts += kVoight2dFactor[3] * arho2_1d_i[3] * arho2_1d_i[3];
    ts += kVoight2dFactor[4] * arho2_1d_i[4] * arho2_1d_i[4];
    ts += kVoight2dFactor[5] * arho2_1d_i[5] * arho2_1d_i[5];
    rho2_[i] = ts;

    double const *const arho3b_1d_i = &arho3b_(i, 0);

    // Eq.4.7d
    ts = arho3b_1d_i[0] * arho3b_1d_i[0];
    ts += arho3b_1d_i[1] * arho3b_1d_i[1];
    ts += arho3b_1d_i[2] * arho3b_1d_i[2];
    ts *= -3.0 / 5.0;

    double const *const arho3_1d_i = &arho3_(i, 0);

    ts += kVoight3dFactor[0] * arho3_1d_i[0] * arho3_1d_i[0];
    ts += kVoight3dFactor[1] * arho3_1d_i[1] * arho3_1d_i[1];
    ts += kVoight3dFactor[2] * arho3_1d_i[2] * arho3_1d_i[2];
    ts += kVoight3dFactor[3] * arho3_1d_i[3] * arho3_1d_i[3];
    ts += kVoight3dFactor[4] * arho3_1d_i[4] * arho3_1d_i[4];
    ts += kVoight3dFactor[5] * arho3_1d_i[5] * arho3_1d_i[5];
    ts += kVoight3dFactor[6] * arho3_1d_i[6] * arho3_1d_i[6];
    ts += kVoight3dFactor[7] * arho3_1d_i[7] * arho3_1d_i[7];
    ts += kVoight3dFactor[8] * arho3_1d_i[8] * arho3_1d_i[8];
    ts += kVoight3dFactor[9] * arho3_1d_i[9] * arho3_1d_i[9];
    rho3_[i] = ts;

    double *const t_ave_1d_i = &t_ave_(i, 0);

    if (rho0_[i] > 0.0) {
      if (ialloy_ == 1) {
        double const *const tsq_ave_1d_i = &tsq_ave_(i, 0);

        t_ave_1d_i[0] = special::FloatDivZero(t_ave_1d_i[0], tsq_ave_1d_i[0]);
        t_ave_1d_i[1] = special::FloatDivZero(t_ave_1d_i[1], tsq_ave_1d_i[1]);
        t_ave_1d_i[2] = special::FloatDivZero(t_ave_1d_i[2], tsq_ave_1d_i[2]);
      } else if (ialloy_ == 2) {
        t_ave_1d_i[0] = element_t1_[species_i];
        t_ave_1d_i[1] = element_t2_[species_i];
        t_ave_1d_i[2] = element_t3_[species_i];
      } else {
        // Eq.4.9
        t_ave_1d_i[0] /= rho0_[i];
        t_ave_1d_i[1] /= rho0_[i];
        t_ave_1d_i[2] /= rho0_[i];
      }
    }

    // Eq.4.4
    gamma_[i] = t_ave_1d_i[0] * rho1_[i] + t_ave_1d_i[1] * rho2_[i] + t_ave_1d_i[2] * rho3_[i];
    if (rho0_[i] > 0.0) {
      gamma_[i] /= (rho0_[i] * rho0_[i]);
    }

    int const ibar_i = element_ibar_[species_i];

    double g_gamma = GGamma(gamma_[i], ibar_i);

    double const z = NumNearestNeighborsInReferenceStructure(element_lattice_(species_i, species_i));

    VectorOfSizeDIM shape_factors;

    GetShapeFactors(element_lattice_(species_i, species_i),
                    element_stheta_(species_i, species_i),
                    element_ctheta_(species_i, species_i),
                    shape_factors);

    double gbar;
    double dgbar;

    if (ibar_i <= 0) {
      gbar = 1.0;
      dgbar = 0.0;
    } else {
      if (mixing_rule_compute_t_ == 1) {
        double const gam = (t_ave_1d_i[0] * shape_factors[0] + t_ave_1d_i[1] * shape_factors[1]
                            + t_ave_1d_i[2] * shape_factors[2])
                           / (z * z);

        gbar = GGamma(gam, ibar_i, dgbar);
      } else {
        double const gam
            = (element_t1_[species_i] * shape_factors[0] + element_t2_[species_i] * shape_factors[1]
               + element_t3_[species_i] * shape_factors[2])
              / (z * z);

        gbar = GGamma(gam, ibar_i);
      }
    }

    double rho_bkgd;

    if (mixing_rule_compute_t_ == 1) {
      // Eq.4.5
      rho_bkgd = element_rho0_[species_i] * z * gbar;
    } else {
      rho_bkgd
          = (dynamo_reference_density_ == 1) ? element_rho0_[species_i] * z : element_ref_rho_[species_i];
    }

    // in Eq.4.3 rho_^bar = rho0_/rho0i * G
    rho_[i] = rho0_[i] * g_gamma;
    // Eq.4.3
    double const rhob = rho_[i] / rho_bkgd;

    double const denom = 1.0 / rho_bkgd;

    double dg_gamma;
    g_gamma = GGamma(gamma_[i], ibar_i, dg_gamma);

    // In Eq.4.36a
    dgamma1_[i] = (g_gamma - 2 * dg_gamma * gamma_[i]) * denom;

    // In Eq.4.36a
    dgamma2_[i] = special::FloatDivZero(dg_gamma * denom, rho0_[i]);

    // dgamma3_ is nonzero only if we are using the "mixed" rule for
    // computing t in the reference system (which is not correct, but
    // included for backward compatibility
    dgamma3_[i] = (mixing_rule_compute_t_ != 1) ? 0.0 : rho0_[i] * g_gamma * dgbar / (gbar * z * z) * denom;

    // Eq.4.2
    embedding = Embedding(element_A_[species_i], element_Ec_(species_i, species_i), rhob, frhop_[i]);
  }

  void MEAMC::ComputeAtomicElectronDensities(int const species_i,
                                             int const species_j,
                                             double const rij,
                                             double &rhoa0_i,
                                             double &drhoa0_i,
                                             double &rhoa1_i,
                                             double &drhoa1_i,
                                             double &rhoa2_i,
                                             double &drhoa2_i,
                                             double &rhoa3_i,
                                             double &drhoa3_i,
                                             double &rhoa0_j,
                                             double &drhoa0_j,
                                             double &rhoa1_j,
                                             double &drhoa1_j,
                                             double &rhoa2_j,
                                             double &drhoa2_j,
                                             double &rhoa3_j,
                                             double &drhoa3_j) {
    // Compute pair densities and derivatives
    double const invrei = 1.0 / element_re_(species_i, species_i);

    // (rij/r0i - 1) in Eq.4.8
    double const ai = rij * invrei - 1.0;

    // Element-dependent density scaling in Eq.4.5 or Eq.4.8
    double const rhoi0 = element_rho0_[species_i];

    // The atomic electron densities. Eq.4.8
    rhoa0_i = rhoi0 * std::exp(-element_beta0_[species_i] * ai);
    drhoa0_i = -element_beta0_[species_i] * invrei * rhoa0_i;

    rhoa1_i = rhoi0 * std::exp(-element_beta1_[species_i] * ai);
    drhoa1_i = -element_beta1_[species_i] * invrei * rhoa1_i;

    rhoa2_i = rhoi0 * std::exp(-element_beta2_[species_i] * ai);
    drhoa2_i = -element_beta2_[species_i] * invrei * rhoa2_i;

    rhoa3_i = rhoi0 * std::exp(-element_beta3_[species_i] * ai);
    drhoa3_i = -element_beta3_[species_i] * invrei * rhoa3_i;

    if (species_i == species_j) {
      rhoa0_j = rhoa0_i;
      drhoa0_j = drhoa0_i;

      rhoa1_j = rhoa1_i;
      drhoa1_j = drhoa1_i;

      rhoa2_j = rhoa2_i;
      drhoa2_j = drhoa2_i;

      rhoa3_j = rhoa3_i;
      drhoa3_j = drhoa3_i;
    } else {
      double const invrej = 1.0 / element_re_(species_j, species_j);

      // (rij/r0i - 1) in Eq.4.8
      double const aj = rij * invrej - 1.0;

      double const rhoj0 = element_rho0_[species_j];

      rhoa0_j = rhoj0 * std::exp(-element_beta0_[species_j] * aj);
      drhoa0_j = -element_beta0_[species_j] * invrej * rhoa0_j;

      rhoa1_j = rhoj0 * std::exp(-element_beta1_[species_j] * aj);
      drhoa1_j = -element_beta1_[species_j] * invrej * rhoa1_j;

      rhoa2_j = rhoj0 * std::exp(-element_beta2_[species_j] * aj);
      drhoa2_j = -element_beta2_[species_j] * invrej * rhoa2_j;

      rhoa3_j = rhoj0 * std::exp(-element_beta3_[species_j] * aj);
      drhoa3_j = -element_beta3_[species_j] * invrej * rhoa3_j;
    }

    if (ialloy_ == 1) {
      double const t1m_i = element_t1_[species_i];
      double const t2m_i = element_t2_[species_i];
      double const t3m_i = element_t3_[species_i];

      rhoa1_i *= t1m_i;
      rhoa2_i *= t2m_i;
      rhoa3_i *= t3m_i;

      drhoa1_i *= t1m_i;
      drhoa2_i *= t2m_i;
      drhoa3_i *= t3m_i;

      double const t1m_j = element_t1_[species_j];
      double const t2m_j = element_t2_[species_j];
      double const t3m_j = element_t3_[species_j];

      rhoa1_j *= t1m_j;
      rhoa2_j *= t2m_j;
      rhoa3_j *= t3m_j;

      drhoa1_j *= t1m_j;
      drhoa2_j *= t2m_j;
      drhoa3_j *= t3m_j;
    }
  }

  bool MEAMC::StringToLattice(std::string const &str, bool const single, Lattice &lat) {
    if (str == "fcc") {
      lat = Lattice::FCC;
    } else if (str == "bcc") {
      lat = Lattice::BCC;
    } else if (str == "hcp") {
      lat = Lattice::HCP;
    } else if (str == "dim") {
      lat = Lattice::DIM;
    } else if (str == "dia") {
      lat = Lattice::DIA;
    } else if (str == "dia3") {
      lat = Lattice::DIA3;
    } else if (str == "lin") {
      lat = Lattice::LIN;
    } else if (str == "zig") {
      lat = Lattice::ZIG;
    } else if (str == "tri") {
      lat = Lattice::TRI;
    } else {
      if (single) {
        return true;
      }
      if (str == "b1") {
        lat = Lattice::B1;
      } else if (str == "c11") {
        lat = Lattice::C11;
      } else if (str == "l12") {
        lat = Lattice::L12;
      } else if (str == "b2") {
        lat = Lattice::B2;
      } else if (str == "ch4") {
        lat = Lattice::CH4;
      } else if (str == "lin") {
        lat = Lattice::LIN;
      } else if (str == "zig") {
        lat = Lattice::ZIG;
      } else if (str == "tri") {
        lat = Lattice::TRI;
      } else {
        return true;
      }
    }

    return false;
  }

  double MEAMC::Rose(double const r,
                     double const re,
                     double const alpha,
                     double const Ec,
                     double const repuls,
                     double const attrac,
                     int const form) {
    if (r <= 0.0) {
      return 0.0;
    }

    // Eq.4.10d
    double const r_re = r / re;
    double const astar = alpha * (r_re - 1.0);

    double const astar_cube = special::Cube(astar);
    double const exp_astar = std::exp(-astar);

    double result;
    if (form == 1) {
      double const a3 = (-attrac + repuls / r);
      result = -Ec * (1 + astar + a3 * astar_cube) * exp_astar;
    } else if (form == 2) {
      double const a3 = (astar >= 0) ? attrac : repuls;
      result = -Ec * (1 + astar + a3 * astar_cube) * exp_astar;
    } else {
      double const a3 = (astar >= 0) ? attrac : repuls;
      result = -Ec * (1 + astar + a3 * astar_cube / r_re) * exp_astar;
    }
    return result;
  }

  double MEAMC::NumNearestNeighborsInReferenceStructure(Lattice const &lat) {
    // Number of nearest neighbors for the reference structure
    switch (lat) {
      case Lattice::FCC: {
        return 12.0;
      }
      case Lattice::BCC: {
        return 8.0;
      }
      case Lattice::HCP: {
        return 12.0;
      }
      case Lattice::DIA:
      case Lattice::DIA3: {
        return 4.0;
      }
      case Lattice::DIM: {
        return 1.0;
      }
      case Lattice::B1: {
        return 6.0;
      }
      case Lattice::C11: {
        return 10.0;
      }
      case Lattice::L12: {
        return 12.0;
      }
      case Lattice::B2: {
        return 8.0;
      }
      // DYNAMO currently implemented this way while it needs
      // two Z values, 4 and 1
      case Lattice::CH4: {
        return 4.0;
      }
      case Lattice::LIN:
      case Lattice::ZIG:
      case Lattice::TRI: {
        return 2.0;
      }
    }
    return 0.0;
  }

  double MEAMC::NumSecondNearestNeighborsInReferenceStructure(Lattice const &lat,
                                                              double const cmin,
                                                              double const cmax,
                                                              double const stheta,
                                                              double &distance_ratio,
                                                              double &second_neighbor_screening) {
    double number_of_second_nearest_neighbors_in_reference_structure = 0.0;
    int number_of_first_neighbors_screening_second_neighbors = 0;

    switch (lat) {
      case Lattice::FCC: {
        number_of_second_nearest_neighbors_in_reference_structure = 6.0;
        distance_ratio = kSqrtTwo;
        number_of_first_neighbors_screening_second_neighbors = 4;
        break;
      }
      case Lattice::BCC: {
        number_of_second_nearest_neighbors_in_reference_structure = 6.0;
        distance_ratio = 2.0 / kSqrtThree;
        number_of_first_neighbors_screening_second_neighbors = 4;
        break;
      }
      case Lattice::HCP: {
        number_of_second_nearest_neighbors_in_reference_structure = 6.0;
        distance_ratio = kSqrtTwo;
        number_of_first_neighbors_screening_second_neighbors = 4;
        break;
      }
      case Lattice::B1: {
        number_of_second_nearest_neighbors_in_reference_structure = 12.0;
        distance_ratio = kSqrtTwo;
        number_of_first_neighbors_screening_second_neighbors = 2;
        break;
      }
      case Lattice::DIA: {  // 2NN
        number_of_second_nearest_neighbors_in_reference_structure = 12.0;
        distance_ratio = std::sqrt(8.0 / 3.0);
        number_of_first_neighbors_screening_second_neighbors = 1;
        break;
      }
      case Lattice::DIA3: {  // 3NN
        number_of_second_nearest_neighbors_in_reference_structure = 12.0;
        distance_ratio = std::sqrt(11.0 / 3.0);
        number_of_first_neighbors_screening_second_neighbors = 4;
        break;
      }
      case Lattice::CH4:  // does not have 2nn structure so it returns 0
      case Lattice::LIN:  // line
      case Lattice::DIM: {
        // this really shouldn't be allowed; make sure screening is zero
        distance_ratio = 1.0;
        second_neighbor_screening = 0.0;
        return 0;
      }
      case Lattice::TRI: {  // TRI
        number_of_second_nearest_neighbors_in_reference_structure = 1.0;
        distance_ratio = 2.0 * stheta;
        number_of_first_neighbors_screening_second_neighbors = 2;
        break;
      }
      case Lattice::ZIG: {  // zig-zag
        number_of_second_nearest_neighbors_in_reference_structure = 2.0;
        distance_ratio = 2.0 * stheta;
        number_of_first_neighbors_screening_second_neighbors = 1;
        break;
      }
      case Lattice::L12: {
        number_of_second_nearest_neighbors_in_reference_structure = 6.0;
        distance_ratio = kSqrtTwo;
        number_of_first_neighbors_screening_second_neighbors = 4;
        break;
      }
      case Lattice::B2: {
        number_of_second_nearest_neighbors_in_reference_structure = 6.0;
        distance_ratio = 2.0 / kSqrtThree;
        number_of_first_neighbors_screening_second_neighbors = 4;
        break;
      }
      case Lattice::C11: {
        // unsupported lattice flag C11
        break;
      }
      default: {
        // unknown lattic flag
        break;
      }
    }

    // Compute screening for each first neighbor
    // special case for 3NN diamond structure
    double const c = (lat == Lattice::DIA3) ? 1.0 : 4.0 / (distance_ratio * distance_ratio) - 1.0;
    double const x = (c - cmin) / (cmax - cmin);
    double const s_ijk = RadialCutoff(x, fcut_function_form_);

    second_neighbor_screening = special::PowInt(s_ijk, number_of_first_neighbors_screening_second_neighbors);

    return number_of_second_nearest_neighbors_in_reference_structure;
  }

  double MEAMC::GGamma(double const gamma, int const ibar) const {
    switch (ibar) {
      case 0:
      case 4: {
        double const switchpoint = -gsmooth_factor_ / (gsmooth_factor_ + 1.0);
        double const g_gamma = (gamma < switchpoint) ? 1.0 / (gsmooth_factor_ + 1.0)
                                                           * std::pow((switchpoint / gamma), gsmooth_factor_)
                                                     : 1.0 + gamma;
        return std::sqrt(g_gamma);
      }
      case 1: {
        return std::exp(gamma * 0.5);
      }
      case 3: {
        double const g_gamma = 1.0 + std::exp(-gamma);
        return 2.0 / g_gamma;
      }
      case -5: {
        double const g_gamma = 1.0 + gamma;
        return (g_gamma >= 0) ? std::sqrt(g_gamma) : -std::sqrt(-g_gamma);
      }
    }
    return 0.0;
  }

  double MEAMC::GGamma(double const gamma, int const ibar, double &dg_gamma) const {
    switch (ibar) {
      case 0:
      case 4: {
        double const switchpoint = -gsmooth_factor_ / (gsmooth_factor_ + 1);
        if (gamma < switchpoint) {
          double g_gamma = 1 / (gsmooth_factor_ + 1);
          g_gamma *= std::pow(switchpoint / gamma, gsmooth_factor_);
          g_gamma = std::sqrt(g_gamma);
          dg_gamma = -gsmooth_factor_ * g_gamma / (2.0 * gamma);
          return g_gamma;
        } else {
          double const g_gamma = std::sqrt(1.0 + gamma);
          dg_gamma = 0.5 / g_gamma;
          return g_gamma;
        }
      }
      case 1: {
        double const g_gamma = std::exp(gamma * 0.5);
        dg_gamma = g_gamma * 0.5;
        return g_gamma;
      }
      case 3: {
        double const g_gamma = 2.0 / (1.0 + std::exp(-gamma));
        dg_gamma = g_gamma * (2.0 - g_gamma) * 0.5;
        return g_gamma;
      }
      case -5: {
        double const gamma_1 = 1.0 + gamma;
        if (gamma_1 >= 0) {
          double const g_gamma = std::sqrt(gamma_1);
          dg_gamma = 0.5 / g_gamma;
          return g_gamma;
        } else {
          double const g_gamma = -std::sqrt(-gamma_1);
          dg_gamma = -0.5 / g_gamma;
          return g_gamma;
        }
      }
    }
    dg_gamma = 1.0;
    return 0.0;
  }

  double MEAMC::Embedding(double const A, double const Ec, double const rhobar) const {
    return (rhobar > 0.0)                                 ? A * Ec * rhobar * std::log(rhobar)
           : (use_rhobar_linear_embedding_function_ == 0) ? 0.0
                                                          : -A * Ec * rhobar;
  }

  double MEAMC::Embedding(double const A, double const Ec, double const rhobar, double &embedding_df) const {
    if (rhobar > 0.0) {
      double const lrb = std::log(rhobar);
      double const AEc = A * Ec;
      embedding_df = AEc * (1.0 + lrb);
      return AEc * rhobar * lrb;
    }

    if (use_rhobar_linear_embedding_function_ == 0) {
      embedding_df = 0.0;
      return 0.0;
    }

    double const AEc = A * Ec;
    embedding_df = -AEc;
    return -AEc * rhobar;
  }

  double MEAMC::ComputePhi(double const r, int const a, int const b) {
    // Init phi
    double phi_m = 0.0;

    // if r = 0, just return 0
    if (special::IsZero(r)) {
      return phi_m;
    }

    double rho0_a;
    double rho1_a;
    double rho2_a;
    double rho3_a;

    double rho0_b;
    double rho1_b;
    double rho2_b;
    double rho3_b;

    // Calculate density functions, assuming reference configuration
    ComputeReferenceConfigurationDensity(
        r, a, b, &rho0_a, &rho1_a, &rho2_a, &rho3_a, &rho0_b, &rho1_b, &rho2_b, &rho3_b);

    // if densities are too small, numerical problems may result; just return zero
    if (special::IsZero(rho0_a, 1e-14) && special::IsZero(rho0_b, 1e-14)) {
      return phi_m;
    }

    Lattice const lat_a = element_lattice_(a, a);
    Lattice const lat_b = element_lattice_(b, b);
    Lattice const lat_ab = element_lattice_(a, b);

    // get number of neighbors in the reference structure
    // Nref[i][j] = # of i's neighbors of type j
    double const z_a = NumNearestNeighborsInReferenceStructure(lat_a);
    double const z_b = NumNearestNeighborsInReferenceStructure(lat_b);
    double const z_ab = NumNearestNeighborsInReferenceStructure(lat_ab);

    double t1_av_a;
    double t2_av_a;
    double t3_av_a;

    double t1_av_b;
    double t2_av_b;
    double t3_av_b;

    double rhobar_a;
    double rhobar_b;

    double g_a;
    double g_b;

    double rho_bkgd_a;
    double rho_bkgd_b;

    // Calculate average weighting factors for the reference structure

    if (ialloy_ == 2) {
      t1_av_a = element_t1_[a];
      t1_av_b = element_t1_[b];

      t2_av_a = element_t2_[a];
      t2_av_b = element_t2_[b];

      t3_av_a = element_t3_[a];
      t3_av_b = element_t3_[b];
    }
    // average weighting factors for the reference structure
    else if (lat_ab == Lattice::C11) {
      double const scalfac = 1.0 / (rho0_a + rho0_b);

      t1_av_a = scalfac * (element_t1_[a] * rho0_a + element_t1_[b] * rho0_b);
      t1_av_b = t1_av_a;

      t2_av_a = scalfac * (element_t2_[a] * rho0_a + element_t2_[b] * rho0_b);
      t2_av_b = t2_av_a;

      t3_av_a = scalfac * (element_t3_[a] * rho0_a + element_t3_[b] * rho0_b);
      t3_av_b = t3_av_a;

    } else if (lat_ab == Lattice::L12) {
      // In Eq.4.8
      double const rnorm_a = r / element_re_(a, a) - 1.0;
      double const rnorm_b = r / element_re_(b, b) - 1.0;

      // Eq.4.8
      double const rhoa0_a = element_rho0_[a] * std::exp(-element_beta0_[a] * rnorm_a);
      double const rhoa0_b = element_rho0_[b] * std::exp(-element_beta0_[b] * rnorm_b);

      double const scalfac = 8 * rhoa0_a + 4 * rhoa0_b;

      t1_av_a = (8 * element_t1_[a] * rhoa0_a + 4 * element_t1_[b] * rhoa0_b) / scalfac;
      t1_av_b = element_t1_[a];

      t2_av_a = (8 * element_t2_[a] * rhoa0_a + 4 * element_t2_[b] * rhoa0_b) / scalfac;
      t2_av_b = element_t2_[a];

      t3_av_a = (8 * element_t3_[a] * rhoa0_a + 4 * element_t3_[b] * rhoa0_b) / scalfac;
      t3_av_b = element_t3_[a];
    }
    // any of FCC, BCC, HCP, DIM, DIA, DIA3, B1, B2, CH4, LIN, ZIG, TRI
    else {
      // all neighbors are of the opposite type
      t1_av_a = element_t1_[b];
      t2_av_a = element_t2_[b];
      t3_av_a = element_t3_[b];

      t1_av_b = element_t1_[a];
      t2_av_b = element_t2_[a];
      t3_av_b = element_t3_[a];
    }

    // for c11b structure, calculate background electron
    // densities Eq.4.3
    if (lat_ab == Lattice::C11) {
      if (lat_a == Lattice::DIA) {
        // The int cast is happening here to be compliaent with LAMMPS
        rhobar_a = special::Square((static_cast<int>(z_ab) / 2) * (rho0_b + rho0_a))
                   + t1_av_a * special::Square(rho1_b - rho1_a)
                   + t2_av_a / 6.0 * special::Square(rho2_b + rho2_a)
                   + 121.0 / 40.0 * t3_av_a * special::Square(rho3_b - rho3_a);
        rhobar_a = std::sqrt(rhobar_a);

        rhobar_b = special::Square(z_ab * rho0_a) + kTwoThird * t2_av_a * special::Square(rho2_a);
        rhobar_b = std::sqrt(rhobar_b);
      } else {
        rhobar_a = special::Square(z_ab * rho0_b) + kTwoThird * t2_av_b * special::Square(rho2_b);
        rhobar_a = std::sqrt(rhobar_a);

        // The int cast is happening here to be compliaent with LAMMPS
        rhobar_b = special::Square((static_cast<int>(z_ab) / 2) * (rho0_a + rho0_b))
                   + t1_av_b * special::Square(rho1_a - rho1_b)
                   + t2_av_b / 6.0 * special::Square(rho2_a + rho2_b)
                   + 121.0 / 40.0 * t3_av_b * special::Square(rho3_a - rho3_b);
        rhobar_b = std::sqrt(rhobar_b);
      }
    } else {
      int const ibar_a = element_ibar_[a];
      int const ibar_b = element_ibar_[b];

      // for other structures, use Eq.4.5
      //
      // Compute the composition dependent electron density scaling Eq.4.5
      //
      // If using mixing rule for t, apply to reference structure;
      // else use precomputed values
      if (mixing_rule_compute_t_ == 1) {
        // If selection parameter for Gamma function for element a is negative
        if (ibar_a <= 0) {
          g_a = 1.0;
        } else {
          VectorOfSizeDIM shape_factors;

          GetShapeFactors(lat_a, element_stheta_(a, a), element_ctheta_(a, a), shape_factors);

          double const gam_ref
              = (shape_factors[0] * t1_av_a + shape_factors[1] * t2_av_a + shape_factors[2] * t3_av_a)
                / (z_a * z_a);

          g_a = GGamma(gam_ref, ibar_a);
        }

        // If selection parameter for Gamma function for element b is negative
        if (ibar_b <= 0) {
          g_b = 1.0;
        } else {
          VectorOfSizeDIM shape_factors;

          GetShapeFactors(lat_b, element_stheta_(b, b), element_ctheta_(b, b), shape_factors);

          double const gam_ref
              = (shape_factors[0] * t1_av_b + shape_factors[1] * t2_av_b + shape_factors[2] * t3_av_b)
                / (z_b * z_b);

          g_b = GGamma(gam_ref, ibar_b);
        }

        // Eq.4.5
        rho_bkgd_a = element_rho0_[a] * z_a * g_a;
        rho_bkgd_b = element_rho0_[b] * z_b * g_b;
      }
      // mixing_rule_compute_t_ != 1
      else {
        if (dynamo_reference_density_ == 1) {
          rho_bkgd_a = element_rho0_[a] * z_a;
          rho_bkgd_b = element_rho0_[b] * z_b;
        } else {
          rho_bkgd_a = element_ref_rho_[a];
          rho_bkgd_b = element_ref_rho_[b];
        }
      }

      if (special::IsZero(rho0_a, 1e-14)) {
        g_a = GGamma(0.0, ibar_a);
      } else {
        double const gam_ref = (t1_av_a * rho1_a + t2_av_a * rho2_a + t3_av_a * rho3_a) / (rho0_a * rho0_a);

        g_a = GGamma(gam_ref, ibar_a);
      }

      if (special::IsZero(rho0_b, 1e-14)) {
        g_b = GGamma(0.0, ibar_b);
      } else {
        double const gam_ref = (t1_av_b * rho1_b + t2_av_b * rho2_b + t3_av_b * rho3_b) / (rho0_b * rho0_b);

        g_b = GGamma(gam_ref, ibar_b);
      }

      // Eq.4.3
      rhobar_a = rho0_a / rho_bkgd_a * g_a;
      rhobar_b = rho0_b / rho_bkgd_b * g_b;
    }

    // compute the embedding energy, Eq.4.2
    double const embedding_a = Embedding(element_A_[a], element_Ec_(a, a), rhobar_a);

    double const embedding_b = Embedding(element_A_[b], element_Ec_(b, b), rhobar_b);

    // compute the energy of the reference state as a function
    // of interatomic spacing. Rose function, Eq.4.10c
    double const rose = Rose(r,
                             element_re_(a, b),
                             element_alpha_(a, b),
                             element_Ec_(a, b),
                             element_repuls_(a, b),
                             element_attrac_(a, b),
                             rose_function_form_);

    // calculate the pair energy
    if (lat_ab == Lattice::C11) {
      if (lat_a == Lattice::DIA) {
        double const phi_aa = ComputePhi(r, a, a);
        phi_m = (3. * rose - embedding_b - 2. * embedding_a - 5. * phi_aa) / z_ab;
      } else {
        double const phi_bb = ComputePhi(r, b, b);
        phi_m = (3. * rose - embedding_a - 2. * embedding_b - 5. * phi_bb) / z_ab;
      }

    } else if (lat_ab == Lattice::L12) {
      double distance_ratio;
      double second_neighbor_screening;

      double phi_aa = ComputePhi(r, a, a);

      // account for second neighbor a-a potential here...
      double const z_nn2 = NumSecondNearestNeighborsInReferenceStructure(lat_a,
                                                                         element_Cmin_(a, a, a),
                                                                         element_Cmax_(a, a, a),
                                                                         element_stheta_(a, b),
                                                                         distance_ratio,
                                                                         second_neighbor_screening);

      phi_aa += ComputePhiSeries(second_neighbor_screening, z_a, z_nn2, r, a, a, distance_ratio);

      phi_m = rose / 3.0 - embedding_a * 0.25 - embedding_b / 12.0 - phi_aa;

    } else if (lat_ab == Lattice::CH4) {
      phi_m = (5.0 * rose - embedding_a - 4.0 * embedding_b) * 0.25;

    } else if (lat_ab == Lattice::ZIG) {
      phi_m = 2.0 * rose - embedding_a - embedding_b;

      if (a != b) {
        double const s_aab = Sijk(1.0, a, a, b);
        double const b11s = -2.0 / z_ab * s_aab;

        double const s_bba = Sijk(1.0, b, b, a);
        double const b22s = -2.0 / z_ab * s_bba;

        double const phi_aa = ComputePhi(2.0 * element_stheta_(a, b) * r, a, a);
        double const phi_bb = ComputePhi(2.0 * element_stheta_(a, b) * r, b, b);

        phi_m += phi_aa * b11s + phi_bb * b22s;
      }

      phi_m /= z_ab;

    } else if (lat_ab == Lattice::TRI) {
      phi_m = 3.0 * rose - 2.0 * embedding_a - embedding_b;

      if (a != b) {
        double const s_aab = Sijk(1.0, a, a, b);
        double const b11s = -2.0 / z_ab * s_aab;

        double const phi_aa = ComputePhi(2.0 * element_stheta_(a, b) * r, a, a);

        phi_m += phi_aa * b11s;
      }

      phi_m /= z_ab;

    } else {
      // potential is computed from Rose function and embedding energy
      phi_m = (2.0 * rose - embedding_a - embedding_b) / z_ab;
    }

    return phi_m;
  }

  double MEAMC::ComputePhiSeries(double const second_neighbor_screening,
                                 double const z1,
                                 double const z2,
                                 double const r,
                                 int const a,
                                 int const b,
                                 double const distance_ratio) {
    double phi_sum = 0.0;
    if (second_neighbor_screening > 0.0) {
      double const b_nn2 = -z2 * second_neighbor_screening / z1;
      for (int n = 1; n <= 10; ++n) {
        double const c1 = special::PowInt(b_nn2, n);
        double const c2 = special::PowInt(distance_ratio, n);
        double const c3 = ComputePhi(r * c2, a, b);
        double const phi_val = c1 * c3;
        if (special::IsZero(phi_val)) {
          // once either term becomes zero at some point, all folliwng will also
          // be zero necessary to avoid numerical error (nan or infty) due to
          // exponential decay in ComputePhi
          break;
        }
        phi_sum += phi_val;
      }
    }
    return phi_sum;
  }

  void MEAMC::CheckIndex(int const num, int const lim, int const nidx, int const *idx, int *ierr) {
    if (nidx < num) {
      *ierr = 2;
      return;
    }

    *ierr = 0;
    for (int i = 0; i < num; ++i) {
      if ((idx[i] < 0) || (idx[i] >= lim)) {
        *ierr = 3;
        return;
      }
    }
  }

  void MEAMC::SetParameter(int const which,
                           double const value,
                           int const nindex,
                           int const *index,
                           int *errorflag) {
    *errorflag = 0;

    switch (which) {
      // 0 = Ec
      case 0: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_Ec_(index[0], index[1]) = value;
        break;
      }
      // 1 = alpha
      case 1: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_alpha_(index[0], index[1]) = value;
        break;
      }
      // 2 = rho0
      case 2: {
        CheckIndex(1, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_rho0_[index[0]] = value;
        break;
      }
      // 3 = delta
      case 3: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_delta_(index[0], index[1]) = value;
        break;
      }
      // 4 = lattce
      case 4: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_lattice_(index[0], index[1]) = static_cast<Lattice>(value);
        break;
      }
      // 5 = attrac
      case 5: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_attrac_(index[0], index[1]) = value;
        break;
      }
      // 6 = repuls
      case 6: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_repuls_(index[0], index[1]) = value;
        break;
      }
      // 7 = nn2
      case 7: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        int const i1 = std::min(index[0], index[1]);
        int const i2 = std::max(index[0], index[1]);
        element_nn2_(i1, i2) = static_cast<int>(value);
        break;
      }
      // 8 = Cmin
      case 8: {
        CheckIndex(3, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_Cmin_(index[0], index[1], index[2]) = value;
        break;
      }
      // 9 = Cmax
      case 9: {
        CheckIndex(3, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_Cmax_(index[0], index[1], index[2]) = value;
        break;
      }
      // 10 = rc
      case 10: {
        cutoff_radius_ = value;
        break;
      }
      // 11 = delr
      case 11: {
        delr_ = value;
        break;
      }
      // 12 = augt1
      case 12: {
        augment_t1_ = static_cast<int>(value);
        break;
      }
      // 13 = gsmooth_factor
      case 13: {
        gsmooth_factor_ = value;
        break;
      }
      // 14 = re
      case 14: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        element_re_(index[0], index[1]) = value;
        break;
      }
      // 15 = ialloy
      case 15: {
        ialloy_ = static_cast<int>(value);
        break;
      }
      // 16 = mixture_ref_t
      case 16: {
        mixing_rule_compute_t_ = static_cast<int>(value);
        break;
      }
      // 17 = erose_form
      case 17: {
        rose_function_form_ = static_cast<int>(value);
        if (!(rose_function_form_ == 0 || rose_function_form_ == 1 || rose_function_form_ == 2)) {
          *errorflag = 5;
          return;
        }
        break;
      }
      // 18 = zbl
      case 18: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        int const i1 = std::min(index[0], index[1]);
        int const i2 = std::max(index[0], index[1]);
        element_zbl_(i1, i2) = static_cast<int>(value);
        break;
      }
      // 19 = emb_lin_neg
      case 19: {
        use_rhobar_linear_embedding_function_ = static_cast<int>(value);
        break;
      }
      // 20 = bkgd_dyn
      case 20: {
        dynamo_reference_density_ = static_cast<int>(value);
        break;
      }
      // 21 = theta
      case 21: {
        CheckIndex(2, number_of_element_types_, nindex, index, errorflag);
        if (*errorflag != 0) {
          return;
        }
        int const i1 = std::min(index[0], index[1]);
        int const i2 = std::max(index[0], index[1]);
        element_theta_(i1, i2) = value;
        break;
      }
      // 22 = fcut_form
      case 22: {
        fcut_function_form_ = static_cast<int>(value);
        if (!(fcut_function_form_ == 0 || fcut_function_form_ == 1)) {
          *errorflag = 4;
          return;
        }
        break;
      }
      default: {
        *errorflag = 1;
      }
    }
  }

  void MEAMC::FillOffDiagonalAlloyParameters() {
    // Loop over pairs
    for (int i = 0; i < number_of_element_types_; ++i) {
      for (int j = 0; j < number_of_element_types_; ++j) {
        // Treat off-diagonal pairs
        // If i>j, set all equal to i<j case (which has already been set,
        // here or in the input file)
        if (i > j) {
          element_re_(i, j) = element_re_(j, i);
          element_Ec_(i, j) = element_Ec_(j, i);
          element_alpha_(i, j) = element_alpha_(j, i);
          element_lattice_(i, j) = element_lattice_(j, i);
          element_nn2_(i, j) = element_nn2_(j, i);
          // theta for lin,tri,zig references
          element_stheta_(i, j) = element_stheta_(j, i);
          element_ctheta_(i, j) = element_ctheta_(j, i);
          // If i<j and term is unset, use default values
          // (e.g. mean of i-i and j-j)
        } else if (i < j) {
          if (special::IsZero(element_Ec_(i, j))) {
            if (element_lattice_(i, j) == Lattice::L12) {
              element_Ec_(i, j) = (3 * element_Ec_(i, i) + element_Ec_(j, j)) * 0.25 - element_delta_(i, j);
            } else if (element_lattice_(i, j) == Lattice::C11) {
              if (element_lattice_(i, i) == Lattice::DIA) {
                element_Ec_(i, j) = (2 * element_Ec_(i, i) + element_Ec_(j, j)) / 3.0 - element_delta_(i, j);
              } else {
                element_Ec_(i, j) = (element_Ec_(i, i) + 2 * element_Ec_(j, j)) / 3.0 - element_delta_(i, j);
              }
            } else {
              element_Ec_(i, j) = (element_Ec_(i, i) + element_Ec_(j, j)) * 0.5 - element_delta_(i, j);
            }
          }
          if (special::IsZero(element_alpha_(i, j))) {
            element_alpha_(i, j) = (element_alpha_(i, i) + element_alpha_(j, j)) * 0.5;
          }
          if (special::IsZero(element_re_(i, j))) {
            element_re_(i, j) = (element_re_(i, i) + element_re_(j, j)) * 0.5;
          }
        }
      }
    }

    // cmin[i][j][k] is symmetric in i-j, but not k.  For all triplets
    // where i>j, set equal to the i<j element.  Likewise for cmax.
    for (int i = 1; i < number_of_element_types_; ++i) {
      for (int j = 0; j < i; ++j) {
        for (int k = 0; k < number_of_element_types_; ++k) {
          element_Cmin_(i, j, k) = element_Cmin_(j, i, k);
          element_Cmax_(i, j, k) = element_Cmax_(j, i, k);
        }
      }
    }

    // element_ebound_ gives the squared distance such that, for rik2 or rjk2 >
    // element_ebound_, atom k definitely lies outside the screening function
    // ellipse (so there is no need to calculate its effects).  Here, compute it
    // for all triplets [i][j][k] so that element_ebound_[i][j] is the maximized
    // over k
    for (int i = 0; i < number_of_element_types_; ++i) {
      for (int j = 0; j < number_of_element_types_; ++j) {
        for (int k = 0; k < number_of_element_types_; ++k) {
          double const cmax = element_Cmax_(i, j, k);
          double const eb = (cmax * cmax) / (4.0 * (cmax - 1.0));
          element_ebound_(i, j) = std::max(element_ebound_(i, j), eb);
        }
      }
    }
  }

  void MEAMC::ComputeScreeningAndDerivative(int const i,
                                            int const number_of_neighbors,
                                            int const *const neighbors_of_particle,
                                            int const offset,
                                            const VectorOfSizeDIM *const coordinates,
                                            int const *const particle_species_codes,
                                            int const *const particle_contributing) {
    double *const scrfcn = &scrfcn_[offset];
    double *const dscrfcn = &dscrfcn_[offset];

    int const species_i = particle_species_codes[i];

    // delr_ controls the distance over which the
    // radial cutoff is smoothed from 1 to 0 near r = r_c
    double const delr_inv = 1.0 / delr_;

    double const xi = coordinates[i][0];
    double const yi = coordinates[i][1];
    double const zi = coordinates[i][2];

    // Setup loop over neighbors of current particle
    for (int jn = 0, nth_neighbor_of_i = -1; jn < number_of_neighbors; ++jn) {
      int const j = neighbors_of_particle[jn];

      // Effective half list
      if (!(particle_contributing[j] && j < i)) {
        ++nth_neighbor_of_i;
        int const species_j = particle_species_codes[j];

        // First compute screening function itself, sij
        double const xj = coordinates[j][0];
        double const yj = coordinates[j][1];
        double const zj = coordinates[j][2];

        double const dx_ij = xj - xi;
        double const dy_ij = yj - yi;
        double const dz_ij = zj - zi;

        double const rij2 = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

        // atoms i and j are outside the cutoff radius
        if (rij2 > cutoff_radius_squared_) {
          scrfcn[nth_neighbor_of_i] = 0.0;
          dscrfcn[nth_neighbor_of_i] = 0.0;
          continue;
        }

        double const rbound = element_ebound_(species_i, species_j) * rij2;

        double const rij = std::sqrt(rij2);

        double const rnorm = (cutoff_radius_ - rij) * delr_inv;

        double sij = 1.0;

        for (int kn = 0; kn < number_of_neighbors; ++kn) {
          int const k = neighbors_of_particle[kn];
          if (k == j) {
            continue;
          }

          double const xk = coordinates[k][0];
          double const yk = coordinates[k][1];
          double const zk = coordinates[k][2];

          double const dx_jk = xk - xj;
          double const dy_jk = yk - yj;
          double const dz_jk = zk - zj;

          double const rjk2 = dx_jk * dx_jk + dy_jk * dy_jk + dz_jk * dz_jk;

          // definitely lies outside the screening function ellipse
          // so there is no need to calculate its effects
          if (rjk2 > rbound) {
            continue;
          }

          double const dx_ik = xk - xi;
          double const dy_ik = yk - yi;
          double const dz_ik = zk - zi;

          double const rik2 = dx_ik * dx_ik + dy_ik * dy_ik + dz_ik * dz_ik;

          // definitely lies outside the screening function ellipse
          // so there is no need to calculate its effects
          if (rik2 > rbound) {
            continue;
          }

          double const x_ik = rik2 / rij2;
          double const x_jk = rjk2 / rij2;

          double const x_ijk = x_ik - x_jk;

          double const a = 1 - x_ijk * x_ijk;

          // if a < 0, then ellipse equation doesn't describe this case
          // and atom k can't possibly screen i-j
          if (a <= 0.0) {
            continue;
          }

          // Eq.4.11d
          double cijk = (2.0 * (x_ik + x_jk) + a - 2.0) / a;

          int const species_k = particle_species_codes[k];

          double const cmax = element_Cmax_(species_i, species_j, species_k);

          // RadialCutoff(x) is 1, Eq.4.11e, x >= 1
          if (cijk >= cmax) {
            continue;
          }

          double const cmin = element_Cmin_(species_i, species_j, species_k);

          // note that cijk may be slightly negative (within numerical
          // tolerance) if atoms are colinear, so don't reject that case here
          // (other negative cijk cases were handled by the test on "a" above)

          // fc(x) is 0, Eq.4.11e, x <= 0
          if (cijk <= cmin) {
            sij = 0.0;
            break;
            // compute fc(x), Eq.4.11e, 0 < x < 1
          } else {
            double const delc = cmax - cmin;
            cijk = (cijk - cmin) / delc;

            // Eq.4.11c
            double const s_ijk = RadialCutoff(cijk, fcut_function_form_);
            sij *= s_ijk;
          }
        }  // Loop over k

        double dfc_ij;
        double const fc_ij = RadialCutoff(rnorm, dfc_ij, fcut_function_form_);
        dfc_ij *= delr_inv;

        double const sijfcij = sij * fc_ij;

        // Screening function, Eq.4.11a
        scrfcn[nth_neighbor_of_i] = sijfcij;

        // Now compute derivatives
        dscrfcn[nth_neighbor_of_i] = 0.0;

        if (special::IsNotZero(sijfcij) && special::IsNotOne(sijfcij)) {
          for (int kn = 0; kn < number_of_neighbors; ++kn) {
            int const k = neighbors_of_particle[kn];
            if (k == j) {
              continue;
            }

            double const xk = coordinates[k][0];
            double const yk = coordinates[k][1];
            double const zk = coordinates[k][2];

            double const dx_jk = xk - xj;
            double const dy_jk = yk - yj;
            double const dz_jk = zk - zj;

            double const rjk2 = dx_jk * dx_jk + dy_jk * dy_jk + dz_jk * dz_jk;

            // definitely lies outside the screening function ellipse
            // so there is no need to calculate its effects
            if (rjk2 > rbound) {
              continue;
            }

            double const dx_ik = xk - xi;
            double const dy_ik = yk - yi;
            double const dz_ik = zk - zi;

            double const rik2 = dx_ik * dx_ik + dy_ik * dy_ik + dz_ik * dz_ik;

            // definitely lies outside the screening function ellipse
            // so there is no need to calculate its effects
            if (rik2 > rbound) {
              continue;
            }

            double const x_ik = rik2 / rij2;
            double const x_jk = rjk2 / rij2;

            double const x_ijk = x_ik - x_jk;

            double const a = 1 - x_ijk * x_ijk;

            // if a < 0, then ellipse equation doesn't describe this case and
            // atom k can't possibly screen i-j
            if (a <= 0.0) {
              continue;
            }

            double cijk = (2.0 * (x_ik + x_jk) + a - 2.0) / a;

            int const species_k = particle_species_codes[k];

            double const cmax = element_Cmax_(species_i, species_j, species_k);

            // fc(x) is 1, Eq.4.11e, x >= 1
            if (cijk >= cmax) {
              continue;
            }

            // Note that cijk may be slightly negative (within numerical
            // tolerance) if atoms are colinear, so don't reject that case
            // here (other negative cijk cases were handled by the test
            // on "a" above)
            //
            // Note that we never have 0 < cijk < cmin here, else sij=0
            // (rejected above)
            {
              double const cmin = element_Cmin_(species_i, species_j, species_k);

              double const delc = cmax - cmin;
              cijk = (cijk - cmin) / delc;

              double df_ijk;
              double const s_ijk = RadialCutoff(cijk, df_ijk, fcut_function_form_);

              double const coef1 = df_ijk / (delc * s_ijk);

              double const dc_ij = DCijkDRij(rij2, rik2, rjk2);

              dscrfcn[nth_neighbor_of_i] += coef1 * dc_ij;
            }
          }  // Loop over k

          double const coef2 = sij * dfc_ij / rij;

          // Eq.4.22a
          dscrfcn[nth_neighbor_of_i] = dscrfcn[nth_neighbor_of_i] * sijfcij - coef2;
        }  // if (special::IsNotZero(sijfcij) && special::IsNotOne(sijfcij))
      }    // Effective half list
    }      // Loop over j
  }

  void MEAMC::ComputeIntermediateDensityTerms(int const i,
                                              int const number_of_neighbors,
                                              int const *const neighbors_of_particle,
                                              int &offset,
                                              const VectorOfSizeDIM *const coordinates,
                                              int const *const particle_species_codes,
                                              int const *const particle_contributing) {
    double const *const scrfcn = &scrfcn_[offset];

    int const species_i = particle_species_codes[i];

    double const xi = coordinates[i][0];
    double const yi = coordinates[i][1];
    double const zi = coordinates[i][2];

    int nth_neighbor_of_i = -1;

    // Setup loop over neighbors of current particle
    for (int jn = 0; jn < number_of_neighbors; ++jn) {
      int const j = neighbors_of_particle[jn];

      // Effective half list
      if (!(particle_contributing[j] && j < i)) {
        ++nth_neighbor_of_i;

        // Screening function.
        double const sij = scrfcn[nth_neighbor_of_i];

        if (special::IsNotZero(sij)) {
          double const dx = coordinates[j][0] - xi;
          double const dy = coordinates[j][1] - yi;
          double const dz = coordinates[j][2] - zi;

          double const dxdx = dx * dx;
          double const dydy = dy * dy;
          double const dzdz = dz * dz;

          double const rij2 = dxdx + dydy + dzdz;

          if (rij2 < cutoff_radius_squared_) {
            int const species_j = particle_species_codes[j];

            double const rij = std::sqrt(rij2);

            // (rij/r0i - 1) in Eq.4.8
            double const aj = rij / element_re_(species_j, species_j) - 1.0;

            // Element-dependent density scaling in Eq.4.5 or Eq.4.8
            double const rhoj0 = element_rho0_[species_j];

            // The atomic electron densities. Eq.4.8
            double rhoa0_j = rhoj0 * std::exp(-element_beta0_[species_j] * aj) * sij;
            double rhoa1_j = rhoj0 * std::exp(-element_beta1_[species_j] * aj) * sij;
            double rhoa2_j = rhoj0 * std::exp(-element_beta2_[species_j] * aj) * sij;
            double rhoa3_j = rhoj0 * std::exp(-element_beta3_[species_j] * aj) * sij;

            // Eq.4.7a
            rho0_[i] += rhoa0_j;

            double const t1_j = element_t1_[species_j];
            double const t2_j = element_t2_[species_j];
            double const t3_j = element_t3_[species_j];

            if (ialloy_ == 1) {
              rhoa1_j *= t1_j;
              rhoa2_j *= t2_j;
              rhoa3_j *= t3_j;

              double *const tsq_ave_1d_i = &tsq_ave_(i, 0);

              tsq_ave_1d_i[0] += t1_j * t1_j * rhoa0_j;
              tsq_ave_1d_i[1] += t2_j * t2_j * rhoa0_j;
              tsq_ave_1d_i[2] += t3_j * t3_j * rhoa0_j;
            }

            // For ialloy_ = 2, use single-element value (not average)
            if (ialloy_ != 2) {
              double *const t_ave_1d_i = &t_ave_(i, 0);

              // Eq.4.9, \sum_j {tk0j * rhoa0_j}
              t_ave_1d_i[0] += t1_j * rhoa0_j;
              t_ave_1d_i[1] += t2_j * rhoa0_j;
              t_ave_1d_i[2] += t3_j * rhoa0_j;
            }

            double *const arho1_1d_i = &arho1_(i, 0);

            // In Eq.4.7b
            double const A1j = rhoa1_j / rij;
            arho1_1d_i[0] += A1j * dx;
            arho1_1d_i[1] += A1j * dy;
            arho1_1d_i[2] += A1j * dz;

            double const dxdy = dx * dy;
            double const dxdz = dx * dz;
            double const dydz = dy * dz;

            double *const arho2_1d_i = &arho2_(i, 0);

            // In Eq.4.7c, the first part
            double const A2j = rhoa2_j / rij2;
            arho2_1d_i[0] += A2j * dxdx;
            arho2_1d_i[1] += A2j * dxdy;
            arho2_1d_i[2] += A2j * dxdz;
            arho2_1d_i[3] += A2j * dydy;
            arho2_1d_i[4] += A2j * dydz;
            arho2_1d_i[5] += A2j * dzdz;

            // In Eq.4.7c, the second part
            arho2b_[i] += rhoa2_j;

            double const dxdxdx = dxdx * dx;
            double const dxdxdy = dxdx * dy;
            double const dxdxdz = dxdx * dz;
            double const dxdydy = dxdy * dy;
            double const dxdydz = dxdy * dz;
            double const dxdzdz = dxdz * dz;
            double const dydydy = dydy * dy;
            double const dydydz = dydy * dz;
            double const dydzdz = dydz * dz;
            double const dzdzdz = dzdz * dz;

            double *const arho3_1d_i = &arho3_(i, 0);

            // In Eq.4.7d, the first part
            double const A3j = rhoa3_j / (rij2 * rij);
            arho3_1d_i[0] += A3j * dxdxdx;
            arho3_1d_i[1] += A3j * dxdxdy;
            arho3_1d_i[2] += A3j * dxdxdz;
            arho3_1d_i[3] += A3j * dxdydy;
            arho3_1d_i[4] += A3j * dxdydz;
            arho3_1d_i[5] += A3j * dxdzdz;
            arho3_1d_i[6] += A3j * dydydy;
            arho3_1d_i[7] += A3j * dydydz;
            arho3_1d_i[8] += A3j * dydzdz;
            arho3_1d_i[9] += A3j * dzdzdz;

            double *const arho3b_1d_i = &arho3b_(i, 0);

            // In Eq.4.7d, the second part
            arho3b_1d_i[0] += rhoa3_j * dx / rij;
            arho3b_1d_i[1] += rhoa3_j * dy / rij;
            arho3b_1d_i[2] += rhoa3_j * dz / rij;

            if (particle_contributing[j]) {
              double ai;
              double rhoa0_i;
              double rhoa1_i;
              double rhoa2_i;
              double rhoa3_i;

              if (species_i == species_j) {
                ai = aj;
                rhoa0_i = rhoa0_j;
                rhoa1_i = rhoa1_j;
                rhoa2_i = rhoa2_j;
                rhoa3_i = rhoa3_j;
              } else {
                // Element-dependent density scaling in Eq.4.5 or Eq.4.8
                double const rhoi0 = element_rho0_[species_i];

                // (rij/r0i - 1) in Eq.4.8
                ai = rij / element_re_(species_i, species_i) - 1.0;

                // The atomic electron densities. Eq.4.8
                rhoa0_i = rhoi0 * std::exp(-element_beta0_[species_i] * ai) * sij;
                rhoa1_i = rhoi0 * std::exp(-element_beta1_[species_i] * ai) * sij;
                rhoa2_i = rhoi0 * std::exp(-element_beta2_[species_i] * ai) * sij;
                rhoa3_i = rhoi0 * std::exp(-element_beta3_[species_i] * ai) * sij;
              }

              // Eq.4.7a
              rho0_[j] += rhoa0_i;

              double const t1_i = element_t1_[species_i];
              double const t2_i = element_t2_[species_i];
              double const t3_i = element_t3_[species_i];

              if (ialloy_ == 1) {
                if (species_i != species_j) {
                  rhoa1_i *= t1_i;
                  rhoa2_i *= t2_i;
                  rhoa3_i *= t3_i;
                }

                double *const tsq_ave_1d_j = &tsq_ave_(j, 0);

                tsq_ave_1d_j[0] += t1_i * t1_i * rhoa0_i;
                tsq_ave_1d_j[1] += t2_i * t2_i * rhoa0_i;
                tsq_ave_1d_j[2] += t3_i * t3_i * rhoa0_i;
              }

              // For ialloy_ = 2, use single-element value (not average)
              if (ialloy_ != 2) {
                double *const t_ave_1d_j = &t_ave_(j, 0);

                // Eq.4.9
                t_ave_1d_j[0] += t1_i * rhoa0_i;
                t_ave_1d_j[1] += t2_i * rhoa0_i;
                t_ave_1d_j[2] += t3_i * rhoa0_i;
              }

              double *const arho1_1d_j = &arho1_(j, 0);

              // In Eq.4.7b
              double const A1i = rhoa1_i / rij;
              arho1_1d_j[0] -= A1i * dx;
              arho1_1d_j[1] -= A1i * dy;
              arho1_1d_j[2] -= A1i * dz;

              double *const arho2_1d_j = &arho2_(j, 0);

              // In Eq.4.7c, the first part
              double const A2i = rhoa2_i / rij2;
              arho2_1d_j[0] += A2i * dxdx;
              arho2_1d_j[1] += A2i * dxdy;
              arho2_1d_j[2] += A2i * dxdz;
              arho2_1d_j[3] += A2i * dydy;
              arho2_1d_j[4] += A2i * dydz;
              arho2_1d_j[5] += A2i * dzdz;

              // In Eq.4.7c, the second part
              arho2b_[j] += rhoa2_i;

              double *const arho3_1d_j = &arho3_(j, 0);

              // In Eq.4.7d, the first part
              double const A3i = rhoa3_i / (rij2 * rij);
              arho3_1d_j[0] -= A3i * dxdxdx;
              arho3_1d_j[1] -= A3i * dxdxdy;
              arho3_1d_j[2] -= A3i * dxdxdz;
              arho3_1d_j[3] -= A3i * dxdydy;
              arho3_1d_j[4] -= A3i * dxdydz;
              arho3_1d_j[5] -= A3i * dxdzdz;
              arho3_1d_j[6] -= A3i * dydydy;
              arho3_1d_j[7] -= A3i * dydydz;
              arho3_1d_j[8] -= A3i * dydzdz;
              arho3_1d_j[9] -= A3i * dzdzdz;

              double *const arho3b_1d_j = &arho3b_(j, 0);

              // In Eq.4.7d, the second part
              arho3b_1d_j[0] -= rhoa3_i * dx / rij;
              arho3b_1d_j[1] -= rhoa3_i * dy / rij;
              arho3b_1d_j[2] -= rhoa3_i * dz / rij;
            }  // particle_contributing[j]
          }    // (rij2 < cutoff_radius_squared_)
        }      // (special::IsNotZero(sij))
      }        // Effective half list
    }          // Loop over j

    if (nth_neighbor_of_i > -1) {
      offset += nth_neighbor_of_i + 1;
    }
  }

  void MEAMC::ResizePairPotentialArrays() {
    std::size_t const n1
        = static_cast<std::size_t>(number_of_element_types_ * (number_of_element_types_ + 1) / 2);

    // nr_ is a pair function discretization parameters
    std::size_t const n2 = static_cast<std::size_t>(nr_);

    // allocate memory for array that defines the potential
    phir_.resize(n1, n2);

    // allocate coeff memory
    phirar1_.resize(n1, n2);
    phirar2_.resize(n1, n2);
    phirar3_.resize(n1, n2);
    phirar4_.resize(n1, n2);
    phirar5_.resize(n1, n2);
    phirar6_.resize(n1, n2);
  }

  void MEAMC::ComputePairPotential() {
    ResizePairPotentialArrays();

    double distance_ratio;
    double second_neighbor_screening;

    double z1;
    double z2;

    // loop over pairs of element types
    for (int a = 0, nv2 = 0; a < number_of_element_types_; ++a) {
      for (int b = a; b < number_of_element_types_; ++b, ++nv2) {
        // Loop over r values and compute
        for (int j = 0; j < nr_; ++j) {
          double const r = j * dr_;

          phir_(nv2, j) = ComputePhi(r, a, b);

          // if using second-nearest neighbor, solve recursive problem
          // (see Lee and Baskes, PRB 62(13):8564 eqn.(21))
          if (element_nn2_(a, b) == 1) {
            Lattice const lat_aa = element_lattice_(a, a);
            Lattice const lat_ab = element_lattice_(a, b);
            Lattice const lat_bb = element_lattice_(b, b);

            z2 = NumSecondNearestNeighborsInReferenceStructure(lat_ab,
                                                               element_Cmin_(a, a, b),
                                                               element_Cmax_(a, a, b),
                                                               element_stheta_(a, b),
                                                               distance_ratio,
                                                               second_neighbor_screening);

            // The B1, B2, and L12 cases with NN2 have a trick to them; we
            // need to compute the contributions from second nearest
            // neighbors, like a-a pairs, but need to include NN2
            // contributions to those pairs as well.
            if (lat_ab == Lattice::B1 || lat_ab == Lattice::B2 || lat_ab == Lattice::L12
                || lat_ab == Lattice::DIA) {
              double const rarat = r * distance_ratio;

              // phi_aa
              double phi_aa = ComputePhi(rarat, a, a);

              z1 = NumNearestNeighborsInReferenceStructure(lat_aa);
              z2 = NumSecondNearestNeighborsInReferenceStructure(lat_aa,
                                                                 element_Cmin_(a, a, a),
                                                                 element_Cmax_(a, a, a),
                                                                 element_stheta_(a, a),
                                                                 distance_ratio,
                                                                 second_neighbor_screening);

              phi_aa += ComputePhiSeries(second_neighbor_screening, z1, z2, rarat, a, a, distance_ratio);

              // phi_bb
              double phi_bb = ComputePhi(rarat, b, b);

              z1 = NumNearestNeighborsInReferenceStructure(lat_bb);
              z2 = NumSecondNearestNeighborsInReferenceStructure(lat_bb,
                                                                 element_Cmin_(b, b, b),
                                                                 element_Cmax_(b, b, b),
                                                                 element_stheta_(b, b),
                                                                 distance_ratio,
                                                                 second_neighbor_screening);

              phi_bb += ComputePhiSeries(second_neighbor_screening, z1, z2, rarat, b, b, distance_ratio);

              if (lat_ab == Lattice::B1 || lat_ab == Lattice::B2 || lat_ab == Lattice::DIA) {
                // Add contributions to the B1 or B2 potential
                z1 = NumNearestNeighborsInReferenceStructure(lat_ab);
                z2 = NumSecondNearestNeighborsInReferenceStructure(lat_ab,
                                                                   element_Cmin_(a, a, b),
                                                                   element_Cmax_(a, a, b),
                                                                   element_stheta_(a, b),
                                                                   distance_ratio,
                                                                   second_neighbor_screening);

                phir_(nv2, j) -= z2 * second_neighbor_screening / (2.0 * z1) * phi_aa;

                z2 = NumSecondNearestNeighborsInReferenceStructure(lat_ab,
                                                                   element_Cmin_(b, b, a),
                                                                   element_Cmax_(b, b, a),
                                                                   element_stheta_(a, b),
                                                                   distance_ratio,
                                                                   second_neighbor_screening);

                phir_(nv2, j) -= z2 * second_neighbor_screening / (2.0 * z1) * phi_bb;

              }
              // lat_ab == Lattice::L12
              else {
                // The L12 case has one last trick; we have to be careful to
                // compute the correct screening between 2nd-neighbor pairs.
                // 1-1 second-neighbor pairs are screened by 2 type 1 atoms
                // and two type 2 atoms. 2-2 second-neighbor pairs are
                // screened by 4 type 1 atoms.
                double const s_aaa = Sijk(1.0, a, a, a);
                double const s_aab = Sijk(1.0, a, a, b);
                double const s_bba = Sijk(1.0, b, b, a);

                double const S11 = s_aaa * s_aaa * s_aab * s_aab;
                double const S22 = special::PowInt(s_bba, 4);

                phir_(nv2, j) -= 0.75 * S11 * phi_aa + 0.25 * S22 * phi_bb;
              }

            } else {
              z1 = NumNearestNeighborsInReferenceStructure(lat_ab);

              phir_(nv2, j) += ComputePhiSeries(second_neighbor_screening, z1, z2, r, a, b, distance_ratio);
            }
          }  // (element_nn2_(a, b) == 1)

          // For Zbl potential:
          if (element_zbl_(a, b) == 1) {
            // Eq.4.10d
            double const astar = element_alpha_(a, b) * (r / element_re_(a, b) - 1.0);

            // potential is Zbl potential
            if (astar <= -3.0) {
              zbl_->SetCoeff(a, b, element_atomic_number_[a], element_atomic_number_[b]);

              phir_(nv2, j) = zbl_->ComputePhi(r, a, b);
            }
            // potential is linear combination with Zbl potential
            else if (astar > -3.0 && astar < -1.0) {
              zbl_->SetCoeff(a, b, element_atomic_number_[a], element_atomic_number_[b]);

              double const frac = RadialCutoff(1 - (astar + 1.0) / (-3.0 + 1.0), fcut_function_form_);
              double const phizbl = zbl_->ComputePhi(r, a, b);
              phir_(nv2, j) = frac * phir_(nv2, j) + (1 - frac) * phizbl;
            }
          }  // (element_zbl_(a, b) == 1)
        }    // Loop over j

        // call interpolation
        SplineInterpolate(nv2);
      }  // Loop over b element
    }    // Loop over a element
  }

  void MEAMC::ComputeCompositionDependentDensityScaling() {
    double g_a;
    VectorOfSizeDIM shape_factors;

    // loop over element types
    for (int a = 0; a < number_of_element_types_; ++a) {
      // Zi0 in Eq.4.5
      double const z1 = NumNearestNeighborsInReferenceStructure(element_lattice_(a, a));

      // The zeroth order density in the reference structure, with
      // equilibrium spacing, is just the number of first neighbors times
      // the element_rho0_ coefficient...
      double rho0a = element_rho0_[a] * z1;

      // If selection parameter for Gamma function for element a is negative
      if (element_ibar_[a] <= 0) {
        g_a = 1.0;
      } else {
        GetShapeFactors(element_lattice_(a, a), element_stheta_(a, a), element_ctheta_(a, a), shape_factors);

        double const gam_ref = (element_t1_[a] * shape_factors[0] + element_t2_[a] * shape_factors[1]
                                + element_t3_[a] * shape_factors[2])
                               / (z1 * z1);

        g_a = GGamma(gam_ref, element_ibar_[a]);
      }

      // ...unless we have unscreened second neighbors, in which case we
      // add on the contribution from those (accounting for partial
      // screening)
      if (element_nn2_(a, a) == 1) {
        double distance_ratio;
        double second_neighbor_screening;

        double const z2 = NumSecondNearestNeighborsInReferenceStructure(element_lattice_(a, a),
                                                                        element_Cmin_(a, a, a),
                                                                        element_Cmax_(a, a, a),
                                                                        element_stheta_(a, a),
                                                                        distance_ratio,
                                                                        second_neighbor_screening);

        double const rho0_2nn = element_rho0_[a] * std::exp(-element_beta0_[a] * (distance_ratio - 1));

        rho0a += z2 * rho0_2nn * second_neighbor_screening;
      }

      // The composition dependent electron density rho0i in Eq.4.5
      element_ref_rho_[a] = rho0a * g_a;
    }
  }

  // Calculate density functions, assuming reference configuration
  void MEAMC::ComputeReferenceConfigurationDensity(double const r,
                                                   int const a,
                                                   int const b,
                                                   double *rho0_a,
                                                   double *rho1_a,
                                                   double *rho2_a,
                                                   double *rho3_a,
                                                   double *rho0_b,
                                                   double *rho1_b,
                                                   double *rho2_b,
                                                   double *rho3_b) {
    // In Eq.4.8
    double rnorm_a = r / element_re_(a, a) - 1.0;
    double rnorm_b = r / element_re_(b, b) - 1.0;

    // Element depdendent density scaling
    double const rhoa0 = element_rho0_[a];
    double const rhob0 = element_rho0_[b];

    // Eq.4.8
    double const rhoa0_a = rhoa0 * std::exp(-element_beta0_[a] * rnorm_a);
    double const rhoa1_a = rhoa0 * std::exp(-element_beta1_[a] * rnorm_a);
    double const rhoa2_a = rhoa0 * std::exp(-element_beta2_[a] * rnorm_a);
    double const rhoa3_a = rhoa0 * std::exp(-element_beta3_[a] * rnorm_a);

    // Eq.4.8
    double const rhoa0_b = rhob0 * std::exp(-element_beta0_[b] * rnorm_b);
    double const rhoa1_b = rhob0 * std::exp(-element_beta1_[b] * rnorm_b);
    double const rhoa2_b = rhob0 * std::exp(-element_beta2_[b] * rnorm_b);
    double const rhoa3_b = rhob0 * std::exp(-element_beta3_[b] * rnorm_b);

    Lattice const lat_ab = element_lattice_(a, b);

    double const z_ab = NumNearestNeighborsInReferenceStructure(lat_ab);

    *rho1_a = 0.0;
    *rho2_a = 0.0;
    *rho3_a = 0.0;

    *rho1_b = 0.0;
    *rho2_b = 0.0;
    *rho3_b = 0.0;

    switch (lat_ab) {
      case Lattice::FCC: {
        *rho0_a = 12.0 * rhoa0_b;
        *rho0_b = 12.0 * rhoa0_a;
        break;
      }
      case Lattice::BCC: {
        *rho0_a = 8.0 * rhoa0_b;
        *rho0_b = 8.0 * rhoa0_a;
        break;
      }
      case Lattice::B1: {
        *rho0_a = 6.0 * rhoa0_b;
        *rho0_b = 6.0 * rhoa0_a;
        break;
      }
      case Lattice::DIA:
      case Lattice::DIA3: {
        *rho0_a = 4.0 * rhoa0_b;
        *rho0_b = 4.0 * rhoa0_a;

        *rho3_a = 32.0 / 9.0 * rhoa3_b * rhoa3_b;
        *rho3_b = 32.0 / 9.0 * rhoa3_a * rhoa3_a;
        break;
      }
      case Lattice::HCP: {
        *rho0_a = 12 * rhoa0_b;
        *rho0_b = 12 * rhoa0_a;

        *rho3_a = kThird * rhoa3_b * rhoa3_b;
        *rho3_b = kThird * rhoa3_a * rhoa3_a;
        break;
      }
      case Lattice::DIM: {
        VectorOfSizeDIM shape_factors;

        GetShapeFactors(lat_ab, 0, 0, shape_factors);

        *rho0_a = rhoa0_b;
        *rho0_b = rhoa0_a;

        *rho1_a = shape_factors[0] * rhoa1_b * rhoa1_b;
        *rho1_b = shape_factors[0] * rhoa1_a * rhoa1_a;

        *rho2_a = shape_factors[1] * rhoa2_b * rhoa2_b;
        *rho2_b = shape_factors[1] * rhoa2_a * rhoa2_a;

        *rho3_a = shape_factors[2] * rhoa3_b * rhoa3_b;
        *rho3_b = shape_factors[2] * rhoa3_a * rhoa3_a;
        break;
      }
      case Lattice::C11: {
        *rho0_a = rhoa0_a;
        *rho0_b = rhoa0_b;

        *rho1_a = rhoa1_a;
        *rho1_b = rhoa1_b;

        *rho2_a = rhoa2_a;
        *rho2_b = rhoa2_b;

        *rho3_a = rhoa3_a;
        *rho3_b = rhoa3_b;
        break;
      }
      case Lattice::L12: {
        *rho0_a = 8 * rhoa0_a + 4 * rhoa0_b;
        *rho0_b = 12 * rhoa0_a;

        if (ialloy_ == 1) {
          *rho2_a = 8.0 / 3.0 * special::Square(rhoa2_a * element_t2_[a] - rhoa2_b * element_t2_[b]);
          double const denom = 8.0 * rhoa0_a * special::Square(element_t2_[a])
                               + 4.0 * rhoa0_b * special::Square(element_t2_[b]);
          if (denom > 0.) {
            *rho2_a = *rho2_a / denom * *rho0_a;
          }

        } else {
          *rho2_a = 8.0 / 3.0 * (rhoa2_a - rhoa2_b) * (rhoa2_a - rhoa2_b);
        }
        break;
      }
      case Lattice::B2: {
        *rho0_a = 8.0 * rhoa0_b;
        *rho0_b = 8.0 * rhoa0_a;
        break;
      }
      case Lattice::CH4: {
        *rho0_a = 4.0 * rhoa0_b;  // in assumption that 'a' represent carbon
        *rho0_b = rhoa0_a;        // in assumption that 'b' represent hydrogen

        VectorOfSizeDIM shape_factors;

        GetShapeFactors(Lattice::DIM, 0, 0, shape_factors);  // H
        *rho1_b = shape_factors[0] * rhoa1_a * rhoa1_a;
        *rho2_b = shape_factors[1] * rhoa2_a * rhoa2_a;
        *rho3_b = shape_factors[2] * rhoa3_a * rhoa3_a;

        GetShapeFactors(lat_ab, 0, 0, shape_factors);  // C
        *rho1_a = shape_factors[0] * rhoa1_b * rhoa1_b;
        *rho2_a = shape_factors[1] * rhoa2_b * rhoa2_b;
        *rho3_a = shape_factors[2] * rhoa3_b * rhoa3_b;
        break;
      }
      case Lattice::LIN: {
        *rho0_a = rhoa0_b * z_ab;
        *rho0_b = rhoa0_a * z_ab;

        VectorOfSizeDIM shape_factors;

        GetShapeFactors(lat_ab, element_stheta_(a, b), element_ctheta_(a, b), shape_factors);

        *rho1_b = shape_factors[0] * rhoa1_a * rhoa1_a;
        *rho2_b = shape_factors[1] * rhoa2_a * rhoa2_a;
        *rho3_b = shape_factors[2] * rhoa3_a * rhoa3_a;

        *rho1_a = shape_factors[0] * rhoa1_b * rhoa1_b;
        *rho2_a = shape_factors[1] * rhoa2_b * rhoa2_b;
        *rho3_a = shape_factors[2] * rhoa3_b * rhoa3_b;
        break;
      }
      case Lattice::ZIG: {
        *rho0_a = rhoa0_b * z_ab;
        *rho0_b = rhoa0_a * z_ab;

        VectorOfSizeDIM shape_factors;

        GetShapeFactors(lat_ab, element_stheta_(a, b), element_ctheta_(a, b), shape_factors);

        *rho1_b = shape_factors[0] * rhoa1_a * rhoa1_a;
        *rho2_b = shape_factors[1] * rhoa2_a * rhoa2_a;
        *rho3_b = shape_factors[2] * rhoa3_a * rhoa3_a;

        *rho1_a = shape_factors[0] * rhoa1_b * rhoa1_b;
        *rho2_a = shape_factors[1] * rhoa2_b * rhoa2_b;
        *rho3_a = shape_factors[2] * rhoa3_b * rhoa3_b;
        break;
      }
      case Lattice::TRI: {
        *rho0_a = rhoa0_b;
        *rho0_b = rhoa0_a * z_ab;

        VectorOfSizeDIM shape_factors;

        GetShapeFactors(lat_ab, element_stheta_(a, b), element_ctheta_(a, b), shape_factors);

        *rho1_b = shape_factors[0] * rhoa1_a * rhoa1_a;
        *rho2_b = shape_factors[1] * rhoa2_a * rhoa2_a;
        *rho3_b = shape_factors[2] * rhoa3_a * rhoa3_a;

        shape_factors[0] = 1.0;
        shape_factors[1] = kTwoThird;
        shape_factors[2] = 1.0 - 0.6 * shape_factors[0];

        *rho1_a = shape_factors[0] * rhoa1_b * rhoa1_b;
        *rho2_a = shape_factors[1] * rhoa2_b * rhoa2_b;
        *rho3_a = shape_factors[2] * rhoa3_b * rhoa3_b;
        break;
      }
    }

    if (element_nn2_(a, b) == 1) {
      double distance_ratio;
      double second_neighbor_screening;

      double z_ab_nn2 = NumSecondNearestNeighborsInReferenceStructure(lat_ab,
                                                                      element_Cmin_(a, a, b),
                                                                      element_Cmax_(a, a, b),
                                                                      element_stheta_(a, b),
                                                                      distance_ratio,
                                                                      second_neighbor_screening);

      rnorm_a = distance_ratio * r / element_re_(a, a) - 1.0;
      rnorm_b = distance_ratio * r / element_re_(b, b) - 1.0;

      double const rhoa0_1_nn2 = element_rho0_[a] * std::exp(-element_beta0_[a] * rnorm_a);
      double const rhoa0_2_nn2 = element_rho0_[b] * std::exp(-element_beta0_[b] * rnorm_b);

      if (lat_ab == Lattice::L12) {
        // As usual, L12 thinks it'shape_factors special; we need
        // to be careful computing the screening functions
        double const s_aaa = Sijk(1.0, a, a, a);
        double const s_aab = Sijk(1.0, a, a, b);
        double const s_bba = Sijk(1.0, b, b, a);

        double const S11 = s_aaa * s_aaa * s_aab * s_aab;
        double const S22 = special::PowInt(s_bba, 4);

        *rho0_a += 6.0 * S11 * rhoa0_1_nn2;
        *rho0_b += 6.0 * S22 * rhoa0_2_nn2;

      } else {
        // For other cases, assume that second neighbor is of the
        // same type, first neighbor may be of different type
        *rho0_a += z_ab_nn2 * second_neighbor_screening * rhoa0_1_nn2;

        // Assume z_ab_nn2 and distance_ratio don't depend on
        // order, but second_neighbor_screening might
        z_ab_nn2 = NumSecondNearestNeighborsInReferenceStructure(lat_ab,
                                                                 element_Cmin_(b, b, a),
                                                                 element_Cmax_(b, b, a),
                                                                 element_stheta_(a, b),
                                                                 distance_ratio,
                                                                 second_neighbor_screening);

        *rho0_b += z_ab_nn2 * second_neighbor_screening * rhoa0_2_nn2;
      }
    }
  }

  void MEAMC::SplineInterpolate(int const ind) {
    double const *const phirar_1d = &phir_(ind, 0);
    double *const phirar1_1d = &phirar1_(ind, 0);

    phirar1_1d[0] = phirar_1d[1] - phirar_1d[0];
    phirar1_1d[1] = 0.5 * (phirar_1d[2] - phirar_1d[0]);
    phirar1_1d[nr_ - 2] = 0.5 * (phirar_1d[nr_ - 1] - phirar_1d[nr_ - 3]);
    phirar1_1d[nr_ - 1] = 0.0;

    for (int j = 2; j < nr_ - 2; ++j) {
      phirar1_1d[j]
          = ((phirar_1d[j - 2] - phirar_1d[j + 2]) + 8.0 * (phirar_1d[j + 1] - phirar_1d[j - 1])) / 12.;
    }

    double *const phirar2_1d = &phirar2_(ind, 0);

    for (int j = 0; j < nr_ - 1; ++j) {
      phirar2_1d[j] = 3.0 * (phirar_1d[j + 1] - phirar_1d[j]) - 2.0 * phirar1_1d[j] - phirar1_1d[j + 1];
    }
    phirar2_1d[nr_ - 1] = 0.0;

    double *const phirar3_1d = &phirar3_(ind, 0);

    for (int j = 0; j < nr_ - 1; ++j) {
      phirar3_1d[j] = phirar1_1d[j] + phirar1_1d[j + 1] - 2.0 * (phirar_1d[j + 1] - phirar_1d[j]);
    }
    phirar3_1d[nr_ - 1] = 0.0;

    double const rdrar = 1.0 / dr_;

    double *const phirar4_1d = &phirar4_(ind, 0);

    for (int j = 0; j < nr_; ++j) {
      phirar4_1d[j] = phirar1_1d[j] * rdrar;
    }

    double *const phirar5_1d = &phirar5_(ind, 0);

    for (int j = 0; j < nr_; ++j) {
      phirar5_1d[j] = 2.0 * phirar2_1d[j] * rdrar;
    }

    double *const phirar6_1d = &phirar6_(ind, 0);

    for (int j = 0; j < nr_; ++j) {
      phirar6_1d[j] = 3.0 * phirar3_1d[j] * rdrar;
    }
  }

}  // namespace openKIM