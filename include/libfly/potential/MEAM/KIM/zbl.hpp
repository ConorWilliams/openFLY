//
// zbl.hpp
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

#ifndef ZBL_HPP
#define ZBL_HPP

#pragma GCC system_header

#include <cmath>
#include <cstddef>

#include "helper.hpp"
#include "special.hpp"

namespace openKIM {

  /*!
   * \brief ZBL pair interaction style
   *
   * It computes the Ziegler-Biersack-Littmark (ZBL) screened nuclear repulsion
   * for describing high-energy collisions between atoms.
   *
   * The ZBL interaction is already smoothed to 0.0 at the cutoff.
   */
  class ZBL {
  public:
    /*!
     * \brief Construct a new ZBL object
     *
     * \param number_of_element_types The number of element types
     */
    ZBL(std::size_t const number_of_element_types);

    /*!
     * \brief Set the coeff
     *
     * \param species_i i atom type
     * \param species_j j atom type
     * \param z_species_i The nuclear charge of the species_i atoms
     * \param z_species_j The nuclear charge of the species_j atoms
     */
    inline void SetCoeff(int const species_i,
                         int const species_j,
                         double const z_species_i,
                         double const z_species_j);

    /*!
     * \brief Convert units of the parameters
     *
     * \param convert_length_factor length unit conversion factor
     * \param convert_energy_factor energy unit conversion factor
     */
    inline void ConvertUnit(double const convert_length_factor, double const convert_energy_factor);

    /*!
     * \brief Compute ZBL pair energy
     *
     * \param rij Pair atoms distance
     * \param species_i i atom type
     * \param species_j j atom type
     * \return double
     */
    inline double ComputePhi(double const rij, int const species_i, int const species_j);

  private:
    /*! Conversion of q^2/r to energy */
    double qqr2e{14.399645354084361};

    /*! ZBL constant with length unit*/
    double a0{0.46850};

    Array2D<double> d1a;
    Array2D<double> d2a;
    Array2D<double> d3a;
    Array2D<double> d4a;
    Array2D<double> zze;
  };

  inline ZBL::ZBL(std::size_t const number_of_element_types) {
    d1a.resize(number_of_element_types, number_of_element_types, 0.0);
    d2a.resize(number_of_element_types, number_of_element_types, 0.0);
    d3a.resize(number_of_element_types, number_of_element_types, 0.0);
    d4a.resize(number_of_element_types, number_of_element_types, 0.0);
    zze.resize(number_of_element_types, number_of_element_types, 0.0);
  }

  inline void ZBL::ConvertUnit(double const convert_length_factor, double const convert_energy_factor) {
    if (special::IsNotOne(convert_length_factor)) {
      qqr2e *= convert_length_factor;
      a0 *= convert_length_factor;
    }
    if (special::IsNotOne(convert_energy_factor)) {
      qqr2e *= convert_energy_factor;
    }
  }

  inline void ZBL::SetCoeff(int const species_i,
                            int const species_j,
                            double const z_species_i,
                            double const z_species_j) {
    double const ainv = (std::pow(z_species_i, 0.230) + std::pow(z_species_j, 0.230)) / a0;
    d1a(species_i, species_j) = 0.201620 * ainv;
    d2a(species_i, species_j) = 0.402900 * ainv;
    d3a(species_i, species_j) = 0.942290 * ainv;
    d4a(species_i, species_j) = 3.199800 * ainv;
    zze(species_i, species_j) = z_species_i * z_species_j * qqr2e;
    if (species_i != species_j) {
      d1a(species_j, species_i) = d1a(species_i, species_j);
      d2a(species_j, species_i) = d2a(species_i, species_j);
      d3a(species_j, species_i) = d3a(species_i, species_j);
      d4a(species_j, species_i) = d4a(species_i, species_j);
      zze(species_j, species_i) = zze(species_i, species_j);
    }
  }

  inline double ZBL::ComputePhi(double const rij, int const species_i, int const species_j) {
    double const d1aij = d1a(species_i, species_j);
    double const d2aij = d2a(species_i, species_j);
    double const d3aij = d3a(species_i, species_j);
    double const d4aij = d4a(species_i, species_j);
    double const zzeij = zze(species_i, species_j);
    double const sum1 = 0.028170 * std::exp(-d1aij * rij);
    double const sum2 = 0.280220 * std::exp(-d2aij * rij);
    double const sum3 = 0.509860 * std::exp(-d3aij * rij);
    double const sum4 = 0.181750 * std::exp(-d4aij * rij);
    double const sum = sum1 + sum2 + sum3 + sum4;
    if (rij > 0.0) {
      double result = zzeij * sum / rij;
      return result;
    }
    return sum;
  }

}  // namespace openKIM

#endif  // ZBL_HPP