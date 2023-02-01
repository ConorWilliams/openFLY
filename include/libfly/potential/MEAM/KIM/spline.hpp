//
// spline.hpp
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

#ifndef SPLINE_HPP
#define SPLINE_HPP

#pragma GCC system_header

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "utility.hpp"

#ifdef SPLINE_LOG_ERROR
#undef SPLINE_LOG_ERROR
#endif

namespace openKIM {

/*!
 * \brief Helper macro for printing error message
 * formats messages, filename, line number and function
 * name into an std::ostringstream object
 *
 */
#define SPLINE_LOG_ERROR(msg)                                                 \
  {                                                                           \
    std::ostringstream ss;                                                    \
    ss << "\nError :" << __FILE__ << ":" << __LINE__ << ":@(" << __FUNCTION__ \
       << ")\n"                                                               \
       << msg << "\n\n";                                                      \
    std::cerr << ss.str();                                                    \
  }

/*! \class Spline
 * \brief Spline-based Modified Embedded Atom method (MEAM) potential.
 *
 */
class Spline {
 public:
  /*!
   * \brief Get the Number Of Knots
   *
   * \return int
   */
  inline int GetNumberOfKnots() const;

  /*!
   * \brief Process and parse the spline knots from a text file.
   *
   * \param parameter_file_pointer FILE pointer to the opened file
   * \param max_line_size maximum line size
   * \param is_new_format if the potential file in a new format or not
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  inline int ProcessParameterFile(std::FILE *const parameter_file_pointer,
                                  int const max_line_size,
                                  bool const is_new_format = false);

  /*!
   * \brief Calculates the second derivatives of the cubic spline.
   *
   */
  inline void PrepareSpline();

  /*!
   * \brief Calculates the second derivatives of the cubic spline after unit
   *        conversion.
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  inline int UpdateSpline();

  /*!
   * \brief Return true if knots are on a uniform grid
   *
   * \return int
   */
  inline int AreKnotsOnRegularGrid() const;

  /*!
   * \brief Returns the cutoff radius of this function.
   *
   * \return double
   */
  inline double GetCutoff() const;

  /*!
   * \brief Evaluates the spline function at position x.
   *
   * \tparam KnotsAreOnRegularGid If true, use only MEAM potentials with spline
   *         knots on a uniform grid. If not use MEAM potentials with
   *         non-uniform spline knots.
   * \param x Position x
   * \return double
   */
  template <bool KnotsAreOnRegularGid>
  inline double Eval(double const x) const;

  /*!
   * \brief Evaluates the spline function and its first derivative at position
   * x.
   *
   * \tparam KnotsAreOnRegularGid If true, use only MEAM potentials with spline
   *         knots on a uniform grid. If not use MEAM potentials with
   *         non-uniform spline knots.
   * \param x Position x
   * \param deriv First derivative at position x.
   * \return double
   */
  template <bool KnotsAreOnRegularGid>
  inline double Eval(double const x, double &deriv) const;

  /*!
   * \brief Convert units of the parameters
   *
   * \param convert_x_factor x unit conversion factor
   * \param convert_value_factor function value unit conversion factor
   */
  inline void ConvertUnit(double const convert_x_factor,
                          double const convert_value_factor);

  /*!
   * \brief Writes a Gnuplot script that plots the spline function.
   *
   * \param filename
   * \param title
   *
   * \note
   * This function is provided for debugging
   */
  inline void WriteGnuplot(const char *filename,
                           const char *title = nullptr) const;

  /*!
   * \brief Returns pointer to the underlying first derivative at the beginning
   *        of spline
   *
   * \return double*
   */
  inline double const *derivative_0_data() const noexcept;
  inline double *derivative_0_data() noexcept;

  /*!
   * \brief Returns pointer to the underlying first derivative at the end of
   *        spline
   *
   * \return double*
   */
  inline double const *derivative_n_data() const noexcept;
  inline double *derivative_n_data() noexcept;

  /*!
   * \brief Returns pointer to the underlying knots x array.
   *
   * \return double*
   */
  inline double const *x_data() const noexcept;
  inline double *x_data() noexcept;

  /*!
   * \brief Returns pointer to the underlying knots value array
   *
   * \return double*
   */
  inline double const *value_data() const noexcept;
  inline double *value_data() noexcept;

  /*!
   * \brief Returns pointer to the underlying second derivative value
   *        at spline knots
   *
   * \return double const*
   */
  inline double const *second_derivative_value_data() const noexcept;
  inline double *second_derivative_value_data() noexcept;

  /*!
   * \brief Return the old format string
   *
   * \return std::string
   */
  inline std::string OldFormatString() const noexcept;

 private:
  /*! Number of spline knots */
  int number_knots_{0};

  /*! Indicates that all spline knots are on a regular grid. */
  int knots_are_on_regular_grid_;

  /*! The beginning of the interval on which the spline function is defined. */
  double x_min_;

  /*! The end of the interval on which the spline function is defined. */
  double x_max_;

  /*! The end of the spline interval after it has been shifted to begin at
   * knots_x_=0. */
  double x_max_shifted_;

  /*! The distance between knots if this is a grid spline with equidistant
   *  knots. */
  double regular_grid_distance_;

  /*! The squared distance between knots if this is a grid spline with
   *  equidistant knots. */
  double regular_grid_distance_square_;

  /*! 1/regular_grid_distance_ */
  double regular_grid_distance_inverse_;

  /*! First derivative of the in terpolating function at knot 0 */
  double knot_1_first_derivative_;

  /*! First derivative of the in terpolating function at knot
   *  (number_knots_-1) */
  double knot_n_first_derivative_;

  /*! String in the old formatted potential file */
  std::string old_format_;

  /*! Positions of spline knots */
  std::vector<double> knots_x_;

  /*! Shifted positions of spline knots */
  std::vector<double> knots_shifted_x_;

  /*! Function values at spline knots */
  std::vector<double> knots_value_;

  /*! Second derivatives of the interpolating function at the spline knots */
  std::vector<double> knots_second_derivative_value_;

  /*! knots_value_delta_[i] =
   *  (knots_value_[i+1]-knots_value_[i])/regular_grid_distance_ */
  std::vector<double> knots_value_delta_;
};

inline int Spline::GetNumberOfKnots() const { return number_knots_; }

inline int Spline::ProcessParameterFile(std::FILE *const parameter_file_pointer,
                                        int const max_line_size,
                                        bool const is_new_format) {
  std::unique_ptr<char> next_line_(new char[max_line_size]);
  char *next_line = next_line_.get();

  // Create a utility object
  Utility ut;

  // If new format, read the spline format.
  if (is_new_format) {
    // should always be "spline3eq" for now.
    if (ut.GetNextLine(parameter_file_pointer, next_line, max_line_size)) {
      std::string msg = "End of file while reading the ";
      msg += "spline knots from the potential file.\n";
      SPLINE_LOG_ERROR(msg);
      return true;
    }
  }

  // Read the next line
  if (ut.GetNextLine(parameter_file_pointer, next_line, max_line_size)) {
    std::string msg = "End of file while reading the ";
    msg += "spline knots from the potential file.\n";
    SPLINE_LOG_ERROR(msg);
    return true;
  }

  // The number of spline knots.
  number_knots_ = std::atoi(next_line);
  if (number_knots_ < 2) {
    std::string msg = "Invalid number of spline knots (n = ";
    msg += std::to_string(number_knots_) + ") in the potential file.\n";
    SPLINE_LOG_ERROR(msg);
    return true;
  }

  // Set the number of elements
  knots_x_.resize(number_knots_);
  knots_shifted_x_.resize(number_knots_);
  knots_value_.resize(number_knots_);
  knots_second_derivative_value_.resize(number_knots_);
  knots_value_delta_.resize(number_knots_);

  // Read the next line
  if (ut.GetNextLine(parameter_file_pointer, next_line, max_line_size)) {
    std::string msg = "End of file while reading the ";
    msg += "spline knots from the potential file.\n";
    SPLINE_LOG_ERROR(msg);
    return true;
  }

  // Parse first derivatives at the beginning and end of spline.
  knot_1_first_derivative_ = std::atof(std::strtok(next_line, " \t\n\r\f"));
  knot_n_first_derivative_ = std::atof(std::strtok(nullptr, " \t\n\r\f"));

  // Read the next line and skip it in the old format
  if (!is_new_format) {
    // Skip a line in the old format
    if (ut.GetNextLine(parameter_file_pointer, next_line, max_line_size)) {
      std::string msg = "End of file while reading the ";
      msg += "spline knots from the potential file.\n";
      SPLINE_LOG_ERROR(msg);
      return true;
    }
    old_format_ = std::string(next_line);
    old_format_ =
        old_format_.substr(old_format_.find_first_not_of(" \t\n\r\f"));
    old_format_ =
        old_format_.substr(0, old_format_.find_last_not_of(" \t\n\r\f") + 1);
  }

  // Parse knot coordinates.
  for (int i = 0; i < number_knots_; ++i) {
    // Read the next line
    if (ut.GetNextLine(parameter_file_pointer, next_line, max_line_size)) {
      std::string msg = "End of file while reading the ";
      msg += "spline knots from the potential file.\n";
      SPLINE_LOG_ERROR(msg);
      return true;
    }

    double knot_x;
    double knot_value;
    double knot_second_derivative_value;

    if (std::sscanf(next_line, "%lg %lg %lg", &knot_x, &knot_value,
                    &knot_second_derivative_value) != 3) {
      std::string msg = "Invalid line with the wrong number ";
      msg += "of input knots in the potential file.\n";
      SPLINE_LOG_ERROR(msg);
      return true;
    }

    knots_x_[i] = knot_x;
    knots_value_[i] = knot_value;
  }

  PrepareSpline();

  // Everything is good
  return false;
}

inline void Spline::PrepareSpline() {
  x_min_ = knots_x_[0];
  x_max_ = knots_x_[number_knots_ - 1];
  x_max_shifted_ = x_max_ - x_min_;

  regular_grid_distance_ = x_max_shifted_ / (number_knots_ - 1);
  regular_grid_distance_square_ =
      regular_grid_distance_ * regular_grid_distance_;
  regular_grid_distance_inverse_ = 1.0 / regular_grid_distance_;

  knots_are_on_regular_grid_ = 1;
  for (int i = 1; i < number_knots_ - 1; ++i) {
    if (std::abs(regular_grid_distance_ * i + x_min_ - knots_x_[i]) > 1e-8) {
      knots_are_on_regular_grid_ = 0;
      break;
    }
  }
}

inline int Spline::UpdateSpline() {
  x_min_ = knots_x_[0];
  x_max_ = knots_x_[number_knots_ - 1];
  x_max_shifted_ = x_max_ - x_min_;

  regular_grid_distance_ = x_max_shifted_ / (number_knots_ - 1);
  regular_grid_distance_square_ =
      regular_grid_distance_ * regular_grid_distance_;
  regular_grid_distance_inverse_ = 1.0 / regular_grid_distance_;

  int knots_are_on_regular_grid = 1;
  for (int i = 1; i < number_knots_ - 1; ++i) {
    if (std::abs(regular_grid_distance_ * i + x_min_ - knots_x_[i]) > 1e-8) {
      knots_are_on_regular_grid = 0;
      break;
    }
  }

  if (knots_are_on_regular_grid != knots_are_on_regular_grid_) {
    std::string msg = "The knots distribution form has changed.\n";
    msg += "The knots are no longer ";
    msg += knots_are_on_regular_grid_ ? "on uniform grid.\n" : "non-uniform.\n";
    SPLINE_LOG_ERROR(msg);
    return true;
  }

  std::vector<double> u(number_knots_);

  {
    double const dx1 = knots_x_[1] - knots_x_[0];
    double const dv1 = knots_value_[1] - knots_value_[0];

    u[0] = (3.0 / dx1) * (dv1 / dx1 - knot_1_first_derivative_);
  }

  knots_second_derivative_value_[0] = -0.5;

  // The decomposition loop of the tridiagonal algorithm.
  // y2 and u are used for temporary storage of the decomposed factors.
  for (int i = 1; i < number_knots_ - 1; ++i) {
    double const dx1 = knots_x_[i] - knots_x_[i - 1];
    double const dx2 = knots_x_[i + 1] - knots_x_[i - 1];
    double const dx3 = knots_x_[i + 1] - knots_x_[i];
    double const sig = dx1 / dx2;
    double const p = sig * knots_second_derivative_value_[i - 1] + 2.0;
    knots_second_derivative_value_[i] = (sig - 1.0) / p;
    double const dv1 = knots_value_[i + 1] - knots_value_[i];
    double const dv2 = knots_value_[i] - knots_value_[i - 1];
    double const ui = dv1 / dx3 - dv2 / dx1;
    u[i] = (6.0 * ui / dx2 - sig * u[i - 1]) / p;
  }

  {
    double const dx1_inv =
        1.0 / (knots_x_[number_knots_ - 1] - knots_x_[number_knots_ - 2]);
    double const dv1 =
        knots_value_[number_knots_ - 1] - knots_value_[number_knots_ - 2];

    double const qn = 0.5;
    double const un =
        3.0 * dx1_inv * (knot_n_first_derivative_ - dv1 * dx1_inv);

    knots_second_derivative_value_[number_knots_ - 1] =
        (un - qn * u[number_knots_ - 2]) /
        (qn * knots_second_derivative_value_[number_knots_ - 2] + 1.0);
  }

  // This is the back substitution loop of the tridiagonal algorithm.
  for (int i = number_knots_ - 2; i > -1; --i) {
    knots_second_derivative_value_[i] =
        knots_second_derivative_value_[i] *
            knots_second_derivative_value_[i + 1] +
        u[i];
  }

  // Shift the spline to knots_x_=0 to speed up interpolation.
  for (int i = 0; i < number_knots_; ++i) {
    knots_shifted_x_[i] = knots_x_[i] - x_min_;
  }

  if (knots_are_on_regular_grid_) {
    for (int i = 0; i < number_knots_ - 1; ++i) {
      double const dv1 = knots_value_[i + 1] - knots_value_[i];
      knots_value_delta_[i] = dv1 / regular_grid_distance_;
    }

    double const c = regular_grid_distance_ * 6.0;

    for (int i = 0; i < number_knots_; ++i) {
      knots_second_derivative_value_[i] /= c;
    }
  }

  // Everything is good
  return false;
}

inline int Spline::AreKnotsOnRegularGrid() const {
  return knots_are_on_regular_grid_;
}

inline double Spline::GetCutoff() const { return knots_x_[number_knots_ - 1]; }

template <bool KnotsAreOnRegularGid>
inline double Spline::Eval(double const x) const {
  double const x_shifted = x - x_min_;
  // Left extrapolation.
  if (x_shifted <= 0.0) {
    double val = knots_value_[0] + knot_1_first_derivative_ * x_shifted;
    return val;
  }

  // Right extrapolation.
  if (x_shifted >= x_max_shifted_) {
    double val = knots_value_[number_knots_ - 1] +
                 knot_n_first_derivative_ * (x_shifted - x_max_shifted_);
    return val;
  }

  if (KnotsAreOnRegularGid) {
    // Calculate the interval knots_x_ is in.
    int lo = static_cast<int>(x_shifted * regular_grid_distance_inverse_);
    int const hi = lo + 1;
    // Do spline interpolation.
    double const a = knots_shifted_x_[hi] - x_shifted;
    double const b = regular_grid_distance_ - a;
    double const aa = (a * a - regular_grid_distance_square_) * a;
    double const bb = (b * b - regular_grid_distance_square_) * b;
    double val = knots_value_[hi] - a * knots_value_delta_[lo] +
                 (aa * knots_second_derivative_value_[lo] +
                  bb * knots_second_derivative_value_[hi]);
    return val;
  } else {
    // knots are not on regular grid
    // Do interval search.
    int lo = 0;
    int hi = number_knots_ - 1;
    while (hi - lo > 1) {
      int const k = (hi + lo) / 2;
      if (knots_shifted_x_[k] > x_shifted) {
        hi = k;
      } else {
        lo = k;
      }
    }
    double const grid_distance = knots_shifted_x_[hi] - knots_shifted_x_[lo];
    double const grid_distance_square = grid_distance * grid_distance;
    // Do spline interpolation.
    double const a = (knots_shifted_x_[hi] - x_shifted) / grid_distance;
    double const b = 1.0 - a;
    double const aa = (a * a - 1.0) * a * knots_second_derivative_value_[lo];
    double const bb = (b * b - 1.0) * b * knots_second_derivative_value_[hi];
    double val = a * knots_value_[lo] + b * knots_value_[hi] +
                 (aa + bb) * grid_distance_square / 6.0;
    return val;
  }
}

template <bool KnotsAreOnRegularGid>
inline double Spline::Eval(double const x, double &deriv) const {
  double const x_shifted = x - x_min_;
  // Left extrapolation.
  if (x_shifted <= 0.0) {
    deriv = knot_1_first_derivative_;
    double val = knots_value_[0] + knot_1_first_derivative_ * x_shifted;
    return val;
  }

  // Right extrapolation.
  if (x_shifted >= x_max_shifted_) {
    deriv = knot_n_first_derivative_;
    double val = knots_value_[number_knots_ - 1] +
                 knot_n_first_derivative_ * (x_shifted - x_max_shifted_);
    return val;
  }

  if (KnotsAreOnRegularGid) {
    // Calculate the interval knots_x_ is in.
    int lo = static_cast<int>(x_shifted * regular_grid_distance_inverse_);
    int const hi = lo + 1;
    // Do spline interpolation.
    double const a = knots_shifted_x_[hi] - x_shifted;
    double const b = regular_grid_distance_ - a;
    double const aad = 3.0 * a * a - regular_grid_distance_square_;
    double const bbd = 3.0 * b * b - regular_grid_distance_square_;
    deriv = knots_value_delta_[lo] + (bbd * knots_second_derivative_value_[hi] -
                                      aad * knots_second_derivative_value_[lo]);
    double const aa = (a * a - regular_grid_distance_square_) * a;
    double const bb = (b * b - regular_grid_distance_square_) * b;
    double val = knots_value_[hi] - a * knots_value_delta_[lo] +
                 (aa * knots_second_derivative_value_[lo] +
                  bb * knots_second_derivative_value_[hi]);
    return val;
  } else {
    // knots are not on regular grid
    // Do interval search.
    int lo = 0;
    int hi = number_knots_ - 1;
    while (hi - lo > 1) {
      int const k = (hi + lo) / 2;
      if (knots_shifted_x_[k] > x_shifted) {
        hi = k;
      } else {
        lo = k;
      }
    }
    double const grid_distance = knots_shifted_x_[hi] - knots_shifted_x_[lo];
    double const grid_distance_square = grid_distance * grid_distance;
    // Do spline interpolation.
    double const a = (knots_shifted_x_[hi] - x_shifted) / grid_distance;
    double const b = 1.0 - a;
    double const c = (knots_value_[hi] - knots_value_[lo]) / grid_distance;
    double const aad = (3.0 * a * a - 1.0) * knots_second_derivative_value_[hi];
    double const bbd = (3.0 * b * b - 1.0) * knots_second_derivative_value_[lo];
    deriv = c + (aad - bbd) * grid_distance / 6.0;
    double const aa = (a * a - 1.0) * a * knots_second_derivative_value_[lo];
    double const bb = (b * b - 1.0) * b * knots_second_derivative_value_[hi];
    double val = a * knots_value_[lo] + b * knots_value_[hi] +
                 (aa + bb) * grid_distance_square / 6.0;
    return val;
  }
}

inline void Spline::ConvertUnit(double const convert_x_factor,
                                double const convert_value_factor) {
  if (special::IsOne(convert_x_factor) &&
      special::IsOne(convert_value_factor)) {
    return;
  }

  if (special::IsNotOne(convert_x_factor)) {
    std::for_each(knots_x_.begin(), knots_x_.end(),
                  [&](double &x) { x *= convert_x_factor; });
  }

  if (special::IsNotOne(convert_value_factor)) {
    std::for_each(knots_value_.begin(), knots_value_.end(),
                  [&](double &v) { v *= convert_value_factor; });
  }

  double const convert_factor =
      special::FloatDivZero(convert_value_factor, convert_x_factor);
  knot_1_first_derivative_ *= convert_factor;
  knot_n_first_derivative_ *= convert_factor;
}

inline void Spline::WriteGnuplot(const char *filename,
                                 const char *title) const {
  std::FILE *fp = std::fopen(filename, "w");

  std::fprintf(fp, "#!/usr/bin/env gnuplot\n");
  if (title) {
    std::fprintf(fp, "set title \"%s\"\n", title);
  }

  double const tmin =
      knots_x_[0] - (knots_x_[number_knots_ - 1] - knots_x_[0]) * 0.05;
  double const tmax = knots_x_[number_knots_ - 1] +
                      (knots_x_[number_knots_ - 1] - knots_x_[0]) * 0.05;

  double const delta = (tmax - tmin) / (number_knots_ * 200);

  std::fprintf(fp, "set xrange [%f:%f]\n", tmin, tmax);
  std::fprintf(fp,
               "plot '-' with lines notitle, '-' with points "
               "notitle pt 3 lc 3\n");

  if (knots_are_on_regular_grid_) {
    for (double x = tmin; x <= tmax + 1e-8; x += delta) {
      double const y = Eval<true>(x);
      std::fprintf(fp, "%f %f\n", x, y);
    }
  } else {
    for (double x = tmin; x <= tmax + 1e-8; x += delta) {
      double const y = Eval<false>(x);
      std::fprintf(fp, "%f %f\n", x, y);
    }
  }

  std::fprintf(fp, "e\n");
  for (int i = 0; i < number_knots_; i++) {
    std::fprintf(fp, "%f %f\n", knots_x_[i], knots_value_[i]);
  }

  std::fprintf(fp, "e\n");
  std::fclose(fp);
}

inline double const *Spline::derivative_0_data() const noexcept {
  return &knot_1_first_derivative_;
}

inline double *Spline::derivative_0_data() noexcept {
  return &knot_1_first_derivative_;
}

inline double const *Spline::derivative_n_data() const noexcept {
  return &knot_n_first_derivative_;
}

inline double *Spline::derivative_n_data() noexcept {
  return &knot_n_first_derivative_;
}

inline double const *Spline::x_data() const noexcept { return knots_x_.data(); }

inline double *Spline::x_data() noexcept { return knots_x_.data(); }

inline double const *Spline::value_data() const noexcept {
  return knots_value_.data();
}

inline double *Spline::value_data() noexcept { return knots_value_.data(); }

inline double const *Spline::second_derivative_value_data() const noexcept {
  return knots_second_derivative_value_.data();
}

inline double *Spline::second_derivative_value_data() noexcept {
  return knots_second_derivative_value_.data();
}

inline std::string Spline::OldFormatString() const noexcept {
  return old_format_;
}

}

#undef SPLINE_LOG_ERROR
#endif  // SPLINE_HPP
