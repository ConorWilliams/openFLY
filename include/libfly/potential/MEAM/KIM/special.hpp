//
// special.hpp
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

#ifndef SPECIAL_HPP
#define SPECIAL_HPP

#pragma GCC system_header

#include <cmath>
#include <type_traits>

/*!
 * \file special.hpp
 *
 * \brief This file contains special helper functions and a Pi constant
 *
 * The special helper functions and a constant are:
 * \c MY_PI Pi constant
 * \c Square x**2, use instead of std::pow(x, 2.0)
 * \c Cube x**3, use instead of std::pow(x, 3.0)
 * \c PowInt optimized version of std::pow(x, n) with n being integer
 * \c IsOne return true if the input is (close to) one within a tolerance
 * \c IsNotOne return true if the input is not (close to) one within a tolerance
 * \c IsZero return true if the input is (close to) zero within a tolerance
 * \c IsNotZero return true if the input is not (close to) zero within a
 *              tolerance
 * \c FloatDivZero returns zero if the denominator (d) is zero in f/d
 */

namespace openKIM {

  namespace special {
    static constexpr double MY_PI = 3.14159265358979323846;

    template <typename DataType>
    static inline constexpr typename std::enable_if<std::is_floating_point<DataType>::value, DataType>::type
    Square(DataType const &x) {
      return x * x;
    }

    template <typename DataType>
    static inline constexpr typename std::enable_if<std::is_floating_point<DataType>::value, DataType>::type
    Cube(DataType const &x) {
      return x * x * x;
    }

    template <typename DataType>
    static inline typename std::enable_if<std::is_floating_point<DataType>::value, DataType>::type PowInt(
        DataType const &x,
        int const n) {
      if (x == 0.0) {
        return 0.0;
      }
      int nn = (n > 0) ? n : -n;
      DataType ww = x;
      DataType yy{1};
      for (; nn != 0; nn >>= 1, ww *= ww) {
        if (nn & 1) {
          yy *= ww;
        }
      }
      return (n > 0) ? yy : 1.0 / yy;
    }

    template <typename DataType>
    static inline constexpr typename std::enable_if<std::is_floating_point<DataType>::value, bool>::type
    IsOne(DataType const f, double const tol = 1e-20) {
      return std::abs(f - 1.0) < tol;
    }

    template <typename DataType>
    static inline constexpr typename std::enable_if<std::is_floating_point<DataType>::value, bool>::type
    IsNotOne(DataType const f, double const tol = 1e-20) {
      return std::abs(f - 1.0) >= tol;
    }

    template <typename DataType>
    static inline constexpr typename std::enable_if<std::is_floating_point<DataType>::value, bool>::type
    IsZero(DataType const f, double const tol = 1e-20) {
      return std::abs(f) < tol;
    }

    template <typename DataType>
    static inline constexpr typename std::enable_if<std::is_floating_point<DataType>::value, bool>::type
    IsNotZero(DataType const f, double const tol = 1e-20) {
      return std::abs(f) >= tol;
    }

    template <typename DataType>
    static inline constexpr typename std::enable_if<std::is_floating_point<DataType>::value, DataType>::type
    FloatDivZero(DataType const f, DataType const d) {
      return IsZero<DataType>(d) ? 0.0 : f / d;
    }
  }  // namespace special

}  // namespace openKIM

#endif  // SPECIAL_HPP