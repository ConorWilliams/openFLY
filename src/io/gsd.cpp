// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include "libfly/io/gsd.hpp"

#include <fmt/core.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <nonstd/span.hpp>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include "external/gsd.h"
#include "io/gsd_helpers.hpp"
#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

namespace fly::io {

  FileGSD::FileGSD(std::string_view fname, Flags flag) : m_fname(fname), m_handle(std::make_unique<gsd_handle>()) {
    //
    if (spatial_dims != 3) {
      throw std::runtime_error("GSD files are currently only supported in 3D");
    }
    //
    auto version = gsd_make_version(1, 4);

    switch (flag) {
      case read:
        call_gsd(m_fname, gsd_open, m_handle.get(), m_fname.c_str(), GSD_OPEN_READONLY);
        break;
      case read_write:
        call_gsd(m_fname, gsd_open, m_handle.get(), m_fname.c_str(), GSD_OPEN_READWRITE);
        break;
      case create:
        call_gsd(m_fname, gsd_create_and_open, m_handle.get(), m_fname.c_str(), "openFLY", "hoomd", version, GSD_OPEN_READWRITE, 0);
        break;
    }
  }

  void FileGSD::clear() { call_gsd(m_fname, gsd_truncate, m_handle.get()); }

  int FileGSD::n_frames() const noexcept { return safe_cast<int>(gsd_get_nframes(m_handle.get())); }

  void FileGSD::commit_frame() { call_gsd(m_fname, gsd_end_frame, m_handle.get()); }

  FileGSD::~FileGSD() noexcept { call_gsd(m_fname, gsd_close, m_handle.get()); }

  //   ////////////////////////////////////////////////////////////////////////////////////////////////

  void FileGSD::dump_span(char const *name, std::uint32_t M, nonstd::span<double const> data) {
    //  This is not a public facing function so we do not throw an error.
    XASSERT(data.size() % M == 0, "{} not divisible by {}", data.size(), M);

    write_chunk(m_handle.get(), name, data.size() / M, M, data);
  }

  void FileGSD::load_span(int i, char const *name, std::uint32_t M, nonstd::span<double> data) const {
    read_chunk(i, m_handle.get(), name, -1, safe_cast<int>(M), data);
  }

  void FileGSD::dump_span(char const *name, std::uint32_t M, nonstd::span<int const> data) {
    //  This is not a public facing function so we do not throw an error.
    XASSERT(data.size() % M == 0, "{} not divisible by {}", data.size(), M);

    write_chunk(m_handle.get(), name, data.size() / M, M, data);
  }

  void FileGSD::load_span(int i, char const *name, std::uint32_t M, nonstd::span<int> data) const {
    read_chunk(i, m_handle.get(), name, -1, safe_cast<int>(M), data);
  }

  void FileGSD::dump_span(char const *name, std::uint32_t M, nonstd::span<std::uint32_t const> data) {
    //  This is not a public facing function so we do not throw an error.
    XASSERT(data.size() % M == 0, "{} not divisible by {}", data.size(), M);

    write_chunk(m_handle.get(), name, data.size() / M, M, data);
  }

  void FileGSD::load_span(int i, char const *name, std::uint32_t M, nonstd::span<std::uint32_t> data) const {
    read_chunk(i, m_handle.get(), name, -1, safe_cast<int>(M), data);
  }

  //   ////////////////////////////////////////////////////////////////////////////////////////////////

  void FileGSD::dump_impl(system::Box const &box) {
    //          |L_x    xy L_y   xz L_z|
    // basis =  |0         L_y   yz L_z|
    //          |0         0        L_z|
    Mat<double> basis = box.basis();

    // L_x, L_y, L_z , xy, xz, yz
    std::array<float, 6> const hoomd_basis = {
        (float)basis(0, 0), (float)basis(1, 1), (float)basis(2, 2), (float)basis(0, 1), (float)basis(0, 2), (float)basis(1, 2),
    };

    write_chunk<float>(m_handle.get(), "configuration/box", 6, 1, hoomd_basis);

    std::array<std::uint8_t, 3> const periodicity = {
        box.periodic(0),
        box.periodic(1),
        box.periodic(2),
    };

    write_chunk<std::uint8_t>(m_handle.get(), "log/periodicity", 3, 1, periodicity);
  }

  void FileGSD::load_impl(int i, system::Box &box) const {
    //

    Mat<double> basis = Mat<double>::Zero();

    {  // Load the configuration

      std::array<double, 6> hoomd_basis;  // L_x, L_y, L_z , xy, xz, yz

      read_chunk<double>(i, m_handle.get(), "configuration/box", 6, 1, hoomd_basis);

      basis(0, 0) = hoomd_basis[0];
      basis(1, 1) = hoomd_basis[1];
      basis(2, 2) = hoomd_basis[2];

      basis(0, 1) = hoomd_basis[3];
      basis(0, 2) = hoomd_basis[4];
      basis(1, 2) = hoomd_basis[5];
    }
    // Load the periodicity

    Arr<bool> periodicity = Arr<bool>::Zero();

    {
      std::array<std::uint8_t, 3> pr;

      read_chunk<std::uint8_t>(i, m_handle.get(), "log/periodicity", 3, 1, pr);

      periodicity[0] = pr[0];
      periodicity[1] = pr[1];
      periodicity[2] = pr[2];
    }

    box = system::Box(basis, periodicity);
  }

}  // namespace fly::io