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

#include "io/gsd_helpers.hpp"
#include "libfly/utility/core.hpp"

namespace fly::io {

  BinaryFile::BinaryFile(std::string_view fname, Flags flag) : m_fname(fname), m_handle(std::make_unique<gsd_handle>()) {
    //
    if (spatial_dims != 3) {
      throw std::runtime_error("GSD files are currently only supported in 3D");
    }
    //
    auto version = gsd_make_version(1, 4);

    switch (flag) {
      case read_only:
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

  void BinaryFile::clear() { call_gsd(m_fname, gsd_truncate, m_handle.get()); }

  std::uint64_t BinaryFile::n_frames() const noexcept { return gsd_get_nframes(m_handle.get()); }

  void BinaryFile::end_frame() { call_gsd(m_fname, gsd_end_frame, m_handle.get()); }

  BinaryFile::~BinaryFile() noexcept { call_gsd(m_fname, gsd_close, m_handle.get()); }

  //   ////////////////////////////////////////////////////////////////////////////////////////////////

  void BinaryFile::write(system::Box const &box) {
    //          |L_x    xy L_y   xz L_z|
    // basis =  |0         L_y   yz L_z|
    //          |0         0        L_z|
    Eigen::Matrix<float, 3, 3> basis = box.basis().cast<float>();

    // L_x, L_y, L_z , xy, xz, yz
    std::array<float, 6> const hoomd_basis = {
        basis(0, 0), basis(1, 1), basis(2, 2), basis(0, 1), basis(0, 2), basis(1, 2),
    };

    write_chunk<float>(m_handle.get(), "configuration/box", 6, 1, hoomd_basis);

    std::array<std::uint8_t, 3> const periodicity = {box.periodic(0), box.periodic(1), box.periodic(2)};

    write_chunk<std::uint8_t>(m_handle.get(), "log/periodicity", 3, 1, periodicity);
  }

  auto BinaryFile::read_box(std::uint64_t i) const -> system::Box {
    //

    Eigen::Matrix<float, 3, 3> basis = Eigen::Matrix<float, 3, 3>::Zero();

    {  // Load the configuration

      std::array<float, 6> hoomd_basis;  // L_x, L_y, L_z , xy, xz, yz

      read_chunk<float>(i, m_handle.get(), "configuration/box", 6, 1, hoomd_basis);

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

    return {basis.cast<double>(), periodicity};
  }

  //   ////////////////////////////////////////////////////////////////////////////////////////////////

#define dump_load(TYPE_NAME)                                                                                 \
  void BinaryFile::dump_span(char const *name, std::uint32_t M, nonstd::span<TYPE_NAME const> data) {        \
    ASSERT(data.size() % M == 0, "Chunk `{}`, {} not divisible by {}", name, data.size(), M);                \
    write_chunk(m_handle.get(), name, data.size() / M, M, data);                                             \
  }                                                                                                          \
                                                                                                             \
  void BinaryFile::load_span(std::uint64_t i, char const *name, int M, nonstd::span<TYPE_NAME> data) const { \
    read_chunk(i, m_handle.get(), name, -1, M, data);                                                        \
  }

  dump_load(std::uint8_t);
  dump_load(std::uint16_t);
  dump_load(std::uint32_t);
  dump_load(std::uint64_t);

  dump_load(std::int8_t);
  dump_load(std::int16_t);
  dump_load(std::int32_t);
  dump_load(std::int64_t);

  dump_load(float);
  dump_load(double);

  dump_load(char);

#undef dump_load

}  // namespace fly::io
