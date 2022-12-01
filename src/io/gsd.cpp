// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

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

#include "gsd/gsd/gsd.h"
#include "libfly/utility/core.hpp"

namespace fly::io {

  template <typename F, typename... Args>
  void call_gsd(std::string_view info, F &&f, Args &&...args) {
    switch (std::invoke(std::forward<F>(f), std::forward<Args>(args)...)) {
      case GSD_SUCCESS:
        return;
      case GSD_ERROR_IO:
        throw error("GSD: {} - {}", std::strerror(errno), info);
      case GSD_ERROR_INVALID_ARGUMENT:
        throw error("GSD: Invalid argument - {}", info);
      case GSD_ERROR_NOT_A_GSD_FILE:
        throw error("GSD: Not a GSD file - {}", info);
      case GSD_ERROR_INVALID_GSD_FILE_VERSION:
        throw error("GSD: Invalid GSD file version - {}", info);
      case GSD_ERROR_FILE_CORRUPT:
        throw error("GSD: File corrupt - {}", info);
      case GSD_ERROR_MEMORY_ALLOCATION_FAILED:
        throw error("GSD: Memory allocation failed - {}", info);
      case GSD_ERROR_NAMELIST_FULL:
        throw error("GSD: Name-list full - {}", info);
      case GSD_ERROR_FILE_MUST_BE_WRITABLE:
        throw error("GSD: File must be writeable - {}", info);
      case GSD_ERROR_FILE_MUST_BE_READABLE:
        throw error("GSD: File must be readable - {}", info);
      default:
        throw error("GSD: Unknown error - {}", info);
    }
  }

  BinaryFile::BinaryFile(std::string_view fname, Flags flag) : m_fname(fname), m_handle(std::make_unique<gsd_handle>()) {
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
    if constexpr (spatial_dims == 3) {
      // Special handling for 3d case for HOOMD Schema

      //          |L_x    xy L_y   xz L_z|
      // basis =  |0         L_y   yz L_z|
      //          |0         0        L_z|
      Eigen::Matrix<float, 3, 3> basis = box.basis().cast<float>();

      // L_x, L_y, L_z , xy, xz, yz
      std::array<float, 6> const hoomd_basis = {
          basis(0, 0),
          basis(1, 1),
          basis(2, 2),
          basis(0, 1),
          basis(0, 2),
          basis(1, 2),
      };

      write_chunk("configuration/box", 6, 1, hoomd_basis.data());

      std::array<std::uint8_t, 3> const periodicity = {box.periodic(0), box.periodic(1), box.periodic(2)};

      write_chunk("log/periodicity", 3, 1, periodicity.data());
    } else {
      // Fall back to non-hoomd shema
      Mat basis = box.basis();

      write_chunk("configuration/box", basis.size(), 1, basis.data());

      Arr<std::uint8_t> p;

      for (int i = 0; i < spatial_dims; i++) {
        p[i] = box.periodic(i);
      }

      write_chunk("log/periodicity", p.size(), 1, p.data());
    }
  }

  auto BinaryFile::read_box(std::uint64_t i) const -> system::Box {
    //
    if constexpr (spatial_dims == 3) {
      Eigen::Matrix<float, 3, 3> basis = Eigen::Matrix<float, 3, 3>::Zero();

      {  // Load the configuration

        std::array<float, 6> hoomd_basis;  // L_x, L_y, L_z , xy, xz, yz

        read_chunk(i, "configuration/box", 6, 1, hoomd_basis.data());

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

        read_chunk(i, "log/periodicity", 3, 1, pr.data());

        periodicity[0] = pr[0];
        periodicity[1] = pr[1];
        periodicity[2] = pr[2];
      }

      return {basis.cast<double>(), periodicity};
    } else {
      // Fall back to non-hoomd shema
      Mat basis;

      read_chunk(i, "configuration/box", basis.size(), 1, basis.data());

      Arr<std::uint8_t> p;

      read_chunk(i, "log/periodicity", p.size(), 1, p.data());

      return {basis, p.cast<bool>()};
    }
  }

  //   ////////////////////////////////////////////////////////////////////////////////////////////////

  gsd_type get_gsd_tag(bool is_floating, bool is_signed, std::size_t size) {
    if (is_floating) {
      if (size == gsd_sizeof_type(GSD_TYPE_DOUBLE)) {
        return GSD_TYPE_DOUBLE;
      }
      if (size == gsd_sizeof_type(GSD_TYPE_FLOAT)) {
        return GSD_TYPE_FLOAT;
      }
      throw error("No GSD floating point type with size {}", size);
    }
    if (is_signed) {
      if (size == gsd_sizeof_type(GSD_TYPE_INT8)) {
        return GSD_TYPE_INT8;
      }
      if (size == gsd_sizeof_type(GSD_TYPE_INT16)) {
        return GSD_TYPE_INT16;
      }
      if (size == gsd_sizeof_type(GSD_TYPE_INT32)) {
        return GSD_TYPE_INT32;
      }
      if (size == gsd_sizeof_type(GSD_TYPE_INT64)) {
        return GSD_TYPE_INT64;
      }
      throw error("No GSD signed integer type type with size {}", size);
    } else {
      if (size == gsd_sizeof_type(GSD_TYPE_UINT8)) {
        return GSD_TYPE_UINT8;
      }
      if (size == gsd_sizeof_type(GSD_TYPE_UINT16)) {
        return GSD_TYPE_UINT16;
      }
      if (size == gsd_sizeof_type(GSD_TYPE_UINT32)) {
        return GSD_TYPE_UINT32;
      }
      if (size == gsd_sizeof_type(GSD_TYPE_UINT64)) {
        return GSD_TYPE_UINT64;
      }
      throw error("No GSD unsigned integer type type with size {}", size);
    }
  }

  void BinaryFile::write_chunk(char const *name, std::uint64_t N, std::uint32_t M, void const *data, type_info x) {
    call_gsd(name, gsd_write_chunk, m_handle.get(), name, get_gsd_tag(x.is_floating, x.is_signed, x.size), N, M, 0, data);
  }

  void BinaryFile::read_chunk(std::uint64_t frame, char const *name, std::uint64_t N, std::uint32_t M, void *data, type_info x) const {
    //
    gsd_index_entry const *chunk = gsd_find_chunk(m_handle.get(), frame, name);

    // Fallback to initial frame
    if (!chunk) {
      if (frame == 0) {
        throw error("GSD: Could not find chunk with name '{}' at frame {}", name, frame);
      } else {
        read_chunk(0, name, N, M, data, x);
      }
    } else {
      if (N != chunk->N) {
        throw error("GSD: Chunk '{}' frame {}, expected {} rows found {}", name, frame, N, chunk->N);
      }

      if (M != chunk->M) {
        throw error("GSD: Chunk '{}' frame {}, expected {} columns found {}", name, frame, M, chunk->M);
      }

      if (gsd_type tag = get_gsd_tag(x.is_floating, x.is_signed, x.size); chunk->type != tag) {
        throw error("GSD: Chunk '{}' frame {}, expecting to read type {} but found type {}", name, frame, tag, chunk->type);
      }

      call_gsd(name, gsd_read_chunk, m_handle.get(), data, chunk);
    }
  }

}  // namespace fly::io
