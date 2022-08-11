// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <fmt/core.h>

#include <cstdint>
#include <cstring>
#include <functional>
#include <nonstd/span.hpp>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include "external/gsd.h"
#include "libfly/utility/core.hpp"

/**
 * \file
 *
 * @brief Type safe GSD abstractions.
 *
 */

namespace fly::io {

  namespace detail {

    template <typename T>
    struct get_gsd_type;

    template <>
    struct get_gsd_type<std::uint8_t> {
      static constexpr gsd_type value = GSD_TYPE_UINT8;
    };
    template <>
    struct get_gsd_type<std::uint16_t> {
      static constexpr gsd_type value = GSD_TYPE_UINT16;
    };
    template <>
    struct get_gsd_type<std::uint32_t> {
      static constexpr gsd_type value = GSD_TYPE_UINT32;
    };
    template <>
    struct get_gsd_type<std::uint64_t> {
      static constexpr gsd_type value = GSD_TYPE_UINT64;
    };
    template <>
    struct get_gsd_type<std::int8_t> {
      static constexpr gsd_type value = GSD_TYPE_INT8;
    };
    template <>
    struct get_gsd_type<std::int16_t> {
      static constexpr gsd_type value = GSD_TYPE_INT16;
    };
    template <>
    struct get_gsd_type<std::int32_t> {
      static constexpr gsd_type value = GSD_TYPE_INT32;
    };
    template <>
    struct get_gsd_type<std::int64_t> {
      static constexpr gsd_type value = GSD_TYPE_INT64;
    };
    template <>
    struct get_gsd_type<float> {
      static constexpr gsd_type value = GSD_TYPE_FLOAT;
    };
    template <>
    struct get_gsd_type<double> {
      static constexpr gsd_type value = GSD_TYPE_DOUBLE;
    };

  }  // namespace detail

  //  Convert a type to its GSD enum/tag.
  template <typename T>
  inline constexpr gsd_type get_type_id = detail::get_gsd_type<T>::value;

  //  Convert a type to its GSD enum/tag.
  template <typename T>
  inline constexpr std::string_view get_type_name = detail::get_gsd_type<T>::name;

  // Call the GSD function f, f(args...) and throw any errors.
  template <typename F, typename... Args>
  void call_gsd(std::string_view info, F&& f, Args&&... args) {
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

  // Write a chunk of data, throw any errors.
  template <typename T, gsd_type Tag = get_type_id<T>>
  void write_chunk(struct gsd_handle* handle, char const* name, uint64_t N, uint32_t M, nonstd::span<T const> data) {
    //
    XASSERT(sizeof(T) == gsd_sizeof_type(Tag), "Platform error!", 0);

    if (auto size = N * M; size != data.size()) {
      throw error("GSD: Chunk '{}', trying to write {} T's but expecting to write {}", name, data.size(), size);
    }

    call_gsd(name, gsd_write_chunk, handle, name, Tag, N, M, 0, data.data());
  }

  // Read a chunk expecting N by M into data, use -1 for unknown M or N.
  template <typename T, gsd_type Tag = get_type_id<T>>
  void read_chunk(uint64_t frame, struct gsd_handle* handle, char const* name, int N, int M, nonstd::span<T> data) {
    //
    gsd_index_entry const* chunk = gsd_find_chunk(handle, frame, name);

    if (!chunk) {
      throw error("GSD: Could not find chunk with name '{}' at frame {}", name, frame);
    }

    if (N >= 0 && safe_cast<uint64_t>(N) != chunk->N) {
      throw error("GSD: Chunk '{}', expected {} rows found {}", name, N, chunk->N);
    }

    if (M >= 0 && safe_cast<uint32_t>(M) != chunk->M) {
      throw error("GSD: Chunk '{}', expected {} columns found {}", name, M, chunk->M);
    }

    if (chunk->type != Tag) {
      throw error("GSD: Chunk '{}', expecting to read type {} but found type {}", name, Tag, chunk->type);
    }

    XASSERT(sizeof(T) == gsd_sizeof_type(Tag), "Platform error!", 0);

    if (auto size = chunk->N * chunk->M; size != data.size()) {
      throw error("Chunk {}, not enough space for {} element in span length {}", name, size, data.size());
    }

    call_gsd(name, gsd_read_chunk, handle, data.data(), chunk);
  }

}  // namespace fly::io