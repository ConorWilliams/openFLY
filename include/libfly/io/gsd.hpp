#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <nonstd/span.hpp>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file gsd.hpp
 *
 * @brief Binary IO support in the GSD format.
 *
 * \rst
 *
 * LibFLY supports writing `HOOMD-blue's <https://github.com/glotzerlab/hoomd-blue>`_ `GSD <https://github.com/glotzerlab/gsd>`_ file
 * format. This file format is widely supported, see the `GSD documentation <https://gsd.readthedocs.io/>`_ for details. In summary GSD
 * is a binary file with efficient random access to frames. GSD allows all particle properties to vary from one frame to the next.
 *
 * The GSD format models a hierarchical data structure, the outermost hierarchy is the **frame** which are index numerically. Each
 * frame may contain one or more data **chunks**, chunks are indexed by strings and map to **values** - ``N`` by ``M`` arrays of data.
 * Values are well defined for all fields at all frames. When a data chunk is present in frame ``i``, it defines the values for the
 * frame. When it is not present, the data chunk of the same name at frame ``0`` defines the values for frame ``i``. When ``N`` is
 * zero, an index entry may be written for a data chunk with no actual data written to the file for that chunk.
 *
 * \endrst
 */

/**
 * @brief External type
 */
struct gsd_handle;

namespace fly::io {

  namespace detail {

    // This type is a work around for OVITO, any chunk named "log/particles/*" is assumed to be a per-particle quantity and thus N must
    // be equal to particles/N. This causes problems when using a TypeMap to store properties that have tags of the form
    // "log/particles/*".
    template <typename Tag>
    struct re_write_tag : Tag {
    private:
      static constexpr auto N = 12 + std::string_view{Tag::tag}.size() + 1;

      static constexpr std::array<char, N> build() {
        //
        std::array<char, N> array{'l', 'o', 'g', '/', 't', 'y', 'p', 'e', 'm', 'a', 'p', '/'};

        char const *p = Tag::tag;

        if (std::string_view{Tag::tag}.size() >= 4) {
          if (p[0] == 'l' && p[1] == 'o' && p[2] == 'g' && p[3] == '/') {
            p = p + 4;
          }
        }

        char *write = array.data() + 12;

        while (*p) {
          *write++ = *p++;
        }

        *write = '\0';

        return array;
      }

      static constexpr auto arr = build();

    public:
      static constexpr char const *tag = arr.data();
    };

    static_assert(re_write_tag<Axis>::tag == std::string_view{"log/typemap/particles/axis"});
    static_assert(re_write_tag<Velocity>::tag == std::string_view{"log/typemap/particles/velocity"});

  }  // namespace detail

  /**
   * @brief Flags for file permissions.
   */
  enum Flags {
    read_only,   ///< Open with read only permissions.
    read_write,  ///< Open with read and write permissions.
    create       ///< Open with read and write permissions, overwrite existing files.
  };

  /**
   * @brief A handle to a GSD formatted binary file.
   */
  class BinaryFile {
  private:
    // Pre C++20 workaround for decltype([]{})
    struct noop {
      constexpr void operator()() const noexcept {}
    };

  public:
    /**
     * @brief Open a GSD file.
     *
     * @param fname The name of the file.
     * @param flag Read/Write/Create permissions.
     */
    explicit BinaryFile(std::string_view fname, Flags flag = read_only);

    /**
     * @brief Get the number of frames in the GSD file.
     */
    auto n_frames() const noexcept -> std::uint64_t;

    /**
     * @brief Clear the current file.
     *
     * After truncating, a file will have no frames and no data chunks. The file size will be that of a newly created gsd file. The
     * application, schema, and schema version metadata will be kept. Clear does not close and reopen the file, so it is suitable
     * for writing restart files on Lustre file.
     */
    auto clear() -> void;

    /**
     * @brief Commit the current frame to disk.
     *
     * \rst
     *
     * Example:
     *
     * .. include:: ../../examples/io/io.cpp
     *    :code:
     *
     * \endrst
     *
     * @param f An optional nullery-invokable to call before committing the writes.
     */
    template <typename F = noop>
    auto commit(F &&f = {}) -> void {
      std::invoke(std::forward<F>(f));
      end_frame();
    }

    /**
     * @brief Write a chunk to the current frame.
     *
     * @param name The chunk label (null terminated).
     * @param value Value to write to the chunk.
     */
    template <typename T>
    auto write(char const *name, T const &value) -> std::enable_if_t<std::is_arithmetic_v<T>> {
      write_chunk(name, 1, 1, &value);
    }

    /**
     * @brief Write a fly::system::Box to the current frame.
     *
     * \rst
     * .. warning::
     *     This is a lossy operation as GSD schema requires the basis vector must be stored as single precision floats.
     * \endrst
     */
    auto write(system::Box const &box) -> void;

    /**
     * @brief Write a fly::system::TypeMap to the current frame.
     */
    template <typename... U>
    auto write(system::TypeMap<U...> const &map) -> void {
      //
      write("log/typemap/N", map.num_types());  // explicitly store for reconstruction.

      write(Type{}, static_cast<typename system::TypeMap<U...>::SOA const &>(map));

      (write(detail::re_write_tag<U>{}, static_cast<typename system::TypeMap<U...>::SOA const &>(map)), ...);
    }

    /**
     * @brief Write a tagged property of a fly::system::SoA to the current frame.
     *
     * @param tag Property of ``in`` you would like to write.
     * @param in A fly::system::SoA containing the data to write to the file.
     */
    template <typename T, typename... U>
    auto write(T tag, system::SoA<U...> const &in) -> std::enable_if_t<std::is_arithmetic_v<typename T::scalar_t>> {
      write_chunk(T::tag, in.size(), T::size(), in[tag].derived().data());
    }

    /**
     * @brief Read and return a value stored in the file.
     *
     * @param i Index of frame to read from.
     * @param name Name of the chunk to read from (null terminated).
     */
    template <typename T>
    auto read(std::uint64_t i, char const *name) -> std::enable_if_t<std::is_arithmetic_v<T>, T> {
      T tmp;
      read_chunk(i, name, 1, 1, &tmp);
      return tmp;
    }

    // ///////////////////////////////////////////////

    /**
     * @brief Read a fly::system::Box stored at frame ``i`` and return it.
     *
     * @param i Index of frame to read from.
     */
    auto read_box(std::uint64_t i) const -> system::Box;

    /**
     * @brief Read a property to a fly::system::SoA.
     *
     * @param i Index of frame to read from.
     * @param tag Property of you would like to read.
     * @param out A fly::system::SoA to write the read data to.
     */
    template <typename T, typename... U>
    auto read_to(std::uint64_t i, T tag, system::SoA<U...> &out) -> std::enable_if_t<std::is_arithmetic_v<typename T::scalar_t>> {
      read_chunk(i, T::tag, out.size(), T::size(), out[tag].data());
    }

    /**
     * @brief Read a Box stored at frame ``i`` and return it.
     *
     * @param i Index of frame to read from.
     */
    template <typename... U>
    auto read_map(std::uint64_t i) -> system::TypeMap<U...> {
      //
      system::TypeMap<U...> map(read<Eigen::Index>(i, "log/typemap/N"));

      read_to(i, Type{}, static_cast<typename system::TypeMap<U...>::SOA &>(map));

      (read_to(i, detail::re_write_tag<U>{}, static_cast<typename system::TypeMap<U...>::SOA &>(map)), ...);

      return map;
    }

    // ////////////////////////////////////////////////////////////

    ~BinaryFile() noexcept;

  private:
    std::string m_fname;
    std::unique_ptr<gsd_handle> m_handle;

    // ///////////////////////////////////////////

    void end_frame();

    struct type_info {
      bool is_floating;
      bool is_signed;
      std::size_t size;
    };

    template <typename T>
    static type_info get_info() noexcept {
      return {std::is_floating_point_v<T>, std::is_signed_v<T>, sizeof(T)};
    }

    template <typename T, typename U, typename V>
    auto write_chunk(char const *name, U N, V M, T const *data) -> void {
      write_chunk(name, safe_cast<std::uint64_t>(N), safe_cast<std::uint32_t>(M), data, get_info<T>());
    }

    template <typename T, typename U, typename V>
    auto read_chunk(std::uint64_t frame, char const *name, U N, V M, T *data) const -> void {
      read_chunk(frame, name, safe_cast<std::uint64_t>(N), safe_cast<std::uint32_t>(M), data, get_info<T>());
    }

    // //////////////////////////////////////////

    auto write_chunk(char const *name, std::uint64_t N, std::uint32_t M, void const *data, type_info info) -> void;

    auto read_chunk(std::uint64_t frame, char const *name, std::uint64_t N, std::uint32_t M, void *data, type_info info) const -> void;
  };

}  // namespace fly::io
