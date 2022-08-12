#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <array>
#include <cstdint>
#include <functional>
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
 * format. This file format is widely supported, see the `GSD <https://gsd.readthedocs.io/>`_ documentation.
 *
 * \endrst
 */

/**
 * @brief External type
 */
struct gsd_handle;

namespace fly::io {

  /**
   * @brief Flags for file permissions
   */
  enum Flags {
    read_only,   ///< Open with read only permissions.
    read_write,  ///< Open with read and write permissions.
    create       ///< Open with read and write permissions, overwrite existing files.
  };

  /**
   * @brief file
   *
   */
  class FileGSD {
  private:
    struct noop {
      constexpr void operator()() const noexcept {}
    };

  public:
    /**
     * @brief Construct a new File G S D object
     *
     * @param fname
     * @param flag
     */
    explicit FileGSD(std::string_view fname, Flags flag = read_only);

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
     * @brief Commit the current writes to disk
     *
     * @tparam F An optional nullery-invokable to call before committing the writes.
     */
    template <typename F = noop>
    auto commit(F &&f = {}) -> void {
      std::invoke(std::forward<F>(f));
      end_frame();
    }

    /**
     * @brief Write a Box to the current frame.
     *
     * \rst
     * .. warning::
     *     Due to the GSD schema specification the basis vector must be stored as single precision floats so writing is a lossy
     *     operation.
     * \endrst
     *
     */
    auto write(system::Box const &box) -> void;

    /**
     * @brief Read a Box stored at frame ``i`` and write it to ``out``.
     */
    auto read(std::uint64_t i, system::Box &out) const -> void;

    /**
     * @brief Write a tagged property to the current frame.
     */
    template <typename T, typename... U>
    auto write(T, system::SoA<U...> const &in) -> std::enable_if_t<std::is_arithmetic_v<typename T::scalar_t>> {
      dump_span(T::tag, T::size(),
                nonstd::span<typename T::scalar_t const>{
                    in[T{}].derived().data(),
                    safe_cast<std::size_t>(in.size()) * T::size(),
                });
    }

    /**
     * @brief Read a tagged Property from the ``i``th frame and write it to ``out``.
     */
    template <typename T, typename... U>
    auto read(std::uint64_t i, T, system::SoA<U...> &out) -> std::enable_if_t<std::is_arithmetic_v<typename T::scalar_t>> {
      load_span(i, T::tag, T::size(),
                nonstd::span<typename T::scalar_t>{
                    out[T{}].derived().data(),
                    safe_cast<std::size_t>(out.size()) * T::size(),
                });
    }

    /**
     * @brief Write a TypeMap to the current frame.
     */
    template <typename T, typename... U>
    auto write(system::TypeMap<U...> const &map) -> void {
      (write(U{}, static_cast<typename system::TypeMap<U...>::SOA const &>(map)), ...);
    }

    /**
     * @brief  Read a tagged Property from the ``i``th frame and write it to a TypeMap.
     */
    template <typename T, typename... U>
    auto read(std::uint64_t i, system::TypeMap<U...> &map) -> void {
      (read(i, U{}, static_cast<typename system::TypeMap<U...>::SOA &>(map)), ...);
    }

    // ////////////////////////////////////////////////////////////

    ~FileGSD() noexcept;

  private:
    std::string m_fname;
    std::unique_ptr<gsd_handle> m_handle;

    // ///////////////////////////////////////////

    void end_frame();

/**
 * @brief Define a [dump/load]_span for type ``type``.
 */
#define dump_load(TYPE_NAME)                                                             \
  void dump_span(char const *name, std::uint32_t M, nonstd::span<TYPE_NAME const> data); \
  void load_span(std::uint64_t i, char const *name, int M, nonstd::span<TYPE_NAME> data) const

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

#undef dump_load
  };

}  // namespace fly::io
