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
 * format. This file format is widely supported, see the `GSD documentation <https://gsd.readthedocs.io/>`_ .
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
    // be equal to the particles/N. This causes problems when using a TypeMap to store properties that have tags of the form
    // "log/particles/*".
    template <typename Tag>
    struct re_write_tag : Tag {
    private:
      static constexpr auto N = 12 + std::string_view{Tag::tag}.size() + 1;

      static constexpr std::array<char, N> build() {
        //
        std::array<char, N> array{'l', 'o', 'g', '/', 't', 'y', 'p', 'e', 'm', 'a', 'p', '/'};

        if (std::string_view{Tag::tag}.size() < 4) {
          throw "Tag is too short to rewrite";
        }

        char const *p = Tag::tag;

        if (p[0] == 'l' && p[1] == 'o' && p[2] == 'g' && p[3] == '/') {
          p = p + 4;
        }

        char *write = array.data() + 12;

        while (*p) {
          *write++ = *p++;
        }

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
    // Pre C++20 workaround for decltype([]{})
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
     * @brief Write the number of atoms ``num_atoms`` to the current frame.
     *
     * Despite this being seemingly redundant as every chunk in GSD stores how many records there are some visualisers (OVITO)
     * still require "particles/N" to be set.
     */
    auto write(std::uint32_t num_atoms) -> void { dump_one("particles/N", num_atoms); }

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
     * @brief Write a TypeMap to the current frame.
     */
    template <typename... U>
    auto write(system::TypeMap<U...> const &map) -> void {
      //
      dump_one("log/typemap/N", safe_cast<std::uint32_t>(map.size()));  // explicitly store for reconstruction.

      write(Type{}, static_cast<typename system::TypeMap<U...>::SOA const &>(map));

      (write(detail::re_write_tag<U>{}, static_cast<typename system::TypeMap<U...>::SOA const &>(map)), ...);
    }

    /**
     * @brief Write a tagged property of a fly::system::SoA to the current frame.
     */
    template <typename T, typename... U>
    auto write(T, system::SoA<U...> const &in) -> std::enable_if_t<std::is_arithmetic_v<typename T::scalar_t>> {
      dump_span(T::tag, T::size(),
                nonstd::span<typename T::scalar_t const>{
                    in[T{}].derived().data(),
                    safe_cast<std::size_t>(in.size() * T::size()),
                });
    }

    /**
     * @brief Read the number of atoms stored at frame ``i`` and return it.
     */
    auto read_num_atoms(std::uint64_t i) -> std::uint32_t { return load_one<std::uint32_t>(i, "particles/N"); }

    // ///////////////////////////////////////////////

    /**
     * @brief Read a Box stored at frame ``i`` and return it.
     */
    auto read_box(std::uint64_t i) const -> system::Box;

    /**
     * @brief Read a tagged Property from the ``i``th frame and write it to ``out``.
     */
    template <typename T, typename... U>
    auto read_to(std::uint64_t i, T, system::SoA<U...> &out) -> std::enable_if_t<std::is_arithmetic_v<typename T::scalar_t>> {
      load_span(i, T::tag, T::size(),
                nonstd::span<typename T::scalar_t>{
                    out[T{}].derived().data(),
                    safe_cast<std::size_t>(out.size() * T::size()),
                });
    }

    /**
     * @brief  Read a tagged Property from the ``i``th frame and write it to a TypeMap.
     */
    template <typename... U>
    auto read_map(std::uint64_t i) -> system::TypeMap<U...> {
      //
      system::TypeMap<U...> map(load_one<std::uint32_t>(i, "log/typemap/N"));

      read_to(i, Type{}, static_cast<typename system::TypeMap<U...>::SOA &>(map));

      (read_to(i, detail::re_write_tag<U>{}, static_cast<typename system::TypeMap<U...>::SOA &>(map)), ...);

      return map;
    }

    // ////////////////////////////////////////////////////////////

    ~FileGSD() noexcept;

  private:
    std::string m_fname;
    std::unique_ptr<gsd_handle> m_handle;

    // ///////////////////////////////////////////

    void end_frame();

    template <typename T>
    void dump_one(char const *name, T const &value) {
      dump_span(name, 1, nonstd::span<std::uint32_t const>{&value, 1});
    }

    template <typename T>
    T load_one(std::uint64_t i, char const *name) {
      std::array<T, 1> buf;
      load_span(i, name, 1, buf);
      return buf[0];
    }

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

    dump_load(char);  // Special handling as not a gsd-native type.

#undef dump_load
  };

}  // namespace fly::io
