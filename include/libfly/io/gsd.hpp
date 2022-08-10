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
#include <memory>
#include <nonstd/span.hpp>
#include <string>
#include <string_view>
#include <utility>

#include "libfly/system/SoA.hpp"
#include "libfly/system/box.hpp"
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
    read,        ///< Open with read only permissions.
    read_write,  ///< Open with read and write permissions.
    create       ///< Open with read and write permissions, overwrite existing files.
  };

  /**
   * @brief file
   *
   */
  class FileGSD {
  public:
    /**
     * @brief Construct a new File G S D object
     *
     * @param fname
     * @param flag
     */
    explicit FileGSD(std::string_view fname, Flags flag = read);

    /**
     * @brief yes
     *
     * @tparam Args
     * @param args
     */
    template <typename... Args>
    void dump(Args const &...args) {
      (static_cast<void>(dump_impl(args)), ...);
      commit_frame();
    }

    /**
     * @brief no
     *
     * @tparam Args
     */
    template <typename... Args>
    void load(std::uint64_t i, Args &&...args) const {
      (static_cast<void>(load_impl(i, std::forward<Args>(args))), ...);
    }

    /**
     * @brief plop
     *
     * @return int
     */
    std::uint64_t n_frames() const noexcept;

    /**
     * @brief Clear the current file.
     *
     * After truncating, a file will have no frames and no data chunks. The file size will be that of a newly created gsd file. The
     * application, schema, and schema version metadata will be kept. Clear does not close and reopen the file, so it is suitable
     * for writing restart files on Lustre file.
     */
    void clear();

    ~FileGSD() noexcept;

  private:
    std::string m_fname;
    std::unique_ptr<gsd_handle> m_handle;

    void dump_impl(system::Box const &);
    void load_impl(std::uint64_t i, system::Box &) const;

    template <typename... T>
    void dump_impl(system::SoA<T...> const &soa) {
      //  ¯\_(ツ)_/¯
      (static_cast<void>(dump_span(remove_cref_t<T>::tag, remove_cref_t<T>::size(),
                                   nonstd::span<typename remove_cref_t<T>::scalar_t const>{
                                       soa[remove_cref_t<T>{}].derived().data(),
                                       safe_cast<std::size_t>(soa.size()) * remove_cref_t<T>::size(),
                                   })),
       ...);

      dump_span("particles/N", 1, std::array{safe_cast<std::uint32_t>(soa.size())});
    }

    // ///////////////////////////////////////////

    // dump a dynamic amount of the data in span in expecting rows of length M.
    void dump_span(char const *name, std::uint32_t M, nonstd::span<double const> data);
    // load a dynamic amount of data from file into span in expecting rows of length M, use -1 for dynamic.
    void load_span(std::uint64_t i, char const *name, int M, nonstd::span<double> data) const;

    void dump_span(char const *name, std::uint32_t M, nonstd::span<int const> data);
    void load_span(std::uint64_t i, char const *name, int M, nonstd::span<int> data) const;

    void dump_span(char const *name, std::uint32_t M, nonstd::span<std::uint32_t const> data);
    void load_span(std::uint64_t i, char const *name, int M, nonstd::span<std::uint32_t> data) const;

    // void load_impl(int i, system::Box &) const;

    void commit_frame();
  };

}  // namespace fly::io
