#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: MPL-2.0

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Core>
#include <cstddef>

#include "libfly/utility/asserts.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file adaptor.hpp
 *
 * @brief Internal helper for SoA.
 */

namespace fly::system::detail {

  template <typename Mem> struct Adaptor {
  public:
    Adaptor() = default;

    Adaptor(Adaptor&&) noexcept = default;

    Adaptor(Adaptor const&) = default;

    // If constructing from a const view then deference their pointer. Pass views by value.
    explicit Adaptor(Adaptor<Mem const&> other) : m_data(*other.m_data_ptr) {}

    // If constructing from a view then deference their pointer. Pass views by value.
    explicit Adaptor(Adaptor<Mem&> other) : m_data(*other.m_data_ptr) {}

    // OwnsAll specific.
    explicit Adaptor(int size) : m_data(size * Mem::size()) { ASSERT(size > 0, "Invalid size"); }

    Adaptor& operator=(Adaptor const&) = default;

    Adaptor& operator=(Adaptor&&) = default;

    Adaptor& operator=(Adaptor<Mem const&> other) {
      m_data = *other.m_data_ptr;
      return *this;
    }

    Adaptor& operator=(Adaptor<Mem&> other) {
      m_data = *other.m_data_ptr;
      return *this;
    }

    /**
     * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
     *
     * This is an owning Adaptor hence, model value const-semantics.
     */
    constexpr typename Mem::matrix_ref_t operator()(Mem, int i) {
      ASSERT(i >= 0 && i < m_data.size() / Mem::size(), "Index out of bounds");

      if constexpr (Mem::is_1x1) {
        return m_data[i];
      } else {
        return typename Mem::matrix_ref_t{m_data.data() + i * Mem::size()};
      }
    }

    /**
     * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
     *
     * This is an owning Adaptor hence, model value const-semantics
     */
    constexpr typename Mem::matrix_cref_t operator()(Mem, int i) const {
      ASSERT(i >= 0 && i < m_data.size() / Mem::size(), "Index out of bounds");

      if constexpr (Mem::is_1x1) {
        return m_data[i];
      } else {
        return typename Mem::matrix_cref_t{m_data.data() + i * Mem::size()};
      }
    }

    /**
     * @brief Fetch a view of the array, tagged dispatch on Mem.
     *
     * This is an owning Adaptor hence, model value const-semantics
     */
    constexpr typename Mem::array_ref_t operator[](Mem) noexcept { return m_data; }

    /**
     * @brief Fetch a view of the array, tagged dispatch on Mem.
     *
     * This is an owning Adaptor hence, model value const-semantics
     */
    constexpr typename Mem::array_cref_t operator[](Mem) const noexcept { return m_data; }

  private:
    typename Mem::array_t m_data;  ///< Owns its own array.

    friend struct Adaptor<Mem&>;
    friend struct Adaptor<Mem const&>;

    /**
     * Possibilities
     *
     *  void foo1(Adaptor<Pos>)
     *  void foo2(Adaptor<Pos> &)
     *  void foo3(Adaptor<Pos> const &)
     *
     *  foo1(Adaptor<Pos>{})          use defaulted copy/move
     *  foo1(Adaptor<Pos &>{})        ok to do implicitly as user wants a copy
     *  foo1(Adaptor<Pos const&>{})   ok to do implicitly as user wants a copy
     *
     *  foo2(Adaptor<Pos>{})          use defaulted copy/move, compiler will not allow ref to temp
     *  foo2(Adaptor<Pos &>{})        compiler will not allow ref to temp
     *  foo2(Adaptor<Pos const&>{})   compiler will not allow ref to temp
     *
     *  foo3(Adaptor<Pos>{})          Just a const ref
     *  foo3(Adaptor<Pos &>{})        must not allow implicit construction from a view else this would create a temporary
     *  foo3(Adaptor<Pos const&>{})   same as above ^
     *
     */
  };

  template <typename Mem> struct Adaptor<Mem&> {
  public:
    Adaptor() = default;

    Adaptor(Adaptor&&) noexcept = default;

    Adaptor(Adaptor const&) = default;

    // If constructing from am owning Adaptor then take address of their m_data member, no explicit for construction of a view.
    Adaptor(Adaptor<Mem>& other) : m_data_ptr(&other.m_data) {}

    // Cannot construct from a const view.
    Adaptor(Adaptor<Mem const&> other) = delete;

    Adaptor& operator=(Adaptor const&) = default;

    Adaptor& operator=(Adaptor&&) = default;

    // Cannot assign to a const view
    Adaptor& operator=(Adaptor<Mem const&> other) = delete;

    // Can only assign to a reference to an owning view
    Adaptor& operator=(Adaptor<Mem>& other) {
      m_data_ptr = &other.m_data;
      return *this;
    }

    /**
     * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
     *
     * This is not an owning Adaptor hence, model pointer const-semantics.
     */
    constexpr typename Mem::matrix_ref_t operator()(Mem, int i) const {
      ASSERT(i >= 0 && i < m_data_ptr->size() / Mem::size(), "Index out of bounds");

      if constexpr (Mem::is_1x1) {
        return (*m_data_ptr)[i];
      } else {
        return typename Mem::matrix_ref_t{m_data_ptr->m_data() + i * Mem::size()};
      }
    }

    /**
     * @brief Fetch a view of the array, tagged dispatch on Mem.
     *
     * This is not an owning Adaptor hence, model pointer const-semantics.
     */
    constexpr typename Mem::array_ref_t operator[](Mem) const noexcept {
      ASSERT(m_data_ptr, "Deferencing an empty view adaptor.");
      return *m_data_ptr;
    }

  private:
    Eigen::ArrayBase<typename Mem::array_t>* m_data_ptr = nullptr;  ///< Pointer to the viewed array.

    friend struct Adaptor<Mem>;
    friend struct Adaptor<Mem const&>;
  };

  template <typename Mem> struct Adaptor<Mem const&> {
  public:
    Adaptor() = default;

    Adaptor(Adaptor&&) noexcept = default;

    Adaptor(Adaptor const&) = default;

    // If constructing from a non const view then just copy the pointer, no explicit for construction of a view.
    Adaptor(Adaptor<Mem&> other) : m_data_ptr(other.m_data_ptr) {}

    // If constructing from am owning Adaptor then take address of their m_data member, no explicit for construction of a view.
    Adaptor(Adaptor<Mem> const& other) : m_data_ptr(&other.m_data) {}

    Adaptor& operator=(Adaptor const&) = default;

    Adaptor& operator=(Adaptor&&) = default;

    Adaptor& operator=(Adaptor<Mem&> other) {
      m_data_ptr = other.m_data_ptr;
      return *this;
    }
    // Can only assign to a reference to an owning view
    Adaptor& operator=(Adaptor<Mem> const& other) {
      m_data_ptr = &other.m_data;
      return *this;
    }

    /**
     * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
     *
     * This is not an owning Adaptor hence, model const-pointer const-semantics.
     */
    constexpr typename Mem::matrix_cref_t operator()(Mem, int i) const {
      ASSERT(i >= 0 && i < m_data_ptr->size() / Mem::size(), "Index out of bounds");

      if constexpr (Mem::is_1x1) {
        return (*m_data_ptr)[i];
      } else {
        return typename Mem::matrix_cref_t{m_data_ptr->m_data() + i * Mem::size()};
      }
    }

    /**
     * @brief Fetch a view of the array, tagged dispatch on Mem.
     *
     * This is not an owning Adaptor hence, model const-pointer const-semantics.
     */
    constexpr typename Mem::array_cref_t operator[](Mem) const noexcept {
      ASSERT(m_data_ptr, "Deferencing an empty view adaptor.");
      return *m_data_ptr;
    }

  private:
    Eigen::ArrayBase<typename Mem::array_t> const* m_data_ptr = nullptr;  ///< Pointer to the viewed array.

    friend struct Adaptor<Mem>;
    friend struct Adaptor<Mem&>;
  };

  template <typename Mem> struct Adaptor<Mem&&> { static_assert(always_false<Mem>, "Rvalue reference to member is illegal."); };

  template <typename Mem> struct Adaptor<Mem const&&> { static_assert(always_false<Mem>, "Rvalue reference to member is illegal."); };

}  // namespace fly::system::detail