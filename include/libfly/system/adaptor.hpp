#pragma once

// Copyright © 2020 Conor Williams <conorwilliams@outlooK.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <Eigen/Core>
#include <cstddef>

#include "libfly/utility/core.hpp"

/**
 * \file adaptor.hpp
 *
 * @brief Internal helper for SoA.
 */

namespace fly::system::detail {

  template <typename T>
  struct Adaptor {
  public:
    Adaptor() = default;

    Adaptor(Adaptor&&) noexcept = default;

    Adaptor(Adaptor const&) = default;

    // If constructing from a const view then deference their pointer. Pass views by value.
    explicit Adaptor(Adaptor<T const&> other) : m_data(*other.m_data_ptr) {}

    // If constructing from a view then deference their pointer. Pass views by value.
    explicit Adaptor(Adaptor<T&> other) : m_data(*other.m_data_ptr) {}

    // OwnsAll specific.
    explicit Adaptor(Eigen::Index size) : m_data(size * T::size()) {
      //
      verify(size > 0, "Invalid size {}", size);
    }

    Adaptor& operator=(Adaptor const&) = default;

    Adaptor& operator=(Adaptor&&) = default;

    Adaptor& operator=(Adaptor<T const&> other) {
      m_data = *other.m_data_ptr;
      return *this;
    }

    Adaptor& operator=(Adaptor<T&> other) {
      m_data = *other.m_data_ptr;
      return *this;
    }

    /**
     * @brief Fetch a view of the ``i``th property stored in the array, tagged dispatch on ``T``.
     *
     * This is an owning Adaptor hence, model value const-semantics.
     */
    constexpr typename T::matrix_ref_t operator()(T, Eigen::Index i) {
      //
      ASSERT(i >= 0 && i < m_data.size() / T::size(), "Index {} is !< {}", i, m_data.size() / T::size());

      if constexpr (T::is_1x1) {
        return m_data[i];
      } else {
        return typename T::matrix_ref_t{m_data.data() + i * T::size()};
      }
    }

    /**
     * @brief Fetch a view of the ``i``th property stored in the array, tagged dispatch on T.
     *
     * This is an owning Adaptor hence, model value const-semantics
     */
    constexpr typename T::matrix_cref_t operator()(T, Eigen::Index i) const {
      //
      ASSERT(i >= 0 && i < m_data.size() / T::size(), "Index {} is !< {}", i, m_data.size() / T::size());

      if constexpr (T::is_1x1) {
        return m_data[i];
      } else {
        return typename T::matrix_cref_t{m_data.data() + i * T::size()};
      }
    }

    /**
     * @brief Fetch a view of the array, tagged dispatch on T.
     *
     * This is an owning Adaptor hence, model value const-semantics
     */
    constexpr typename T::array_ref_t operator[](T) noexcept { return {m_data.data(), m_data.rows()}; }

    /**
     * @brief Fetch a view of the array, tagged dispatch on T.
     *
     * This is an owning Adaptor hence, model value const-semantics
     */
    constexpr typename T::array_cref_t operator[](T) const noexcept { return m_data; }

  private:
    typename T::array_t m_data;  ///< Owns its own array.

    friend struct Adaptor<T&>;
    friend struct Adaptor<T const&>;

    /**
     * Possibilities
     *
     *  void foo1(Adaptor<Position>)
     *  void foo2(Adaptor<Position> &)
     *  void foo3(Adaptor<Position> const &)
     *
     *  foo1(Adaptor<Position>{})          use defaulted copy/move
     *  foo1(Adaptor<Position &>{})        OK to do implicitly as user wants a copy
     *  foo1(Adaptor<Position const&>{})   OK to do implicitly as user wants a copy
     *
     *  foo2(Adaptor<Position>{})          use defaulted copy/move, compiler will not allow ref to temp
     *  foo2(Adaptor<Position &>{})        compiler will not allow ref to temp
     *  foo2(Adaptor<Position const&>{})   compiler will not allow ref to temp
     *
     *  foo3(Adaptor<Position>{})          Just a const ref
     *  foo3(Adaptor<Position &>{})        must not allow implicit construction from a view else this would create a temporary
     *  foo3(Adaptor<Position const&>{})   same as above ^
     *
     */
  };

  template <typename T>
  struct Adaptor<T&> {
  public:
    Adaptor() = default;

    Adaptor(Adaptor&&) noexcept = default;

    Adaptor(Adaptor const&) = default;

    // If constructing from am owning Adaptor then take address of their m_data property, no explicit for construction of a view.
    Adaptor(Adaptor<T>& other) : m_data_ptr(&other.m_data) {}

    // Cannot construct from a const view.
    Adaptor(Adaptor<T const&> other) = delete;

    Adaptor& operator=(Adaptor const&) = default;

    Adaptor& operator=(Adaptor&&) = default;

    // Cannot assign to a const view
    Adaptor& operator=(Adaptor<T const&> other) = delete;

    // Can only assign to a reference to an owning view
    Adaptor& operator=(Adaptor<T>& other) {
      m_data_ptr = &other.m_data;
      return *this;
    }

    /**
     * @brief Fetch a view of the ``i``th property stored in the array, tagged dispatch on T.
     *
     * This is not an owning Adaptor hence, model pointer const-semantics.
     */
    constexpr typename T::matrix_ref_t operator()(T, Eigen::Index i) const {
      //
      ASSERT(i >= 0 && i < m_data_ptr->size() / T::size(), "Index {} is !< {}", i, m_data_ptr->size() / T::size());

      if constexpr (T::is_1x1) {
        return (*m_data_ptr)[i];
      } else {
        return typename T::matrix_ref_t{m_data_ptr->data() + i * T::size()};
      }
    }

    /**
     * @brief Fetch a view of the array, tagged dispatch on T.
     *
     * This is not an owning Adaptor hence, model pointer const-semantics.
     */
    constexpr typename T::array_ref_t operator[](T) const {
      //
      ASSERT(m_data_ptr, "Dereferencing an empty view adaptor.", 0);

      return typename T::array_ref_t{m_data_ptr->data(), m_data_ptr->rows()};
    }

  private:
    typename T::array_t* m_data_ptr = nullptr;  ///< Pointer to the viewed array.

    friend struct Adaptor<T>;
    friend struct Adaptor<T const&>;
  };

  template <typename T>
  struct Adaptor<T const&> {
  public:
    Adaptor() = default;

    Adaptor(Adaptor&&) noexcept = default;

    Adaptor(Adaptor const&) = default;

    // If constructing from a non const view then just copy the pointer, no explicit for construction of a view.
    Adaptor(Adaptor<T&> other) : m_data_ptr(other.m_data_ptr) {}

    // If constructing from am owning Adaptor then take address of their m_data property, no explicit for construction of a view.
    Adaptor(Adaptor<T> const& other) : m_data_ptr(&other.m_data) {}

    Adaptor& operator=(Adaptor const&) = default;

    Adaptor& operator=(Adaptor&&) = default;

    Adaptor& operator=(Adaptor<T&> other) {
      m_data_ptr = other.m_data_ptr;
      return *this;
    }
    // Can only assign to a reference to an owning view
    Adaptor& operator=(Adaptor<T> const& other) {
      m_data_ptr = &other.m_data;
      return *this;
    }

    /**
     * @brief Fetch a view of the ``i``th property stored in the array, tagged dispatch on T.
     *
     * This is not an owning Adaptor hence, model const-pointer const-semantics.
     */
    constexpr typename T::matrix_cref_t operator()(T, Eigen::Index i) const {
      //
      ASSERT(i >= 0 && i < m_data_ptr->size() / T::size(), "Index {} is !< {}", i, m_data_ptr->size() / T::size());

      if constexpr (T::is_1x1) {
        return (*m_data_ptr)[i];
      } else {
        return typename T::matrix_cref_t{m_data_ptr->data() + i * T::size()};
      }
    }

    /**
     * @brief Fetch a view of the array, tagged dispatch on T.
     *
     * This is not an owning Adaptor hence, model const-pointer const-semantics.
     */
    constexpr typename T::array_cref_t operator[](T) const {
      ASSERT(m_data_ptr, "Dereferencing an empty view adaptor.", 0);
      return *m_data_ptr;
    }

  private:
    typename T::array_t const* m_data_ptr = nullptr;  ///< Pointer to the viewed array.

    friend struct Adaptor<T>;
    friend struct Adaptor<T&>;
  };

  template <typename T>
  struct Adaptor<T&&> {
    static_assert(always_false<T>, "Rvalue reference to property is illegal.");
  };

  template <typename T>
  struct Adaptor<T const&&> {
    static_assert(always_false<T>, "Rvalue reference to property is illegal.");
  };

}  // namespace fly::system::detail