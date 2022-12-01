#pragma once

// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see <https://www.gnu.org/licenses/>.

#include <Eigen/Core>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "libfly/system/adaptor.hpp"
#include "libfly/utility/core.hpp"

/**
 * \file SoA.hpp
 */

namespace fly::system {

  template <typename... Pr>
  class SoA;  // Required for detail below.

  namespace detail {

    template <typename, typename>
    struct different_SoA : std::false_type {};

    template <typename... Pr>
    struct different_SoA<SoA<Pr...>, SoA<Pr...>> : std::false_type {};

    template <typename... Pr, typename... Mx>
    struct different_SoA<SoA<Pr...>, SoA<Mx...>> : std::true_type {};

  }  // namespace detail

  /**
   * @brief A container for atoms that stores each property in its own array.
   *
   * The default container type used in the libfly; an ``SoA`` models an array of ``Atom`` types
   * but decomposes the atom type and stores each property in a separate array. This enables efficient
   * cache use. Like ``Atom``, the properties  of the "atom" are described through a series of template parameters
   * which should inherit from ``Property``. The properties  of each atom can be accessed by the index of the
   * atom or as an ``Eigen::Array`` to enable collective operations.
   *
   * SoA also supports slicing and reference properties  which transform that property into a view, this enables SoA to act as a
   * concrete type in interfaces whilst allowing implicit conversions from any compatible SoA.
   *
   * \rst
   *
   * Example:
   *
   * .. include:: ../../examples/system/SoA.cpp
   *    :code:
   *
   * \endrst
   *
   * @tparam Pr a series of empty types, derived from ``Property``, to describe each member.
   */
  template <typename... Pr>
  class SoA : public detail::Adaptor<Pr>... {
  public:
    /**
     * @brief True if this SoA is a pure view i.e. all its properties are reference types.
     */
    static constexpr bool owns_none = (std::is_reference_v<Pr> && ...);

    /**
     * @brief True if this SoA is purely owning i.e. non of its properties are reference types.
     */
    static constexpr bool owns_all = (!std::is_reference_v<Pr> && ...);

    /**
     * @brief Detect if a type is a specialization of a SoA but different from this specialization.
     */
    template <typename T>
    static constexpr bool different_SoA_v = detail::different_SoA<SoA<Pr...>, remove_cref_t<T>>::value;

    /**
     * @brief Construct a new empty SoA.
     */
    SoA() = default;

    /**
     * @brief Defaulted move constructor.
     */
    SoA(SoA&&) = default;

    /**
     * @brief Defaulted copy constructor.
     */
    SoA(SoA const&) = default;

    /**
     * @brief Construct a new SoA containing 'size' atoms.
     *
     * \rst
     * .. warning::
     *    New atoms are uninitialized, non-owning SoA's contain ``nullptr``.
     * \endrst
     */
    explicit SoA(Eigen::Index size) : detail::Adaptor<Pr>(size)..., m_size(size) {}

  private:
    /**
     * @brief Check Other is a specialization of a SoA but different from this SoA and that this SoA's Adaptors are (move)
     * constructable from Other's adaptors.
     */
    template <typename Other>
    static constexpr bool constructible = different_SoA_v<Other> && (std::is_constructible_v<detail::Adaptor<Pr>, Other> && ...);

  public:
    /**
     * @brief Implicitly construct a new SoA object from SoA 'other' with different properties .
     *
     * Only SFINE enabled if this SoA owns non of its arrays and this SoA is actually constructable from other.
     *
     */
    template <typename... T, typename = std::enable_if_t<owns_none && constructible<SoA<T...> const&>>>
    SoA(SoA<T...> const& other) : detail::Adaptor<Pr>(other)..., m_size(other.size()) {}

    /**
     * @brief Implicitly construct a new SoA object from SoA 'other' with different properties .
     *
     * Only SFINE enabled if this SoA owns non of its arrays and this SoA is actually constructable from other.
     *
     */
    template <typename... T, typename = std::enable_if_t<owns_none && constructible<SoA<T...>&>>>
    SoA(SoA<T...>& other) : detail::Adaptor<Pr>(other)..., m_size(other.size()) {}

    /**
     * @brief Implicitly construct a new SoA object from SoA 'other' with different properties .
     *
     * Only SFINE enabled if this SoA owns non of its arrays and this SoA is actually constructable from other.
     *
     */
    template <typename... T, typename = std::enable_if_t<owns_none && constructible<SoA<T...>&&>>>
    SoA(SoA<T...>&& other) : detail::Adaptor<Pr>(std::move(other))..., m_size(other.size()) {
      // OK to ``std::move`` ``other`` multiple times as the detail::Adaptor constructor will only move its
      // corresponding base slice.
    }

    /**
     * @brief Explicitly construct a new SoA object from SoA 'other' with different properties .
     *
     * SFINE enabled if this SoA owns some of its arrays and this SoA is actually constructable from other.
     *
     */
    template <typename... T, typename = std::enable_if_t<!owns_none && constructible<SoA<T...> const&>>, typename = void>
    explicit SoA(SoA<T...> const& other) : detail::Adaptor<Pr>(other)..., m_size(other.size()) {
      // The unnamed template parameter is a dummy to distinguish this from the implicit version.
    }

    /**
     * @brief Explicitly construct a new SoA object from SoA 'other' with different properties .
     *
     * SFINE enabled if this SoA owns some of its arrays and this SoA is actually constructable from other.
     *
     */
    template <typename... T, typename = std::enable_if_t<!owns_none && constructible<SoA<T...>&>>, typename = void>
    explicit SoA(SoA<T...>& other) : detail::Adaptor<Pr>(other)..., m_size(other.size()) {
      // The unnamed template parameter is a dummy to distinguish this from the implicit version.
    }

    /**
     * @brief Explicitly construct a new SoA object from SoA 'other' with different properties .
     *
     * SFINE enabled if this SoA owns some of its arrays and this SoA is actually constructable from other.
     *
     */
    template <typename... T, typename = std::enable_if_t<!owns_none && constructible<SoA<T...>&&>>, typename = void>
    explicit SoA(SoA<T...>&& other) : detail::Adaptor<Pr>(std::move(other))..., m_size(other.size()) {
      // OK to ``std::move`` ``other`` multiple times as the detail::Adaptor constructor will only move its
      // corresponding base slice. The unnamed template parameter is a dummy to distinguish this from the implicit version.
    }

    //   private:
    //     template <typename Head, typename... Tail>
    //     static constexpr Head&& first(Head&& head, Tail const&...) noexcept {
    //       return std::forward<Head>(head);
    //     }

    //     template <typename T>
    //     using add_ref = std::conditional_t<std::is_reference_v<T>, T, T const&>;

    //   public:
    //     /**
    //      * @brief Construct a SoA specifying a source for each property.
    //      *
    //      * @param args The ith source initializes the ith property.
    //      */
    //     template <bool Cond = (sizeof...(Pr) > 0), typename = std::enable_if_t<Cond>>
    //     explicit SoA(SoA<add_ref<Pr>>... args) : detail::Adaptor<Pr>(args)..., m_size(first(args...).size()) {
    //       // T -> SoA<T const&>
    //       // T& -> SoA<T&>
    //       // T const & -> SoA<T const &>
    //       ASSERT(((args.size() == m_size) && ...), "SoA constructed from different sizes", 0);
    //     }

    //     /**
    //      * @brief Construct a SoA specifying a source for each property.
    //      *
    //      * @param args The ith source initializes the ith property.
    //      */
    //     template <typename... T, typename = std::enable_if_t<(sizeof...(T) > 1 && sizeof...(T) == sizeof...(Pr))>>
    //     explicit SoA(T&&... args) : detail::Adaptor<Pr>(std::forward<T>(args))..., m_size(first(args...).size()) {
    //       // T -> SoA<T const&>
    //       // T& -> SoA<T&>
    //       // T const & -> SoA<T const &>
    //       ASSERT(((args.size() == m_size) && ...), "SoA constructed from different sizes", 0);
    //     }

    /**
     * @brief Defaulted copy assignment operator.
     *
     * \rst
     *
     * See the `note on assignment to references`_.
     *
     * \endrst
     */
    SoA& operator=(SoA const&) = default;

    /**
     * @brief Defaulted move assignment operator.
     *
     * \rst
     *
     * See the `note on assignment to references`_.
     *
     * \endrst
     */
    SoA& operator=(SoA&&) = default;

    /**
     * @brief Assign to a SoA with with different properties .
     *
     * \rst
     *
     * .. _`note on assignment to references`:
     *
     * .. note::
     *    Assignment to reference properties  follows pointer-semantics:
     *
     *    .. include:: ../../examples/system/SoA_assign.cpp
     *       :code:
     *
     * \endrst
     */
    template <typename T, typename = std::enable_if_t<different_SoA_v<T>>>
    SoA& operator=(T&& other) {
      m_size = other.m_size;
      (static_cast<void>(static_cast<detail::Adaptor<Pr>&>(*this) = std::forward<T>(other)), ...);
      return *this;
    }

    /**
     * @brief Rebind a reference property to point at ''other''.
     *
     * @param other Property tag
     * @return void
     */
    template <typename T, typename U>
    auto rebind(T, U const& other) -> std::enable_if_t<(std::is_same_v<T const&, Pr> || ...)> {
      ASSERT(size() == other.size(), "rebinding to a differently sized SoA: {} != {}", size(), other.size());
      static_cast<detail::Adaptor<T const&>&>(*this) = other;
    }

    /**
     * @brief Rebind a reference property to point at ''other''.
     *
     * @param other Property tag
     * @return void
     */
    template <typename T, typename U>
    auto rebind(T, U& other) -> std::enable_if_t<(std::is_same_v<T&, Pr> || ...) || (std::is_same_v<T const&, Pr> || ...)> {
      ASSERT(size() == other.size(), "rebinding to a differently sized SoA: {} != {}", size(), other.size());
      if constexpr ((std::is_same_v<T&, Pr> || ...)) {
        static_cast<detail::Adaptor<T&>&>(*this) = other;
      } else {
        static_cast<detail::Adaptor<T const&>&>(*this) = other;
      }
    }

    // Inherit operators

    using detail::Adaptor<Pr>::operator()...;

    using detail::Adaptor<Pr>::operator[]...;

    /**
     * @brief Get the number of atoms in the SoA.
     */
    Eigen::Index size() const noexcept { return m_size; }

    /**
     * @brief Resize the SoA, no-op if ``size() == new_size``.
     *
     * \rst
     * .. warning::
     *    Destroys all the data in the SoA, new atoms are all uninitialized.
     * \endrst
     *
     * \return ``true`` if resize occurred, ``false`` otherwise.
     */
    template <bool OwnsAll = owns_all>
    auto destructive_resize(Eigen::Index new_size) -> std::enable_if_t<OwnsAll, bool> {
      if (std::exchange(m_size, new_size) != new_size) {
        *this = SoA(new_size);
        return true;
      }
      return false;
    }

  private:
    Eigen::Index m_size = 0;

    /* clang-format off */ template <typename...>  friend class SoA; /* clang-format on */
  };

}  // namespace fly::system