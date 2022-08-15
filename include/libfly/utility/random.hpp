#pragma once

// The code in this file is adapted from the original implementation: http://prng.di.unimi.it/xoshiro256starstar.c

// Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

// To the extent possible under law, the author has dedicated all copyright
// and related and neighbouring rights to this software to the public domain
// worldwide. This software is distributed without any warranty.

// See <http://creativecommons.org/publicdomain/zero/1.0/>.

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "libfly/utility/core.hpp"

/**
 * \file random.hpp
 *
 * @brief Pseudo random number generators (PRNG).
 */

namespace fly {

  /**
   * @brief A \<random\> compatible implementation of the xoshiro256** 1.0 PRNG
   *
   * \rst
   *
   * From `the original <https://prng.di.unimi.it/>`_:
   *
   *    This is xoshiro256** 1.0, one of our all-purpose, rock-solid generators. It has excellent (sub-ns) speed, a state (256 bits)
   *    that is large enough for any parallel application, and it passes all tests we are aware of.
   *
   * \endrst
   */
  class Xoshiro {
  private:
    /**
     * @brief Utility function.
     */
    [[nodiscard]] static constexpr std::uint64_t rotl(std::uint64_t const x, int const k) noexcept {
      return (x << k) | (x >> (64 - k));
    }

  public:
    /**
     * @brief Construct and seed the PRNG.
     *
     * The state must be seeded so that it is not everywhere zero.
     */
    explicit constexpr Xoshiro(std::array<std::uint64_t, 4> const& seed) : m_state{seed} {
      if (seed == std::array<std::uint64_t, 4>{0, 0, 0, 0}) {
        throw error("{} is a bad seed!", seed);
      }
    }

    /**
     * @brief Get the minimum value of the generator.
     */
    static constexpr auto min() noexcept -> std::uint64_t { return std::numeric_limits<std::uint64_t>::lowest(); }

    /**
     * @brief Get the maximum value of the generator.
     */
    static constexpr auto max() noexcept -> std::uint64_t { return std::numeric_limits<std::uint64_t>::max(); }

    /**
     * @brief Generate a random bit sequence and advance the state of the generator.
     */
    constexpr auto operator()() noexcept -> std::uint64_t {
      std::uint64_t const result = rotl(m_state[1] * 5, 7) * 9;

      std::uint64_t const t = m_state[1] << 17;

      m_state[2] ^= m_state[0];
      m_state[3] ^= m_state[1];
      m_state[1] ^= m_state[2];
      m_state[0] ^= m_state[3];

      m_state[2] ^= t;

      m_state[3] = rotl(m_state[3], 45);

      return result;
    }

    /**
     * @brief This is the jump function for the generator.
     *
     * It is equivalent to 2^128 calls to operator(); it can be used to generate 2^128 non-overlapping sub-sequences for parallel
     * computations.
     */
    constexpr auto jump() noexcept -> void {
      jump_impl({0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c});
    }

    /**
     * @brief This is the long-jump function for the generator.
     *
     * It is equivalent to 2^192 calls to operator(); it can be used to generate 2^64 starting points, from each of which jump() will
     * generate 2^64 non-overlapping sub-sequences for parallel distributed computations.
     */
    constexpr auto long_jump() noexcept -> void {
      jump_impl({0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635});
    }

  private:
    std::array<std::uint64_t, 4> m_state;

    void jump_impl(std::array<std::uint64_t, 4> const& JUMP) noexcept {
      //
      std::array<std::uint64_t, 4> s = {0, 0, 0, 0};

      for (std::uint64_t jump : JUMP) {
        for (int b = 0; b < 64; ++b) {
          if (jump & std::uint64_t{1} << b) {
            s[0] ^= m_state[0];
            s[1] ^= m_state[1];
            s[2] ^= m_state[2];
            s[3] ^= m_state[3];
          }
          operator()();
        }
      }
      m_state = s;
    }
  };

}  // namespace fly
