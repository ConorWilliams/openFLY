// #pragma once

// #include <cstddef>
// #include <exception>
// #include <memory>

// // #include "libatom/neighbour/list.hpp"
// // #include "libatom/sim_cell.hpp"

// namespace otf::potentials {

//   /**
//    * @brief thown from parent classes if a pure virtual function is not supported.
//    */
//   struct unsupported : std::exception {};

//   /**
//    * @brief Specifies the virtual-interface for potentials in OLKMC.
//    */
//   class Base {
//   public:
//     // Copy this object
//     virtual std::unique_ptr<Base> clone() const = 0;

//     /**
//      * @brief Get this potentials cut-off radius that the neighbour lists should be configured for.
//      */
//     virtual floating rcut() const = 0;

//     /**
//      * @brief Compute the energy, assumes the neighbour list are ready.
//      *
//      * Ignores contribution from frozen atoms.
//      *
//      * A conforming potential must only use the positions of the atoms as defined in the nighbour
//      * list which *may* be different from those in the SimCell
//      */
//     virtual floating energy(SimCell const &, neighbour::List &, std::size_t threads) = 0;

//     /**
//      * @brief Compute gradient, assumes the neighbour list are ready.
//      *
//      * Force on frozen atoms will be zero.
//      *
//      * A conforming potential must only use the positions of the atoms as defined in the nighbour
//      * list which *may* be different from those in the SimCell
//      *
//      * @return May optionally return curvature along the minimum mode.
//      */
//     virtual std::optional<floating> gradient(SimCell &, neighbour::List &, std::size_t threads) = 0;

//     /**
//      * @brief Compute hessian matrix, assumes the neighbour list are ready.
//      *
//      * The resulting hessian will be m by m and only include contributions from the m active atoms.
//      *
//      * A conforming potential must only use the positions of the atoms as defined in the nighbour
//      * list which *may* be different from those in the SimCell
//      */
//     virtual void hessian(SimCell &, neighbour::List &, std::size_t threads) = 0;

//     /**
//      * @brief Call parent destructor.
//      */
//     virtual ~Base() {}

//   protected:
//     /**
//      * @brief Protected constructor as this is an interface class.
//      */
//     constexpr Base() noexcept = default;
//   };

// }  // namespace otf::potentials