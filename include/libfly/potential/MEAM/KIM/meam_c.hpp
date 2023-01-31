//
// meam_c.hpp
//
// LGPL Version 2.1 HEADER START
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
//
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301  USA
//
// LGPL Version 2.1 HEADER END
//

//
// Copyright (c) 2020--2021, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//
// Brief: This file is adapted from the LAMMPS software package
//        `lammps/src/USER-MEAMC/meam_dens_final.cpp`
//        `lammps/src/USER-MEAMC/meam_dens_init.cpp`
//        `lammps/src/USER-MEAMC/meam_force.cpp`
//        `lammps/src/USER-MEAMC/meam_funcs.cpp`
//        `lammps/src/USER-MEAMC/meam_impl.cpp`
//        `lammps/src/USER-MEAMC/meam_setup_done.cpp`
//        `lammps/src/USER-MEAMC/meam_setup_global.cpp`
//        `lammps/src/USER-MEAMC/meam_setup_param.cpp`
//        `lammps/src/USER-MEAMC/meam.h`
//        `lammps/src/USER-MEAMC/pair_meamc.cpp`
//        `lammps/src/USER-MEAMC/pair_meamc.h`
//        and it is rewritten and updated for the KIM-API by
//        Yaser Afshar
//

#ifndef MEAM_C_HPP
#define MEAM_C_HPP

#pragma GCC system_header

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "helper.hpp"
#include "zbl.hpp"

namespace openKIM {

  /*!
   * \brief Lattice types
   *
   */
  enum class Lattice : int { FCC, BCC, HCP, DIM, DIA, DIA3, B1, C11, L12, B2, CH4, LIN, ZIG, TRI };

  /*! \class MEAMC
   * \brief Modified Embedded Atom Method
   *
   */
  class MEAMC {
  public:
    /*!
     * \brief Construct a new MEAMC object
     *
     */
    MEAMC() = default;

    /*!
     * \brief Destroy the MEAMC object
     *
     */
    ~MEAMC() = default;

  public:
    /*!
     * \brief MEAM cutoff function, \c Eq.4.11e.
     *
     * The cutoff function describes the smooth, gradually decreasing distance
     * effect on the interactions between the atoms. Normally, the MEAM cut-off
     * function has the \f$ (1-x) \f$ term to the fourth power, and in some
     * cases to the sixth power.
     *
     * \c form=0
       \f[
        f_c(x)=\left\{\begin{matrix}
        1 && x\geqslant 1\\
        \left [1-\left(1-x\right)^4\right]^2 && 0<x<1\\
        0 && x\leqslant 0
        \end{matrix}\right.
       \f]
     *
     * \c form=1
       \f[
        f_c(x)=\left\{\begin{matrix}
        1 && x\geqslant 1\\
        \left [1-\left(1-x\right)^6\right]^2 && 0<x<1\\
        0 && x\leqslant 0
        \end{matrix}\right.
       \f]
     *
     * \param x input x
     * \param df derivative of the radial cutoff function
     * \param form different MEAM cut-off function form {one of 0, 1}
     * \return double
     */
    static inline double RadialCutoff(double const x, int const form);
    static inline double RadialCutoff(double const x, double &df, int const form);

    /*!
     * \brief Derivative of Cikj w.r.t. rik and rjk, Cikj,ik, and Cikj,jk
     *        \c Eq.4.17b & \c Eq.4.17c, respectively.
     *
     * \param rij2 squared distance between I, J
     * \param rik2 squared distance between I, K
     * \param rjk2 squared distance between J, K
     * \param dcijk_drik \f$ \frac{\partial C_{ijk}}{\partial r_{ik}} \f$
     * \param dcijk_drjk \f$ \frac{\partial C_{ijk}}{\partial r_{jk}} \f$
     */
    static inline void DCijkDRikDRjk(double const rij2,
                                     double const rik2,
                                     double const rjk2,
                                     double &dcijk_drik,
                                     double &dcijk_drjk);

    /*!
     * \brief Get the shape factors for various configurations
     *
     * \param lat lattice type
     * \param stheta
     * \param ctheta
     * \param shape_factors
     */
    static void GetShapeFactors(Lattice const &lat,
                                double const stheta,
                                double const ctheta,
                                VectorOfSizeDIM &shape_factors);

    /*!
     * \brief Convert lattice spec to string
     *
     * \param lat lattice
     * \return std::string
     */
    static std::string LatticeToString(Lattice const &lat);

    /*!
     * \brief Process and parse the meam/c parameter file
     *
     * \param library_file_pointer FILE pointer to the opened file
     * \param max_line_size maximum line size
     * \return int
     */
    int ProcessLibraryFile(std::FILE *const library_file_pointer,
                           int const max_line_size,
                           std::vector<std::string> const &element_name);

    /*!
     * \brief Process and parse the meam/c parameter file
     *
     * \param parameter_file_pointer FILE pointer to the opened file
     * \param max_line_size maximum line size
     * \return int
     */
    int ProcessParameterFile(std::FILE *const parameter_file_pointer, int const max_line_size);

    /*!
     * \brief Get the phi and its derivative with respect to rij
     *
     * phi & dphi are defined in \c Eq.4.10b
     *
     * \param species_i I type
     * \param species_j J type
     * \param rij distance
     * \param dphi derivative of phi
     * \return double
     */
    inline double GetPhiAndDerivative(int const species_i,
                                      int const species_j,
                                      double const rij,
                                      double &dphi) const;

    /*!
     * \brief Convert units of the parameters
     *
     * \param convert_length_factor length unit conversion factor
     * \param convert_energy_factor energy unit conversion factor
     */
    void ConvertUnit(double const convert_length_factor, double const convert_energy_factor);

    /*!
     * \brief Complete the set up
     *
     * \param max_cutoff
     */
    void CompleteSetup(double *max_cutoff);

    /*!
     * \brief Grow element arrays and initialize them
     *
     */
    void ResizeElementArrays();

    /*!
     * \brief Grow local arrays if necessary and initialize them to zero
     *
     * \param nall Total number of atoms contributing and non-contributing
     */
    void ResizeDenistyArrays(std::size_t const nall);

    /*!
     * \brief Grow local arrays \c scrfcn_, and \c dscrfcn_ if necessary
     *
     * \param n_neigh Total number of neighbors of all contributing atoms
     */
    void ResizeScreeningArrays(std::size_t const n_neigh);

    /*!
     * \brief First stage in MEAMC density calculation
     *
     * \param i I atom
     * \param number_of_neighbors # of J neighbors for each I atom
     * \param neighbors_of_particle ptr to 1st J int value of each I atom
     * \param offset Offset to the half list number of neighbors
     * \param coordinates Atoms coordinates
     * \param particle_species_codes Particle species code
     * \param particle_contributing Particle contirubuting flag list
     */
    void InitializeDensityCalculation(int const i,
                                      int const number_of_neighbors,
                                      int const *const neighbors_of_particle,
                                      int &offset,
                                      const VectorOfSizeDIM *const coordinates,
                                      int const *const particle_species_codes,
                                      int const *const particle_contributing);

    /*!
     * \brief Final stage in MEAMC density calculation
     *
     * \param i I atom
     * \param species_i I type
     * \param embedding Energy
     */
    void FinalizeDensityCalculation(int const i, int const species_i, double &embedding);

    /*!
     * \brief Compute the atomic electron densities and derivatives \c Eq.4.8
     *
     * \param species_i i atom type
     * \param species_j j atom type
     * \param rij r
     * \param rhoa0_i atomic electron density \c Eq.4.8
     * \param drhoa0_i atomic electron density derivative
     * \param rhoa1_i
     * \param drhoa1_i
     * \param rhoa2_i
     * \param drhoa2_i
     * \param rhoa3_i
     * \param drhoa3_i
     * \param rhoa0_j
     * \param drhoa0_j
     * \param rhoa1_j
     * \param drhoa1_j
     * \param rhoa2_j
     * \param drhoa2_j
     * \param rhoa3_j
     * \param drhoa3_j
     */
    void ComputeAtomicElectronDensities(int const species_i,
                                        int const species_j,
                                        double const rij,
                                        double &rhoa0_i,
                                        double &drhoa0_i,
                                        double &rhoa1_i,
                                        double &drhoa1_i,
                                        double &rhoa2_i,
                                        double &drhoa2_i,
                                        double &rhoa3_i,
                                        double &drhoa3_i,
                                        double &rhoa0_j,
                                        double &drhoa0_j,
                                        double &rhoa1_j,
                                        double &drhoa1_j,
                                        double &rhoa2_j,
                                        double &drhoa2_j,
                                        double &rhoa3_j,
                                        double &drhoa3_j);

  protected:
    /*!
     * \brief Convert lattice spec to Lattice only use single-element lattices
     *        if single=true
     *
     * \param str lattice spec string
     * \param single single-element or not
     * \param lat lattice type
     * \return false and set lat on success
     * \return true on failure
     */
    static bool StringToLattice(std::string const &str, bool const single, Lattice &lat);

    /*!
     * \brief Compute the Sijk \c Eq.4.11c
     *        \f$ S_{ijk} = f_c\left(\frac{(C-Cmin)}{(Camx-Cmin)} \right) \f$
     *
     * \param cijk C for the i-j-k triplet
     * \param i element type i
     * \param j element type j
     * \param k element type k
     * \return double
     */
    inline double Sijk(double const cijk, int const i, int const j, int const k) const;

    /*!
     * \brief Derivative of Cikj w.r.t. rij, Cikj,ij \c Eq.4.17a
     *
     * \param rij2 squared distance between I, J
     * \param rik2 squared distance between I, K
     * \param rjk2 squared distance between J, K
     * \return double \f$ \frac{\partial C_{ijk}}{\partial r_{ij}} \f$
     */
    static inline double DCijkDRij(double const rij2, double const rik2, double const rjk2);

    /*!
     * \brief Compute Rose energy function
     *
     * This function gives the energy of the reference state as a function of
     * interatomic spacing. The form of this function is:
     *
     * \f$ a^* = \alpha * (\frac{r}{r_e} - 1) \f$
     *
     * \c form=0
       \f[
       \text{form}=0\rightarrow E=\left\{\begin{matrix}
       -E_c\left(1+a^*+ {a_\text{repuls}}~ \frac{{a^*}^3}{r/re} \right)exp(-a^*)
       & a^* < 0\\ -E_c\left(1+a^*+ {a_\text{attrac}}~ \frac{{a^*}^3}{r/re}
       \right)exp(-a^*) & a^* \geqslant 0 \end{matrix}\right.
       \f]
     *
     * \c form=1
       \f[
       \text{form}=1\rightarrow E=-E_c\left(1+a^*+ \left(-{a_\text{attrac}}~+
       \frac{{a_\text{repuls}}}{r}\right){a^*}^3 \right)exp(-a^*)
       \f]
     * \c form=2
       \f[
       \text{form}=2\rightarrow E=\left\{\begin{matrix}
       -E_c\left(1+a^*+ {a_\text{repuls}}~ {a^*}^3 \right)exp(-a^*) & a^* < 0\\
       -E_c\left(1+a^*+ {a_\text{attrac}}~ {a^*}^3 \right)exp(-a^*) & a^*
       \geqslant 0 \end{matrix}\right.
       \f]
     *
     * \param r interatomic spacing
     * \param re equilibrium distance between I and J in the reference structure
     * \param alpha alpha parameter for pair potential between I and J (can be
     *              be computed from bulk modulus of reference structure
     * \param Ec cohesive energy of reference structure for I-J mixture
     * \param repuls additional cubic repulsive term in the Rose energy I-J pair
     *               potential
     * \param attrac additional cubic attraction term in the Rose energy I-J pair
     *               potential
     * \param form different Rose energy function form {one of 0, 1, or 2}
     *
     * \return double computed Rose energy function value
     */
    static double Rose(double const r,
                       double const re,
                       double const alpha,
                       double const Ec,
                       double const repuls,
                       double const attrac,
                       int const form);

    /*!
     * \brief Number of nearest neighbors in the reference structure.
     *
     * Get the Zij object (Number of nearest neighbors in the
     * reference structure)
     *
     * \param lat lattice type
     * \return int
     */
    static double NumNearestNeighborsInReferenceStructure(Lattice const &lat);

    /*!
     * \brief Number of second neighbors for the reference structure.
     *
     * numscr = number of atoms that screen the 2NN bond
     *
     * \param lat Lattice type
     * \param cmin Cmin screening parameter
     * \param cmax Cmax screening parameter
     * \param stheta \c sin of theta angle between three atoms
     * \param distance_ratio distance ratio R1/R2
     * \param second_neighbor_screening second neighbor screening function
     * \return int number of second neighbors for the reference structure
     */
    double NumSecondNearestNeighborsInReferenceStructure(Lattice const &lat,
                                                         double const cmin,
                                                         double const cmax,
                                                         double const stheta,
                                                         double &distance_ratio,
                                                         double &second_neighbor_screening);

    /*!
     * \brief Compute the electron density
     *
     * Available forms of the function G(Gamma):
     *
     * ibar =  0 => G = sqrt(1+gamma)
     * ibar =  1 => G = exp(gamma/2)
     * ibar =  2 => not implemented
     * ibar =  3 => G = 2/(1+exp(-gamma))
     * ibar =  4 => G = sqrt(1+gamma)
     * ibar = -5 => G = +-sqrt(abs(1+gamma))
     *
     * \param gamma computed from \c Eq.4.4 or \c Eq.4.6
     * \param ibar selects the form of the function G(Gamma)
     * \param dg_gamma derivative of the gamma
     * \return double
     */
    double GGamma(double const gamma, int const ibar) const;
    double GGamma(double const gamma, int const ibar, double &dg_gamma) const;

    /*!
     * \brief Compute embedding function F(rhobar) and derivative F'(rhobar),
     *        \c Eq.4.2 & \c Eq.4.42 respectively.
     *
     *
     * \param A adjustable parameter
     * \param Ec cohesive energy of reference structure for I-J mixture
     * \param rhobar background electron density for a reference structure
     * \param embedding_df derivative of the Embedding function
     * \return double computed embedding function value
     */
    double Embedding(double const A, double const Ec, double const rhobar) const;
    double Embedding(double const A, double const Ec, double const rhobar, double &embedding_df) const;

    /*!
     * \brief Compute MEAMC pair potential for distance r, element types a and b
     *
     * \param r distance between atoms i and j (element types a and b)
     * \param a element type a
     * \param b element type b
     * \return double
     */
    double ComputePhi(double const r, int const a, int const b);

    /*!
     * \brief Compute 2NN series terms for phi
     *
     * To avoid nan values of phir_ due to rapid decrease of
     * b2nn^n or/and argument of ComputePhi,
     * i.e. r*distance_ratio^n, in some cases (3NN dia with low Cmin value)
     *
     * \param second_neighbor_screening second neighbor screening function
     * \param Z1
     * \param Z2
     * \param r distance between atoms i and j (element types a and b)
     * \param a element type a
     * \param b element type b
     * \param distance_ratio distance ratio
     * \return double
     */
    double ComputePhiSeries(double const second_neighbor_screening,
                            double const Z1,
                            double const Z2,
                            double const r,
                            int const a,
                            int const b,
                            double const distance_ratio);

    /*!
     * \brief A sanity check on index parameters
     *
     * \param num
     * \param lim
     * \param nidx
     * \param idx
     * \param ierr error number
     */
    void CheckIndex(int const num, int const lim, int const nidx, int const *idx, int *ierr);

    /*!
     * \brief Set up the parameters
     *
     * Set up the parameters for any of the 22 keywords:
     * \c  0-->Ec
     * \c  1-->alpha
     * \c  2-->rho0
     * \c  3-->delta
     * \c  4-->lattice
     * \c  5-->attrac
     * \c  6-->repuls
     * \c  7-->nn2
     * \c  8-->Cmin
     * \c  9-->Cmax
     * \c 10-->rc
     * \c 11-->delr
     * \c 12-->augt1
     * \c 13-->gsmooth_factor
     * \c 14-->element_re
     * \c 15-->ialloy
     * \c 16-->mixture_ref_t
     * \c 17-->erose_form
     * \c 18-->element_zbl
     * \c 19-->emb_lin_neg
     * \c 20-->bkgd_dyn
     * \c 21-->theta
     * \c 22-->fcut_form
     *
     * \param which corresponds to the index of the "keyword" array
     * \param value the parameter value to set
     * \param nindex number of indices
     * \param index array of indices
     * \param errorflag The returned errorflag has the following meanings:
     *                  0 = no error
     *                  1 = "which" out of range / invalid keyword
     *                  2 = not enough indices given
     *                  3 = an element index is out of range
     */
    void SetParameter(int const which,
                      double const value,
                      int const nindex,
                      int const *index,
                      int *errorflag);

    /*!
     * \brief Fill off-diagonal alloy parameters
     *
     */
    void FillOffDiagonalAlloyParameters();

    /*!
     * \brief Compute screening function and its derivative
     *
     * \param i I atom
     * \param number_of_neighbors # of J neighbors for each I atom
     * \param neighbors_of_particle pointer to 1st J int value of each I atom
     * \param offset offset to the half list number of neighbors
     * \param coordinates atoms coordinates
     * \param particle_species_codes particle species code
     * \param particle_contributing particle contirubuting flag list
     */
    void ComputeScreeningAndDerivative(int const i,
                                       int const number_of_neighbors,
                                       int const *const neighbors_of_particle,
                                       int const offset,
                                       const VectorOfSizeDIM *const coordinates,
                                       int const *const particle_species_codes,
                                       int const *const particle_contributing);

    /*!
     * \brief Compute intermediate density terms
     *
     * \param i I atom
     * \param number_of_neighbors # of J neighbors for each I atom
     * \param neighbors_of_particle pointer to 1st J int value of each I atom
     * \param offset offset to the half list number of neighbors
     * \param coordinates atoms coordinates
     * \param particle_species_codes particle species code
     * \param particle_contributing particle contirubuting flag list
     */
    void ComputeIntermediateDensityTerms(int const i,
                                         int const number_of_neighbors,
                                         int const *const neighbors_of_particle,
                                         int &offset,
                                         const VectorOfSizeDIM *const coordinates,
                                         int const *const particle_species_codes,
                                         int const *const particle_contributing);

    /*!
     * \brief Grow local pair potential arrays
     *
     */
    void ResizePairPotentialArrays();

    /*!
     * \brief Compute MEAMC pair potential for each pair of element types
     *
     */
    void ComputePairPotential();

    /*!
     * \brief Compute the composition dependent electron density scaling Eq.4.5
     *
     */
    void ComputeCompositionDependentDensityScaling();

    /*!
     * \brief Get the densref object
     *
     * Calculate density functions, assuming reference configuration
     *
     * \param r distance between atoms i and j (element types a and b)
     * \param a element type a
     * \param b element type b
     * \param rho0_a
     * \param rho1_a
     * \param rho2_a
     * \param rho3_a
     * \param rho0_b
     * \param rho1_b
     * \param rho2_b
     * \param rho3_b
     */
    void ComputeReferenceConfigurationDensity(double const r,
                                              int const a,
                                              int const b,
                                              double *rho0_a,
                                              double *rho1_a,
                                              double *rho2_a,
                                              double *rho3_a,
                                              double *rho0_b,
                                              double *rho1_b,
                                              double *rho2_b,
                                              double *rho3_b);
    /*!
     * \brief interpolation
     *
     * \param ind
     */
    void SplineInterpolate(int const ind);

  public:
    /*!
     * \brief The integer flag indicates whether to augment the \c t1 parameter
     *        by 3/5 * t3 to account for old vs. new meam formulations
     *
     *        0 = don't augment t1
     *        1 = augment t1
     */
    int augment_t1_{1};

    /*!
     * \brief The integer flag indicates whether to use an alternative averaging
     *        rule for the \c t parameters
     *
     * For comparison with the \c DYNAMO MEAMC code
     * \c ialloy_=0, standard averaging (matches ialloy=0 in \c DYNAMO code)
     * \c ialloy_=1, alternative averaging (matches ialloy=1 in \c DYNAMO code)
     * \c ialloy_=2, no averaging of t (use single-element values)
     */
    int ialloy_{0};

    /*!
     * \brief The integer flag indicates whether to use the mixture averaging
     *        rule for the \c t parameters
     *
     *        0 -> don't use mixture averaging for \c t in the reference density
     *        1 -> use mixture averaging for \c t in the reference density
     */
    int mixing_rule_compute_t_{0};

    /*!
     * \brief The integer value used to select the form of the Rose energy
     *        function. \sa Rose
     *
     *        0 -> form 0
     *        1 -> form 1
     *        2 -> form 2
     */
    int rose_function_form_{0};

    /*!
     * \brief The integer value used to select the form of the MEAM cut-off
     *        function. \sa RadialCutoff
     *
     *        0 -> has the \f$ (1-x) \f$ term to the fourth power
     *        1 -> has the \f$ (1-x) \f$ term to the sixth power
     *
     */
    int fcut_function_form_{0};

    /*!
     * \brief The integer value used to select embedding function for negative
     *        densities
     *
     *        0 -> F(rho)=0
     *        1 -> F(rho) = -asub*esub*rho (linear in rho, matches \c DYNAMO code)
     */
    int use_rhobar_linear_embedding_function_{0};

    /*!
     * \brief The integer value used to select background density formula
     *
     *        0 -> rho_ref_meam(a) (as in the reference structure)
     *        1 -> rho0_meam(a) * element_rho0_(a) (matches \c DYNAMO code)
     */
    int dynamo_reference_density_{0};

    /*! Pair function discretization parameters / spline coeff array parameters */
    int nr_{1000};

    /*! Cutoff radius */
    double cutoff_radius_{4.0};

    /*!
     * \brief Length of the smoothing distance for cutoff function
     *
     * It controls the distance over which the radial cutoff is
     * smoothed from 1 to 0 near r = cutoff_radius_.
     */
    double delr_{0.1};

    /*! Pair function discretization parameters / spline coeff array parameters */
    double dr_;

    /*!
     * \brief Factor determining the length of the G-function smoothing
     *        factor determining the length of the G-function smoothing
     *
     *        99.0 -> short smoothing region, sharp step
     *        0.5  -> long smoothing region, smooth step
     */
    double gsmooth_factor_{99.0};

    /*! Selection parameter for Gamma function for different element */
    std::vector<int> element_ibar_;

    /*! Atomic number of element */
    std::vector<double> element_atomic_number_;

    /*! Electron density constants */
    std::vector<double> element_beta0_;
    /*! Electron density constants */
    std::vector<double> element_beta1_;
    /*! Electron density constants */
    std::vector<double> element_beta2_;
    /*! Electron density constants */
    std::vector<double> element_beta3_;

    /*! Adjustable parameter */
    std::vector<double> element_A_;

    /*! The average weighting factors Eq.4.9 */
    std::vector<double> element_t1_;
    /*! The average weighting factors Eq.4.9 */
    std::vector<double> element_t2_;
    /*! The average weighting factors Eq.4.9 */
    std::vector<double> element_t3_;

    /*!
     * \brief Relative density for element I.
     *
     * Element-dependent density scaling.
     */
    std::vector<double> element_rho0_;

    /*!
     * \brief lattice structure of I-J reference structure:
     *        \c fcc  -> face centered cubic
     *        \c bcc  -> body centered cubic
     *        \c hcp  -> hexagonal close-packed
     *        \c dim  -> dimer
     *        \c dia  -> diamond (interlaced fcc for alloy)
     *        \c dia3 -> diamond structure with primary 1NN and secondary 3NN
     *                   interaction
     *        \c b1   -> rock salt (NaCl structure)
     *        \c c11  -> MoSi2 structure
     *        \c l12  -> Cu3Au structure (lower case L, followed by 12)
     *        \c b2   -> CsCl structure (interpenetrating simple cubic)
     *        \c ch4  -> methane-like structure, only for binary system
     *        \c lin  -> linear structure (180 degree angle)
     *        \c zig  -> zigzag structure with a uniform angle
     *        \c tri  -> H2O-like structure that has an angle
     */
    Array2D<Lattice> element_lattice_;

    /*!
     * \brief Turn on second-nearest neighbor MEAMC formulation for I-J pair.
     *        (1 if second nearest neighbors are to be computed, else 0)
     *
     *        \c 0 -> second-nearest neighbor formulation off
     *        \c 1 -> second-nearest neighbor formulation on
     */
    Array2D<int> element_nn2_;

    /*!
     * \brief blend the MEAMC I-J pair potential with the ZBL potential for
     *        small atom separations (1 if Zbl potential for small r to be used,
     *        else 0)
     */
    Array2D<int> element_zbl_;

    /*!
     * \brief Alpha parameter for pair potential between I and J
     *       (can be computed from bulk modulus of reference structure)
     */
    Array2D<double> element_alpha_;

    /*!
     * \brief Equilibrium distance between I and J in the reference
     *        structure.
     *
     * Nearest-neighbor distance in the single-element reference structure
     */
    Array2D<double> element_re_;

    /*! Cohesive energy */
    Array2D<double> element_Ec_;

    /*!
     * \brief Heat of formation for I-J alloy.
     *
     * If element_Ec_IJ is zero on input, then it sets
     * element_Ec_IJ = (element_Ec_II + element_Ec_JJ)/2 - element_delta_
     */
    Array2D<double> element_delta_;

    /*! Additional cubic attraction term in Rose energy I-J pair potential */
    Array2D<double> element_attrac_;

    /*! Additional cubic repulsive term in Rose energy I-J pair potential */
    Array2D<double> element_repuls_;

    /*!
     * \brief angle between three atoms in line, zigzag, and trimer reference
     *        structures in degrees. (default = 180)
     */
    Array2D<double> element_theta_;

    /*!
     * \brief sin(theta/2) in radian used in line, zigzag, and trimer reference
     *        structures
     *
     * sin(theta/2) in radian used in line, zigzag, and trimer reference
     * structures theta = angle between three atoms in line, zigzag, and trimer
     * reference structures
     */
    Array2D<double> element_stheta_;

    /*!
     * \brief cos(theta/2) in radian used in line, zigzag, and trimer reference
     *        structures
     *
     * cos(theta/2) in radian used in line, zigzag, and trimer reference
     * structures theta = angle between three atoms in line, zigzag, and trimer
     * reference structures
     */
    Array2D<double> element_ctheta_;

    /*! factor giving maximum boundary of screen function ellipse */
    Array2D<double> element_ebound_;

    /*! Screening function */
    std::vector<double> scrfcn_;

    /*! Derivative of the screening function */
    std::vector<double> dscrfcn_;

    std::vector<double> rho_;
    std::vector<double> frhop_;

    std::vector<double> rho0_;
    std::vector<double> rho1_;
    std::vector<double> rho2_;
    std::vector<double> rho3_;

    std::vector<double> gamma_;
    std::vector<double> dgamma1_;
    std::vector<double> dgamma2_;
    std::vector<double> dgamma3_;

    std::vector<double> arho2b_;
    Array2D<double> arho1_;
    Array2D<double> arho2_;
    Array2D<double> arho3_;
    Array2D<double> arho3b_;
    Array2D<double> t_ave_;
    Array2D<double> tsq_ave_;

    /*!
     * \brief Cmin screening parameter when I-J pair is screened by K (I<=J)
     *        (default = 2.0)
     *
     * The min values in screening cutoff
     */
    Array3D<double> element_Cmin_;

    /*!
     * \brief Cmax screening parameter when I-J pair is screened by K (I<=J)
     *        (default = 2.8)
     *
     * The max values in screening cutoff
     */
    Array3D<double> element_Cmax_;

  private:
    /*! The number of element types */
    int number_of_element_types_{0};

    /*! force cutoff squared */
    double cutoff_radius_squared_;

    /*! Heat of formation for alloys */

    /*! Composition dependent electron density scaling Eq.4.5 */
    std::vector<double> element_ref_rho_;

    /*! index number of pair (similar to Voigt notation; ij = ji) */
    Array2D<int> element_pair_index_;

    /*! pair potential function array */
    Array2D<double> phir_;
    /*! spline coeffs */
    Array2D<double> phirar1_;
    /*! spline coeffs */
    Array2D<double> phirar2_;
    /*! spline coeffs */
    Array2D<double> phirar3_;
    /*! spline coeffs */
    Array2D<double> phirar4_;
    /*! spline coeffs */
    Array2D<double> phirar5_;
    /*! spline coeffs */
    Array2D<double> phirar6_;

    /*! ZBL object */
    std::unique_ptr<ZBL> zbl_;
  };

  // Eq.4.11e
  inline double MEAMC::RadialCutoff(double const x, int const form) {
    if (x >= 1.0) {
      return 1.0;
    }
    if (x <= 0.0) {
      return 0.0;
    }
    double a = 1.0 - x;
    if (form == 0) {
      a *= a;
      a *= a;
      a = 1.0 - a;
      return a * a;
    }
    // form == 1
    double const a3 = a * a * a;
    a = a3 * a3;
    a = 1.0 - a;
    return a * a;
  }

  // Eq.4.11e
  inline double MEAMC::RadialCutoff(double const x, double &df, int const form) {
    if (x >= 1.0) {
      df = 0.0;
      return 1.0;
    }
    if (x <= 0.0) {
      df = 0.0;
      return 0.0;
    }
    double const a = 1.0 - x;
    if (form == 0) {
      double const a3 = a * a * a;
      double const a4 = a * a3;
      double const a1m4 = 1.0 - a4;
      df = 8 * a1m4 * a3;
      return a1m4 * a1m4;
    }
    // form == 1
    double const a2 = a * a;
    double const a5 = a * a2 * a2;
    double const a6 = a * a5;
    double const a1m6 = 1.0 - a6;
    df = 12 * a1m6 * a5;
    return a1m6 * a1m6;
  }

  inline void MEAMC::DCijkDRikDRjk(double const rij2,
                                   double const rik2,
                                   double const rjk2,
                                   double &dcijk_drik,
                                   double &dcijk_drjk) {
    double const rij4 = rij2 * rij2;
    double const rik4 = rik2 * rik2;
    double const rjk4 = rjk2 * rjk2;
    double const a = rik2 - rjk2;
    double denom = rij4 - a * a;
    denom *= denom;
    // Eq.4.17b
    dcijk_drik = 4 * rij2 * (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a);
    dcijk_drik /= denom;
    // Eq.4.17c
    dcijk_drjk = 4 * rij2 * (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a);
    dcijk_drjk /= denom;
  }

  // Eq.4.11c
  inline double MEAMC::Sijk(double const cijk, int const i, int const j, int const k) const {
    double const cmin = element_Cmin_(i, j, k);
    double const x = (cijk - cmin) / (element_Cmax_(i, j, k) - cmin);
    double sijk = RadialCutoff(x, fcut_function_form_);
    return sijk;
  }

  inline double MEAMC::DCijkDRij(double const rij2, double const rik2, double const rjk2) {
    double const rij4 = rij2 * rij2;
    double const a = rik2 - rjk2;
    double const b = rik2 + rjk2;
    double const asq = a * a;
    double denom = rij4 - asq;
    denom *= denom;
    // Eq.4.17a
    double dcijk_drij = -4 * (-2 * rij2 * asq + rij4 * b + asq * b);
    dcijk_drij /= denom;
    return dcijk_drij;
  }

  inline double MEAMC::GetPhiAndDerivative(int const species_i,
                                           int const species_j,
                                           double const rij,
                                           double &dphi) const {
    // Compute phi and phip
    int const ind = element_pair_index_(species_i, species_j);
    double pp = rij / dr_;
    int const kk = std::min(static_cast<int>(pp), nr_ - 2);
    pp -= kk;
    pp = std::min(pp, 1.0);
    // Derivative of Eq.4.10b with respect to rij
    dphi = (phirar6_(ind, kk) * pp + phirar5_(ind, kk)) * pp + phirar4_(ind, kk);
    // Eq.4.10b
    double phi
        = ((phirar3_(ind, kk) * pp + phirar2_(ind, kk)) * pp + phirar1_(ind, kk)) * pp + phir_(ind, kk);
    return phi;
  }

}  // namespace openKIM

#endif  // MEAM_C_HPP
