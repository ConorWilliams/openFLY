
#include <array>
#include <string_view>

/**
 * \file data.hpp
 *
 * @brief Arrays of atomic data.
 */

namespace fly::data {

  /**
   * @brief Named struct for storing atomic data.
   */
  struct sym_mass_num {
    std::string_view symbol;  ///< The atomic symbol.
    double mass;              ///< The average atomic mass.
    int number;               ///< The atomic/proton number.
  };

  /**
   * @brief Atomic data for the first 118 elements.
   */
  constexpr std::array atoms = {
      sym_mass_num{"H", 1.008, 1},           sym_mass_num{"He", 4.0026022, 2},
      sym_mass_num{"Li", 6.94, 3},           sym_mass_num{"Be", 9.01218315, 4},
      sym_mass_num{"B", 10.81, 5},           sym_mass_num{"C", 12.011, 6},
      sym_mass_num{"N", 14.007, 7},          sym_mass_num{"O", 15.999, 8},
      sym_mass_num{"F", 18.9984031636, 9},   sym_mass_num{"Ne", 20.17976, 10},
      sym_mass_num{"Na", 22.989769282, 11},  sym_mass_num{"Mg", 24.305, 12},
      sym_mass_num{"Al", 26.98153857, 13},   sym_mass_num{"Si", 28.085, 14},
      sym_mass_num{"P", 30.9737619985, 15},  sym_mass_num{"S", 32.06, 16},
      sym_mass_num{"Cl", 35.45, 17},         sym_mass_num{"Ar", 39.9481, 18},
      sym_mass_num{"K", 39.09831, 19},       sym_mass_num{"Ca", 40.0784, 20},
      sym_mass_num{"Sc", 44.9559085, 21},    sym_mass_num{"Ti", 47.8671, 22},
      sym_mass_num{"V", 50.94151, 23},       sym_mass_num{"Cr", 51.99616, 24},
      sym_mass_num{"Mn", 54.9380443, 25},    sym_mass_num{"Fe", 55.8452, 26},
      sym_mass_num{"Co", 58.9331944, 27},    sym_mass_num{"Ni", 58.69344, 28},
      sym_mass_num{"Cu", 63.5463, 29},       sym_mass_num{"Zn", 65.382, 30},
      sym_mass_num{"Ga", 69.7231, 31},       sym_mass_num{"Ge", 72.6308, 32},
      sym_mass_num{"As", 74.9215956, 33},    sym_mass_num{"Se", 78.9718, 34},
      sym_mass_num{"Br", 79.904, 35},        sym_mass_num{"Kr", 83.7982, 36},
      sym_mass_num{"Rb", 85.46783, 37},      sym_mass_num{"Sr", 87.621, 38},
      sym_mass_num{"Y", 88.905842, 39},      sym_mass_num{"Zr", 91.2242, 40},
      sym_mass_num{"Nb", 92.906372, 41},     sym_mass_num{"Mo", 95.951, 42},
      sym_mass_num{"Tc", 98.0, 43},          sym_mass_num{"Ru", 101.072, 44},
      sym_mass_num{"Rh", 102.905502, 45},    sym_mass_num{"Pd", 106.421, 46},
      sym_mass_num{"Ag", 107.86822, 47},     sym_mass_num{"Cd", 112.4144, 48},
      sym_mass_num{"In", 114.8181, 49},      sym_mass_num{"Sn", 118.7107, 50},
      sym_mass_num{"Sb", 121.7601, 51},      sym_mass_num{"Te", 127.603, 52},
      sym_mass_num{"I", 126.904473, 53},     sym_mass_num{"Xe", 131.2936, 54},
      sym_mass_num{"Cs", 132.905451966, 55}, sym_mass_num{"Ba", 137.3277, 56},
      sym_mass_num{"La", 138.905477, 57},    sym_mass_num{"Ce", 140.1161, 58},
      sym_mass_num{"Pr", 140.907662, 59},    sym_mass_num{"Nd", 144.2423, 60},
      sym_mass_num{"Pm", 145.0, 61},         sym_mass_num{"Sm", 150.362, 62},
      sym_mass_num{"Eu", 151.9641, 63},      sym_mass_num{"Gd", 157.253, 64},
      sym_mass_num{"Tb", 158.925352, 65},    sym_mass_num{"Dy", 162.5001, 66},
      sym_mass_num{"Ho", 164.930332, 67},    sym_mass_num{"Er", 167.2593, 68},
      sym_mass_num{"Tm", 168.934222, 69},    sym_mass_num{"Yb", 173.0451, 70},
      sym_mass_num{"Lu", 174.96681, 71},     sym_mass_num{"Hf", 178.492, 72},
      sym_mass_num{"Ta", 180.947882, 73},    sym_mass_num{"W", 183.841, 74},
      sym_mass_num{"Re", 186.2071, 75},      sym_mass_num{"Os", 190.233, 76},
      sym_mass_num{"Ir", 192.2173, 77},      sym_mass_num{"Pt", 195.0849, 78},
      sym_mass_num{"Au", 196.9665695, 79},   sym_mass_num{"Hg", 200.5923, 80},
      sym_mass_num{"Tl", 204.38, 81},        sym_mass_num{"Pb", 207.21, 82},
      sym_mass_num{"Bi", 208.980401, 83},    sym_mass_num{"Po", 209.0, 84},
      sym_mass_num{"At", 210.0, 85},         sym_mass_num{"Rn", 222.0, 86},
      sym_mass_num{"Fr", 223.0, 87},         sym_mass_num{"Ra", 226.0, 88},
      sym_mass_num{"Ac", 227.0, 89},         sym_mass_num{"Th", 232.03774, 90},
      sym_mass_num{"Pa", 231.035882, 91},    sym_mass_num{"U", 238.028913, 92},
      sym_mass_num{"Np", 237.0, 93},         sym_mass_num{"Pu", 244.0, 94},
      sym_mass_num{"Am", 243.0, 95},         sym_mass_num{"Cm", 247.0, 96},
      sym_mass_num{"Bk", 247.0, 97},         sym_mass_num{"Cf", 251.0, 98},
      sym_mass_num{"Es", 252.0, 99},         sym_mass_num{"Fm", 257.0, 100},
      sym_mass_num{"Md", 258.0, 101},        sym_mass_num{"No", 259.0, 102},
      sym_mass_num{"Lr", 266.0, 103},        sym_mass_num{"Rf", 267.0, 104},
      sym_mass_num{"Db", 268.0, 105},        sym_mass_num{"Sg", 269.0, 106},
      sym_mass_num{"Bh", 270.0, 107},        sym_mass_num{"Hs", 269.0, 108},
      sym_mass_num{"Mt", 278.0, 109},        sym_mass_num{"Ds", 281.0, 110},
      sym_mass_num{"Rg", 282.0, 111},        sym_mass_num{"Cn", 285.0, 112},
      sym_mass_num{"Nh", 286.0, 113},        sym_mass_num{"Fl", 289.0, 114},
      sym_mass_num{"Mc", 289.0, 115},        sym_mass_num{"Lv", 293.0, 116},
      sym_mass_num{"Ts", 294.0, 117},        sym_mass_num{"Og", 294.0, 118},
  };
}  // namespace fly::data