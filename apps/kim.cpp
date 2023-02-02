
#include <fmt/core.h>

#include <KIM_ChargeUnit.hpp>
#include <KIM_FunctionTypes.hpp>
#include <KIM_ModelRoutineName.hpp>
#include <KIM_SpeciesName.hpp>
#include <KIM_TemperatureUnit.hpp>
#include <KIM_TimeUnit.hpp>
#include <cstddef>
#include <exception>

#include "KIM_SimulatorHeaders.hpp"
#include "KIM_SupportedExtensions.hpp"

//

#include "libfly/env/catalogue.hpp"
#include "libfly/io/gsd.hpp"
#include "libfly/kinetic/basin.hpp"
#include "libfly/kinetic/skmc.hpp"
#include "libfly/kinetic/superbasin.hpp"
#include "libfly/minimise/LBFGS/lbfgs.hpp"
#include "libfly/neigh/list.hpp"
#include "libfly/neigh/sort.hpp"
#include "libfly/potential/generic.hpp"
#include "libfly/saddle//find.hpp"
#include "libfly/saddle/dimer.hpp"
#include "libfly/saddle/perturb.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/atom.hpp"
#include "libfly/system/box.hpp"
#include "libfly/system/boxes/orthorhombic.hpp"
#include "libfly/system/property.hpp"
#include "libfly/system/supercell.hpp"
#include "libfly/system/typemap.hpp"
#include "libfly/utility/core.hpp"
#include "libfly/utility/lattice.hpp"
#include "libfly/utility/random.hpp"

namespace fly::potential {

  class KIM_API {
  public:
    struct Options {
      std::string model_name = "";  // The exact name of the KIM model
    };

    explicit KIM_API(Options const& opt, system::TypeMap<> const& map) : m_opt(opt) {
      int units_accepted = false;
      int err = KIM::Model::Create(KIM::NUMBERING::zeroBased,
                                   KIM::LENGTH_UNIT::A,
                                   KIM::ENERGY_UNIT::eV,
                                   KIM::CHARGE_UNIT::unused,
                                   KIM::TEMPERATURE_UNIT::unused,
                                   KIM::TIME_UNIT::unused,
                                   m_opt.model_name,
                                   &units_accepted,
                                   &m_model);

      fly::verify(!err, "Kim model not found");
      fly::verify(units_accepted, "Kim model rejected our units");
      fly::verify(m_model, "Model pointer not assigned");

      // Check that we know about all required routines

      int num_names = 0;

      KIM::MODEL_ROUTINE_NAME::GetNumberOfModelRoutineNames(&num_names);

      for (int i = 0; i < num_names; ++i) {
        //
        KIM::ModelRoutineName name;

        if (KIM::MODEL_ROUTINE_NAME::GetModelRoutineName(i, &name)) {
          throw fly::error("Unable to get ModelRoutineName");
        }

        int present = 0;
        int required = 0;

        if (m_model->IsRoutinePresent(name, &present, &required)) {
          throw fly::error("Unable to get routine present/required.");
        }

        fmt::print("Model routine name \"{}\": present = {} and required = {}\n",
                   name.ToString(),
                   bool(present),
                   bool(required));

        if ((present == true) && (required == true)) {
          using namespace KIM::MODEL_ROUTINE_NAME;
          if (!((name == Create) || (name == ComputeArgumentsCreate) || (name == Compute) || (name == Refresh)
                || (name == ComputeArgumentsDestroy) || (name == Destroy))) {
            throw fly::error("Unknown Routine \"" + name.ToString() + "\" is required by model.");
          }
        }
      }

      // Print model units
      KIM::LengthUnit lengthUnit;
      KIM::EnergyUnit energyUnit;
      KIM::ChargeUnit chargeUnit;
      KIM::TemperatureUnit temperatureUnit;
      KIM::TimeUnit timeUnit;

      m_model->GetUnits(&lengthUnit, &energyUnit, &chargeUnit, &temperatureUnit, &timeUnit);

      fmt::print("Length unit is \"{}\"\n", lengthUnit.ToString());
      fmt::print("Energy unit is \"{}\"\n", energyUnit.ToString());
      fmt::print("Charge unit is \"{}\"\n", chargeUnit.ToString());
      fmt::print("Temperature unit is \"{}\"\n", temperatureUnit.ToString());
      fmt::print("Time unit is \"{}\"\n", timeUnit.ToString());

      // Check species

      for (Eigen::Index i = 0; i < map.num_types(); i++) {
        int supported;
        int model_code;

        KIM::SpeciesName name = std::string{
            map.get(i, tp_),
        };

        if (m_model->GetSpeciesSupportAndCode(name, &supported, &model_code) || !supported) {
          throw error(
              "Species \"{}\" unsupported by KIM potential \"{}\"", name.ToString(), m_opt.model_name);
        }

        verify(model_code == i, "Species {} uses KIM code {} != {}", name.ToString(), model_code, i);
      }

      // Set up compute request

      if (m_model->ComputeArgumentsCreate(&m_args)) {
        throw error("Unable to create a m_args object");
      }

      fly::verify(m_args, "Compute arguments pointer not assigned");

      // check compute arguments
      int num_comp_names = 0;

      namespace KAN = KIM::COMPUTE_ARGUMENT_NAME;

      KAN::GetNumberOfComputeArgumentNames(&num_comp_names);

      for (int i = 0; i < num_comp_names; ++i) {
        KIM::ComputeArgumentName name;
        KIM::SupportStatus supportStatus;
        KAN::GetComputeArgumentName(i, &name);
        KIM::DataType dataType;
        KAN::GetComputeArgumentDataType(name, &dataType);

        if (m_args->GetArgumentSupportStatus(name, &supportStatus)) {
          throw error("unable to get ComputeArgument SupportStatus");
        }

        fmt::print("Compute argument \"{}\" is of type \"{}\" and has support status \"{}\"\n",
                   name.ToString(),
                   dataType.ToString(),
                   supportStatus.ToString());

        // Everything must be optional
        if (supportStatus == KIM::SUPPORT_STATUS::required) {
          throw error("unsupported required ComputeArgument: \"{}\"", name.ToString());
        }

        // Must have energy and forces as optional arguments
        if ((name == KAN::partialEnergy) || (name == KAN::partialForces)) {
          if (supportStatus != KIM::SUPPORT_STATUS::optional) {
            throw error("Energy or Forces are not available");
          }
        }
      }

      // Check compute callbacks
      int num_call_names = 0;
      KIM::COMPUTE_CALLBACK_NAME::GetNumberOfComputeCallbackNames(&num_call_names);

      for (int i = 0; i < num_call_names; ++i) {
        KIM::ComputeCallbackName name;
        KIM::COMPUTE_CALLBACK_NAME::GetComputeCallbackName(i, &name);
        KIM::SupportStatus supportStatus;
        m_args->GetCallbackSupportStatus(name, &supportStatus);

        fmt::print("Callback \"{}\" has support status \"{}\"\n", name.ToString(), supportStatus.ToString());

        // Cannot handle any "required" callbacks
        if (supportStatus == KIM::SUPPORT_STATUS::required) {
          throw error("unsupported required ComputeCallback");
        }
      }

      // Check supported extensions, if any
      int present = 0;

      if (m_model->IsRoutinePresent(KIM::MODEL_ROUTINE_NAME::Extension, &present, NULL)) {
        throw error("Unable to get Extension present/required.");
      }

      if (present) {
        KIM::SupportedExtensions supportedExtensions;

        if (m_model->Extension(KIM_SUPPORTED_EXTENSIONS_ID, &supportedExtensions)) {
          throw error("Error returned from KIM::Model::Extension().");
        }
      }

      // Check neighbour list requirements

      int const* modelWillNotRequestNeighborsOfNoncontributingParticles;
      double const* cutoff_cluster_model;
      int number_of_neighbor_lists;

      m_model->GetNeighborListPointers(&number_of_neighbor_lists,
                                       &cutoff_cluster_model,
                                       &modelWillNotRequestNeighborsOfNoncontributingParticles);

      fmt::print("Model will request {} neighbour lists\n", number_of_neighbor_lists);

      for (int i = 0; i < number_of_neighbor_lists; ++i) {
        fmt::print("\t List {} has cut-off {} and will-not-request-neigh-of-non-contributing={}\n",
                   i,
                   cutoff_cluster_model[i],
                   bool(modelWillNotRequestNeighborsOfNoncontributingParticles[i]));

        if (!modelWillNotRequestNeighborsOfNoncontributingParticles[i]) {
          throw error("Cannot provide neigh-lists for this KIM potential");
        }
      }
      // ignoring hints from here on...
      if (number_of_neighbor_lists != 1) {
        throw error("Kim model needs {} neighbour lists!", number_of_neighbor_lists);
      }

      // We're compatible with the model. Let's do it.
    }

    ~KIM_API() noexcept {
      if (m_model) {
        if (m_args) {
          if (m_model->ComputeArgumentsDestroy(&m_args)) {
            throw error("KIM failed to clean-up model");
          }
          m_args = nullptr;
        }

        KIM::Model::Destroy(&m_model);

        m_model = nullptr;
      } else {
        verify(m_args == nullptr, "KIM_API in invalid state for clean-up");
      }
    }

    /**
     * @brief Get this potentials cut-off radius.
     *
     * This is the maximum distance two atom can interact. The neighbour::List passed to the other functions
     * should be configured with a cut-off equal or greater than this.
     */
    auto r_cut() const noexcept -> double {
      double influence_distance_cluster_model = 0;
      m_model->GetInfluenceDistance(&influence_distance_cluster_model);
      return influence_distance_cluster_model;
    }

    /**
     * @brief Compute potential energy gradient.
     *
     * Assumes the neighbour list are ready, force on frozen atoms will be zero.
     *
     * @param in Per-atom TypeID's and Frozen properties
     * @param out Potential gradient written to this.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     */
    auto gradient(system::SoA<PotentialGradient&> out,
                  system::SoA<TypeID const&, Frozen const&> in,
                  neigh::List const& nl,
                  int) -> void {
      //

      struct Code : system::Property<int> {};

      system::SoA<Code> codes(in.size());

      for (Eigen::Index i = 0; i < in.size(); i++) {
        codes(Code{}, i) = safe_cast<int>(in(id_, i));
      }

      int num_plus_ghosts = nl.m_num_plus_ghosts;

      {
        namespace KAN = KIM::COMPUTE_ARGUMENT_NAME;

        if (m_args->SetArgumentPointer(KAN::numberOfParticles, &num_plus_ghosts)
            || m_args->SetArgumentPointer(KAN::particleSpeciesCodes, codes[Code{}].data())
            || m_args->SetArgumentPointer(KAN::particleContributing,
                                          nl.m_atoms[neigh::List::Contrib{}].data())
            || m_args->SetArgumentPointer(KAN::coordinates, nl.m_atoms[r_].data())
            || m_args->SetArgumentPointer(KAN::partialEnergy, static_cast<double*>(nullptr))
            || m_args->SetArgumentPointer(KAN::partialForces, out[g_].data())) {
          throw("KIM_API_set_data");
        }
      }

      //   if (m_args->SetCallbackPointer(KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
      //                                  KIM::LANGUAGE_NAME::cpp,
      //                                  (KIM::Function*)&neigh::List::get_cluster_neigh,
      //                                  (void*)&nl)) {
      //     throw error("KIM_API_set_call_back");
      //   }

      //   kim_cluster_model->GetInfluenceDistance(&influence_distance_cluster_model);
      //   int const* modelWillNotRequestNeighborsOfNoncontributingParticles;
      //   kim_cluster_model->GetNeighborListPointers(&number_of_neighbor_lists,
      //                                              &cutoff_cluster_model,
      //                                              &modelWillNotRequestNeighborsOfNoncontributingParticles);
      //   std::cout << "Model has influence distance of : " << influence_distance_cluster_model << std::endl;
      //   std::cout << "Model has numberOfNeighborLists : " << number_of_neighbor_lists << std::endl;
      //   for (int i = 0; i < number_of_neighbor_lists; ++i) {
      //     std::cout << "\t"
      //               << "Neighbor list " << i << " has cutoff " << cutoff_cluster_model[i]
      //               << " with "
      //                  "modelWillNotRequestNeighborsOfNoncontributingParticles "
      //               << modelWillNotRequestNeighborsOfNoncontributingParticles[i] << std::endl;
      //   }
      //   // ignoring hints from here on...
      //   if (number_of_neighbor_lists != 1) MY_ERROR("too many neighbor lists");

      // /* setup particleSpecies */
      // int isSpeciesSupported;
      // error = kim_cluster_model->GetSpeciesSupportAndCode(
      //     KIM::SPECIES_NAME::Ar, &isSpeciesSupported, &(particleSpecies_cluster_model[0]));
      // if (error) MY_ERROR("get_species_code");
      // for (i = 1; i < NCLUSTERPARTS; ++i) particleSpecies_cluster_model[i] =
      // particleSpecies_cluster_model[0];
      // /* setup particleContributing */
      // for (i = 0; i < NCLUSTERPARTS; ++i)
      //   particleContributing_cluster_model[i] = 1; /* every particle contributes */

      // /* setup neighbor lists */
      // /* allocate memory for list */
      // nl_cluster_model.numberOfParticles = NCLUSTERPARTS;
      // nl_cluster_model.NNeighbors = new int[NCLUSTERPARTS];
      // if (NULL == nl_cluster_model.NNeighbors) MY_ERROR("new unsuccessful");

      // nl_cluster_model.neighborList = new int[NCLUSTERPARTS * NCLUSTERPARTS];
      // if (NULL == nl_cluster_model.neighborList) MY_ERROR("new unsuccessful");
    }

    /**
     * @brief Compute the potential energy.
     *
     * Assumes the neighbour list are ready, ignores contributions from the frozen atoms.
     *
     * @param in Per-atom TypeID's and Frozen properties.
     * @param nl Neighbour list (in ready state i.e. neigh::List::update() or neigh::List::rebuild() called).
     * @param threads Number of openMP threads to use.
     * @return double The potential energy of the system of atoms.
     */
    auto energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int threads = 1)
        -> double {}

  private:
    Options m_opt;

    std::vector<int> m_contrib;

    KIM::Model* m_model = nullptr;
    KIM::ComputeArguments* m_args = nullptr;
  };

}  // namespace fly::potential

using namespace fly;

template <typename... T>
system::Supercell<system::TypeMap<>, Position, Frozen, T...> bcc_iron_motif(double a = 2.855300) {
  //
  system::TypeMap<> FeH(2);

  FeH.set(0, tp_, "Fe");
  FeH.set(1, tp_, "H");

  Mat basis{
      {a, 0, 0},
      {0, a, 0},
      {0, 0, a},
  };

  system::Supercell motif
      = system::make_supercell<Position, Frozen, T...>({basis, Arr<bool>::Constant(true)}, FeH, 2);

  motif[fzn_] = false;
  motif[id_] = 0;

  motif(r_, 0) = Vec::Zero();
  motif(r_, 1) = Vec::Constant(0.5);

  return motif;
}

int main(int, const char**) {
  //

  system::Supercell perfect = motif_to_lattice(bcc_iron_motif(), {6, 6, 6});

  system::Supercell cell = remove_atoms(perfect, {1, 3});

  Vec r_H = {2.857 / 2 + 3.14, 2.857 / 2 + 3.14, 2.857 / 4 + 3.14};

  cell = add_atoms(cell, {system::Atom<TypeID, Position, Frozen>(1, r_H, false)});

  ///

  potential::KIM_API::Options opt{
      .model_name = "EAM_Dynamo_Wen_2021_FeH__MO_634187028437_000",
  };

  potential::KIM_API model(opt, cell.map());

  neigh::List nl(cell.box(), model.r_cut());

  nl.rebuild(cell, omp_get_max_threads());

  system::SoA<PotentialGradient> grad(cell.size());

  model.gradient(grad, cell, nl, 0);

  fmt::print("r_cut = {}\n", model.r_cut());

  fmt::print("Working\n");

  return 0;
}