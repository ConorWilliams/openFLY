// Copyright © 2020-2022 Conor Williams <conorwilliams@outlook.com>

// SPDX-License-Identifier: GPL-3.0-or-later

// This file is part of openFLY.

// OpenFLY is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// OpenFLY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License along with openFLY. If not, see
// <https://www.gnu.org/licenses/>.

#include <algorithm>
#include <libfly/potential/KIM/kim.hpp>
//
#include "KIM_SimulatorHeaders.hpp"
#include "KIM_SupportedExtensions.hpp"
#include "libfly/system/SoA.hpp"
#include "libfly/system/property.hpp"
#include "libfly/utility/core.hpp"

namespace fly::potential {

  KIM_API::KIM_API(Options const& opt, system::TypeMap<> const& map) : m_opt(opt), m_map(map.num_types()) {
    //
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

    ASSERT(m_model, "Model pointer not assigned", 0);

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

      dprint(m_opt.debug,
             "KIM: Model routine name \"{}\": present = {} and required = {}\n",
             name.ToString(),
             bool(present),
             bool(required));

      if (present == true && required == true) {
        //
        std::array names = {
            KIM::MODEL_ROUTINE_NAME::Create,
            KIM::MODEL_ROUTINE_NAME::ComputeArgumentsCreate,
            KIM::MODEL_ROUTINE_NAME::Compute,
            KIM::MODEL_ROUTINE_NAME::Refresh,
            KIM::MODEL_ROUTINE_NAME::ComputeArgumentsDestroy,
            KIM::MODEL_ROUTINE_NAME::Destroy,
        };

        if (std::find(names.begin(), names.end(), name) == names.end()) {
          throw fly::error("Unknown Routine \"{}\" is required by model.", name.ToString());
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

    dprint(m_opt.debug, "KIM: Length unit is \"{}\"\n", lengthUnit.ToString());
    dprint(m_opt.debug, "KIM: Energy unit is \"{}\"\n", energyUnit.ToString());
    dprint(m_opt.debug, "KIM: Charge unit is \"{}\"\n", chargeUnit.ToString());
    dprint(m_opt.debug, "KIM: Temperature unit is \"{}\"\n", temperatureUnit.ToString());
    dprint(m_opt.debug, "KIM: Time unit is \"{}\"\n", timeUnit.ToString());

    // Check species

    for (Eigen::Index i = 0; i < map.num_types(); i++) {
      int supported;
      int model_code;

      auto j = safe_cast<std::uint32_t>(i);

      KIM::SpeciesName name = std::string{map.get(j, tp_)};

      if (m_model->GetSpeciesSupportAndCode(name, &supported, &model_code) || !supported) {
        throw error("Species \"{}\" unsupported by KIM potential \"{}\"", name.ToString(), m_opt.model_name);
      }

      m_map.set(j, tp_, map.get(j, tp_));
      m_map.set(j, KIM_code{}, model_code);
    }

    // Set up compute request

    if (m_model->ComputeArgumentsCreate(&m_args)) {
      throw error("Unable to create a m_args object");
    }

    ASSERT(m_args, "Compute arguments pointer not assigned", 0);

    // Check compute arguments.
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

      dprint(m_opt.debug,
             "KIM: Compute argument \"{}\" is of type \"{}\" and has support status \"{}\"\n",
             name.ToString(),
             dataType.ToString(),
             supportStatus.ToString());

      // Everything must be optional
      if (supportStatus == KIM::SUPPORT_STATUS::required) {
        throw error("Unsupported required ComputeArgument: \"{}\"", name.ToString());
      }

      // Must have energy and forces as optional arguments
      if ((name == KAN::partialEnergy) || (name == KAN::partialForces)) {
        if (supportStatus != KIM::SUPPORT_STATUS::optional) {
          throw error("KIM - Energy or Forces are not available");
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

      dprint(m_opt.debug,
             "KIM: Callback \"{}\" has support status \"{}\"\n",
             name.ToString(),
             supportStatus.ToString());

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

    dprint(m_opt.debug, "KIM: Model will request {} neighbour lists\n", number_of_neighbor_lists);

    for (int i = 0; i < number_of_neighbor_lists; ++i) {
      dprint(m_opt.debug,
             "\tKIM: List {} has cut-off {} and will-not-request-neigh-of-non-contributing={}\n",
             i,
             cutoff_cluster_model[i],
             bool(modelWillNotRequestNeighborsOfNoncontributingParticles[i]));

      if (!modelWillNotRequestNeighborsOfNoncontributingParticles[i]) {
        throw error("Cannot provide ghost neigh-lists for this KIM potential");
      }
    }

    if (number_of_neighbor_lists != 1) {
      throw error("Kim model needs {} neighbour lists!", number_of_neighbor_lists);
    }

    // We're compatible with the model. Let's do it.
  }

  KIM_API::~KIM_API() noexcept {
    if (m_model) {
      if (m_args) {
        if (m_model->ComputeArgumentsDestroy(&m_args)) {
          throw error("KIM failed to clean-up compute arguments");
        }
      }
      KIM::Model::Destroy(&m_model);
    }
  }

  auto KIM_API::r_cut() const noexcept -> double {
    double influence_distance_cluster_model = 0;
    m_model->GetInfluenceDistance(&influence_distance_cluster_model);
    return influence_distance_cluster_model;
  }

  double KIM_API::energy(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl, int) {
    //
    int num_plus_ghosts = safe_cast<int>(nl.m_num_plus_ghosts);

    if (m_kim_io.size() != num_plus_ghosts) {
      m_kim_io.destructive_resize(num_plus_ghosts);
    }

    for (Eigen::Index i = 0; i < num_plus_ghosts; i++) {
      m_kim_io(KIM_code{}, i) = m_map.get(in(id_, nl.image_to_real(i)), KIM_code{});
      m_kim_io(KIM_contrib{}, i) = i < nl.size();
    }

    namespace KAN = KIM::COMPUTE_ARGUMENT_NAME;

    if (m_args->SetArgumentPointer(KAN::numberOfParticles, &num_plus_ghosts)
        || m_args->SetArgumentPointer(KAN::particleSpeciesCodes, m_kim_io[KIM_code{}].data())
        || m_args->SetArgumentPointer(KAN::particleContributing, m_kim_io[KIM_contrib{}].data())
        || m_args->SetArgumentPointer(KAN::coordinates, nl.m_atoms[r_].data())
        || m_args->SetArgumentPointer(KAN::partialForces, static_cast<double*>(nullptr))
        || m_args->SetArgumentPointer(KAN::partialEnergy, m_kim_io[KIM_energy{}].data())) {
      throw("KIM_API_set_argument_pointer");
    }

    if (m_args->SetCallbackPointer(KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
                                   KIM::LANGUAGE_NAME::cpp,
                                   reinterpret_cast<void (*)()>(&neigh::List::get_neigh),
                                   const_cast<neigh::List*>(&nl)  // C api has no const-correctness
                                   )) {
      throw error("KIM_API_set_callback");
    }

    int res = 0;

    m_args->AreAllRequiredArgumentsAndCallbacksPresent(&res);

    verify(res, "Failed to set all required arguments and callbacks!");

    if (m_model->Compute(m_args)) {
      throw error("KIM_API_compute");
    }

    // Sum up the partial energies

    double energy = 0;

    for (int i = 0; i < num_plus_ghosts; i++) {
      energy += m_kim_io(KIM_energy{}, i);
    }

    return energy;
  }

  void KIM_API::force(system::SoA<TypeID const&, Frozen const&> in, neigh::List const& nl) {
    int num_plus_ghosts = safe_cast<int>(nl.m_num_plus_ghosts);

    if (m_kim_io.size() != num_plus_ghosts) {
      m_kim_io.destructive_resize(num_plus_ghosts);
    }

    for (Eigen::Index i = 0; i < num_plus_ghosts; i++) {
      m_kim_io(KIM_code{}, i) = m_map.get(in(id_, nl.image_to_real(i)), KIM_code{});
      m_kim_io(KIM_contrib{}, i) = i < nl.size();
    }

    namespace KAN = KIM::COMPUTE_ARGUMENT_NAME;

    if (m_args->SetArgumentPointer(KAN::numberOfParticles, &num_plus_ghosts)
        || m_args->SetArgumentPointer(KAN::particleSpeciesCodes, m_kim_io[KIM_code{}].data())
        || m_args->SetArgumentPointer(KAN::particleContributing, m_kim_io[KIM_contrib{}].data())
        || m_args->SetArgumentPointer(KAN::coordinates, nl.m_atoms[r_].data())
        || m_args->SetArgumentPointer(KAN::partialEnergy, static_cast<double*>(nullptr))
        || m_args->SetArgumentPointer(KAN::partialForces, m_kim_io[KIM_force{}].data())) {
      throw("KIM_API_set_argument_pointer");
    }

    if (m_args->SetCallbackPointer(KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
                                   KIM::LANGUAGE_NAME::cpp,
                                   reinterpret_cast<void (*)()>(&neigh::List::get_neigh),
                                   const_cast<neigh::List*>(&nl)  // C api has no const-correctness
                                   )) {
      throw error("KIM_API_set_callback");
    }

    int res = 0;

    m_args->AreAllRequiredArgumentsAndCallbacksPresent(&res);

    verify(res, "Failed to set all required arguments and callbacks!");

    if (m_model->Compute(m_args)) {
      throw error("KIM_API_compute");
    }
  }

  auto KIM_API::gradient(system::SoA<PotentialGradient&> out,
                         system::SoA<TypeID const&, Frozen const&> in,
                         neigh::List const& nl,
                         int) -> void {
    //

    verify(out.size() == in.size(), "size mismatch! {} != {}", out.size(), in.size());

    force(in, nl);

    out[g_] = 0;

    for (Eigen::Index i = 0; i < nl.m_num_plus_ghosts; ++i) {
      if (auto j = nl.image_to_real(i); !in(fzn_, j)) {
        out(g_, j) -= m_kim_io(KIM_force{}, i);
      }
    }
  }

  auto KIM_API::hessian(system::Hessian& out,
                        system::SoA<TypeID const&, Frozen const&> in,
                        neigh::List const& nl_initial,
                        int threads) -> void {
    // Comput the hessian using the symmetric finite difference of the gradient.

    out.zero_for(in.size());

    system::SoA<Position, Delta> r(in.size());

    Eigen::Index const N = in.size() * Position::size();

    r[r_] = nl_initial.m_atoms[r_].head(N);
    r[del_] = 0;

    system::SoA<PotentialGradient> g_plus(in.size());
    system::SoA<PotentialGradient> g_minus(in.size());

    neigh::List nl(nl_initial.m_box, r_cut() + 2 * m_opt.epsilon);

    nl.rebuild(r, threads);

    for (Eigen::Index i = 0; i < N; ++i) {
      r[del_][i] = -m_opt.epsilon;
      nl.update(r);
      gradient(g_plus, in, nl, threads);

      r[del_][i] = 2 * m_opt.epsilon;
      nl.update(r);
      gradient(g_minus, in, nl, threads);

      r[del_][i] = -m_opt.epsilon;
      nl.update(r);
      r[del_][i] = 0;

      out.get().col(i) = (g_plus[g_] - g_minus[g_]) / (2 * m_opt.epsilon);
    }

    // Symmetrise hessian
    // system::Hessian::Matrix H = out.get().transpose();
    // out.get() = 0.5 * (H + out.get());
  }

}  // namespace fly::potential