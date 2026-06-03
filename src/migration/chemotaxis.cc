// =============================================================================
//
//   Copyright (C) 2020-2024 Vasileios Vavourakis (vasvav@gmail.com)
//   All Rights Reserved.
//
//   Licensed under the GNU General Public License v3.0 (the "License").
//   See the LICENSE file provided in this project details the License.
//   You cannot use this file except in compliance with the License.
//
// =============================================================================

#include "migration/chemotaxis.h"
#include "migration/candidate_selection.h"
#include "migration/field_sampling.h"
#include "biological_cell.h"
#include "global.h"
#include <functional>
#include "core/environment/uniform_grid_environment.h"
#include "core/execution_context/execution_context.h"
#include "core/resource_manager.h"
#include "core/simulation.h"
// =============================================================================
namespace bdm {
namespace migration {
// -----------------------------------------------------------------------------
bool RunChemotaxis(BiologicalCell* cell,
                   const std::string& CP_name,
                   const double step_length,
                   Double3& active_displacement,
                   bool& has_migrated)
{
  if (!cell) return false;
  auto* params = cell->params();
  if (!params) return false;
  auto* rm = Simulation::GetActive()->GetResourceManager();
  auto* rg = Simulation::GetActive()->GetRandom();
  const double minCOORD = params->get<double>("min_boundary"),
               maxCOORD = params->get<double>("max_boundary");
  const double tol = params->get<double>("domain_tolerance");
  const std::vector<std::string>& substances =
    params->get<std::vector<std::string>>("substances");
  const bool simulation_domain_is_2D =
    params->get<bool>("simulation_domain_is_2D");
  const Double3 current_position = cell->GetPosition();
  const double migration_tolerance =
    params->get<double>("migration_tolerance");
  const double self_diameter = cell->GetDiameter();
  auto* env = Simulation::GetActive()->GetEnvironment();
  const auto* uniform_env = dynamic_cast<UniformGridEnvironment*>(env);
  const std::array<int32_t, 6> env_dims =
    uniform_env ? uniform_env->GetDimensions()
                : std::array<int32_t, 6>{0, 0, 0, 0, 0, 0};
  const double env_box_length =
    uniform_env ? static_cast<double>(uniform_env->GetBoxLength()) : 0.0;
  const double neighbor_search_radius =
    uniform_env ? env_box_length : env->GetLargestAgentSize();
  auto* ctxt = Simulation::GetActive()->GetExecutionContext();
  const std::vector<Double3> local_directions =
    BuildShuffledLocalDirections(simulation_domain_is_2D, rg);
  if (!check_agent_position_in_domain(minCOORD, maxCOORD, current_position, tol))
    return false;
  if (step_length <= migration_tolerance) return false;
  bool stop_substances = false;
  auto is_candidate_feasible = [&](const Double3& pt) -> bool {
    if (!check_agent_position_in_domain(minCOORD, maxCOORD, pt, tol)) return false;
    if (uniform_env)
      {
        if (pt[0] < env_dims[0] + env_box_length ||
            pt[0] > env_dims[1] - env_box_length ||
            pt[1] < env_dims[2] + env_box_length ||
            pt[1] > env_dims[3] - env_box_length ||
            pt[2] < env_dims[4] + env_box_length ||
            pt[2] > env_dims[5] - env_box_length)
          return false;
      }
    bool feasible = true;
    auto has_free_space = L2F([&](Agent* neighbor, real_t) {
      if (!feasible) return;
      if (neighbor->GetUid() == cell->GetUid()) return;
      auto* other_cell = dynamic_cast<BiologicalCell*>(neighbor);
      if (!other_cell) return;
      const double min_distance =
        0.5 * (self_diameter + other_cell->GetDiameter()) - 1.0e-6;
      if (L2norm(pt - other_cell->GetPosition()) < min_distance)
        feasible = false;
    });
    ctxt->ForEachNeighbor(has_free_space, pt, pow2(neighbor_search_radius));
    return feasible;
  };
  auto crowding_score_adj = [&](const Double3& pt, double& adj) -> bool {
    adj = 0.0;
    if (!params->have_parameter<bool>(CP_name + "/can_migrate/use_crowding"))
      return true;
    if (!params->get<bool>(CP_name + "/can_migrate/use_crowding"))
      return true;
    const double cr =
      params->get<double>(CP_name + "/can_migrate/crowding_influence_ratio");
    if (cr <= 0.0) return true;
    const double occupancy = cell->ComputeLocalOccupancyRatio(pt, cr);
    const double max_occ =
      params->get<double>(CP_name + "/can_migrate/max_candidate_occupancy");
    if (occupancy >= max_occ) return false;
    const double penalty =
      params->get<double>(CP_name + "/can_migrate/crowding_penalty");
    adj = -penalty * occupancy;
    return true;
  };
  const std::function<bool(const Double3&, double&)> evaluate_candidate =
    [&](const Double3& pt, double& adj) -> bool {
      if (!is_candidate_feasible(pt)) return false;
      return crowding_score_adj(pt, adj);
    };
  for (std::vector<std::string>::const_iterator ci = substances.begin();
       ci != substances.end(); ci++)
    {
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string& BC_name = dg->GetContinuumName();
      if (!params->have_parameter<double>(CP_name + "/can_migrate/chemotaxis/" + BC_name))
        continue;
      const double chemotaxis =
        params->get<double>(CP_name + "/can_migrate/chemotaxis/" + BC_name);
      if (!chemotaxis) continue;
      const double concentration =
        GetInterpolatedValue(dg, current_position, params);
      const double threshold =
        params->have_parameter<double>(
          CP_name + "/can_migrate/chemotaxis/" + BC_name + "/threshold")
        ? params->get<double>(
            CP_name + "/can_migrate/chemotaxis/" + BC_name + "/threshold")
        : 0.0;
      if (!IsChemotaxisActive(concentration, threshold)) continue;
      const double probability =
        params->have_parameter<double>(
          CP_name + "/can_migrate/chemotaxis/" + BC_name + "/probability")
        ? params->get<double>(
            CP_name + "/can_migrate/chemotaxis/" + BC_name + "/probability")
        : 1.0;
      if (rg->Uniform(0.0, 1.0) > probability) continue;
      const std::string chemo_mode =
        params->have_parameter<std::string>(
          CP_name + "/can_migrate/chemotaxis/" + BC_name + "/mode")
        ? params->get<std::string>(
            CP_name + "/can_migrate/chemotaxis/" + BC_name + "/mode")
        : (params->have_parameter<std::string>(CP_name + "/can_migrate/chemotaxis_mode")
           ? params->get<std::string>(CP_name + "/can_migrate/chemotaxis_mode")
           : std::string("local_best"));
      const double chemotaxis_sign = (chemotaxis > 0.0 ? +1.0 : -1.0);
      const std::function<double(const Double3&)> sample_concentration =
        [&](const Double3& pt) {
          return GetInterpolatedValue(dg, pt, params);
        };
      Double3 move_dvec = {0.0, 0.0, 0.0};
      bool applied = false;
      if (chemo_mode == "gradient")
        {
          Double3 grad = GetInterpolatedGradient(dg, current_position, params);
          if (simulation_domain_is_2D) grad[2] = 0.0;
          const double grad_mag = L2norm(grad);
          if (grad_mag >= migration_tolerance)
            {
              const Double3 step_dir = grad * (chemotaxis_sign / grad_mag);
              Double3 point = current_position + step_dir * step_length;
              if (simulation_domain_is_2D) point[2] = current_position[2];
              if (is_candidate_feasible(point))
                {
                  move_dvec = point - current_position;
                  if (L2norm(move_dvec) > migration_tolerance)
                    applied = true;
                }
            }
        }
      else if (chemo_mode == "local_best")
        {
          const double tie_epsilon =
            params->have_parameter<double>(
              CP_name + "/can_migrate/chemotaxis/" + BC_name + "/tie_epsilon")
            ? params->get<double>(
                CP_name + "/can_migrate/chemotaxis/" + BC_name + "/tie_epsilon")
            : 1.0e-6 * std::max(1.0, fabs(concentration));
          applied = SelectLocalBestCandidate(
            local_directions, current_position, step_length,
            simulation_domain_is_2D, chemotaxis_sign, concentration, tie_epsilon,
            migration_tolerance, evaluate_candidate, sample_concentration, move_dvec);
        }
      else if (chemo_mode == "stochastic_local")
        {
          const double beta =
            params->have_parameter<double>(
              CP_name + "/can_migrate/chemotaxis/" + BC_name + "/beta")
            ? params->get<double>(
                CP_name + "/can_migrate/chemotaxis/" + BC_name + "/beta")
            : 10.0;
          applied = SelectStochasticLocalCandidate(
            local_directions, current_position, step_length,
            simulation_domain_is_2D, chemotaxis_sign, concentration, beta,
            migration_tolerance, rg, evaluate_candidate, sample_concentration, move_dvec);
        }
      else
        ABORT_("unrecognized chemotaxis mode \"" + chemo_mode
               + "\" for substance " + BC_name
               + ": use gradient/local_best/stochastic_local");
      if (applied)
        {
          active_displacement += move_dvec;
          has_migrated = true;
          if (!params->get<bool>(CP_name + "/can_migrate/accumulate_path"))
            {
              stop_substances = true;
              break;
            }
        }
    }
  return stop_substances;
}
// =============================================================================
}  // namespace migration
}  // namespace bdm
