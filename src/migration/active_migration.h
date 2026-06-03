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

#ifndef _MIGRATION_ACTIVE_MIGRATION_H_
#define _MIGRATION_ACTIVE_MIGRATION_H_
// =============================================================================
#include "core/container/math_array.h"
#include "neuroscience/param.h"
#include "global.h"
#include <cmath>
#include <string>
// =============================================================================
namespace bdm {
class BiologicalCell;
class Random;
namespace migration {
// -----------------------------------------------------------------------------
/// Per-timestep active migration kinematics derived from phenotype parameters.
struct ActiveMigrationKinematics {
  double step_length = 0.0;           ///< speed * dt, optionally diameter-capped
  double event_probability = 0.0;     ///< P(attempt migration this step)
  double reorientation_probability = 1.0; ///< P(redraw direction / run chemotaxis)
};
// -----------------------------------------------------------------------------
/// step = speed * time_step, or legacy half_range when speed is absent.
/// Optional cap: min(step, max_step_fraction_diameter * diameter).
inline
double ComputeMigrationStepLength(const ::Parameters& params,
                                  const std::string& phenotype_name,
                                  const double time_step,
                                  const double cell_diameter)
{
  const std::string speed_name = phenotype_name + "/can_migrate/speed";
  const std::string half_range_name = phenotype_name + "/can_migrate/half_range";
  double step = 0.0;
  if (params.have_parameter<double>(speed_name))
    step = params.get<double>(speed_name) * time_step;
  else if (params.have_parameter<double>(half_range_name))
    step = params.get<double>(half_range_name);
  const std::string cap_name =
    phenotype_name + "/can_migrate/max_step_fraction_diameter";
  if (params.have_parameter<double>(cap_name) && cell_diameter > 0.0)
    {
      const double cap = params.get<double>(cap_name) * cell_diameter;
      if (cap > 0.0) step = std::min(step, cap);
    }
  return step;
}
// -----------------------------------------------------------------------------
/// Per-step migration attempt probability.
/// Uses can_migrate/probability_rate [1/time] when set: 1 - exp(-rate * dt).
/// Otherwise falls back to legacy can_migrate/probability per step.
inline
double ComputeMigrationEventProbability(const ::Parameters& params,
                                        const std::string& phenotype_name,
                                        const double time_step)
{
  const std::string rate_name = phenotype_name + "/can_migrate/probability_rate";
  if (params.have_parameter<double>(rate_name))
    {
      const double rate = params.get<double>(rate_name);
      if (rate <= 0.0) return 0.0;
      return 1.0 - std::exp(-rate * time_step);
    }
  const std::string prob_name = phenotype_name + "/can_migrate/probability";
  if (params.have_parameter<double>(prob_name))
    return std::clamp(params.get<double>(prob_name), 0.0, 1.0);
  return 0.0;
}
// -----------------------------------------------------------------------------
/// Stochastic reorientation probability for this step.
/// persistence_time [time]: 1 - exp(-dt / persistence_time).
/// Legacy strength_of_time (integer >= 1): equivalent to persisting (N-1)/N of
/// steps — migrate only when (index_time-1) % strength_of_time == 0.
inline
double ComputeReorientationProbability(const ::Parameters& params,
                                       const std::string& phenotype_name,
                                       const double time_step)
{
  const std::string persistence_name =
    phenotype_name + "/can_migrate/persistence_time";
  if (params.have_parameter<double>(persistence_name))
    {
      const double persistence_time = params.get<double>(persistence_name);
      if (persistence_time <= 0.0) return 1.0;
      return 1.0 - std::exp(-time_step / persistence_time);
    }
  const std::string strength_name = phenotype_name + "/can_migrate/strength_of_time";
  if (params.have_parameter<int>(strength_name))
    {
      const int strength_time = std::max(1, params.get<int>(strength_name));
      if (strength_time <= 1) return 1.0;
      const int index_time = params.get<int>("index time");
      return ((index_time - 1) % strength_time == 0) ? 1.0 : 0.0;
    }
  return 1.0;
}
// -----------------------------------------------------------------------------
inline
ActiveMigrationKinematics ResolveActiveMigrationKinematics(
    const ::Parameters& params,
    const std::string& phenotype_name,
    const double time_step,
    const double cell_diameter)
{
  ActiveMigrationKinematics k;
  k.step_length = ComputeMigrationStepLength(
    params, phenotype_name, time_step, cell_diameter);
  k.event_probability = ComputeMigrationEventProbability(
    params, phenotype_name, time_step);
  k.reorientation_probability = ComputeReorientationProbability(
    params, phenotype_name, time_step);
  return k;
}
// -----------------------------------------------------------------------------
inline
bool ShouldAttemptMigration(Random* rg, const double event_probability)
{
  return rg->Uniform(0.0, 1.0) <= event_probability;
}
// -----------------------------------------------------------------------------
inline
bool ShouldReorient(Random* rg, const double reorientation_probability)
{
  if (reorientation_probability >= 1.0) return true;
  if (reorientation_probability <= 0.0) return false;
  return rg->Uniform(0.0, 1.0) <= reorientation_probability;
}
// -----------------------------------------------------------------------------
/// Random unit direction for undirected migration (8/26-neighbour symmetry).
inline
Double3 SampleRandomMigrationDirection(const bool simulation_domain_is_2d,
                                       Random* rg)
{
  Double3 direction = {
    rg->Uniform(-1.0, 1.0),
    rg->Uniform(-1.0, 1.0),
    simulation_domain_is_2d ? 0.0 : rg->Uniform(-1.0, 1.0)
  };
  if (!normalize(direction, direction))
    direction = {1.0, 0.0, 0.0};
  return direction;
}
// -----------------------------------------------------------------------------
/// Timestep-consistent random walk + chemotaxis for one cell.
bool RunActiveMigration(BiologicalCell* cell,
                        const std::string& phenotype_name,
                        Double3& active_displacement,
                        bool& has_migrated);
// =============================================================================
}  // namespace migration
}  // namespace bdm
// =============================================================================
#endif  // _MIGRATION_ACTIVE_MIGRATION_H_
