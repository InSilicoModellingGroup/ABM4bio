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

#include "migration/active_migration.h"
#include "migration/chemotaxis.h"
#include "biological_cell.h"
#include "global.h"
#include "core/simulation.h"
// =============================================================================
namespace bdm {
namespace migration {
// -----------------------------------------------------------------------------
bool RunActiveMigration(BiologicalCell* cell,
                        const std::string& phenotype_name,
                        Double3& active_displacement,
                        bool& has_migrated)
{
  if (!cell) return false;
  auto* params = cell->params();
  if (!params) return false;
  auto* rg = Simulation::GetActive()->GetRandom();
  const double time_step = params->get<double>("time_step");
  const double migration_tolerance =
    params->get<double>("migration_tolerance");
  const ActiveMigrationKinematics kinematics = ResolveActiveMigrationKinematics(
    *params, phenotype_name, time_step, cell->GetDiameter());
  if (kinematics.step_length <= migration_tolerance) return false;
  if (!ShouldAttemptMigration(rg, kinematics.event_probability)) return false;
  if (!ShouldReorient(rg, kinematics.reorientation_probability)) return false;
  active_displacement = {0.0, 0.0, 0.0};
  const bool simulation_domain_is_2d =
    params->get<bool>("simulation_domain_is_2D");
  const Double3 random_direction =
    SampleRandomMigrationDirection(simulation_domain_is_2d, rg);
  const Double3 random_dvec = random_direction * kinematics.step_length;
  if (L2norm(random_dvec) > migration_tolerance)
    {
      active_displacement += random_dvec;
      has_migrated = true;
    }
  return RunChemotaxis(
    cell, phenotype_name, kinematics.step_length, active_displacement, has_migrated);
}
// =============================================================================
}  // namespace migration
}  // namespace bdm
