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

#ifndef _MIGRATION_CHEMOTAXIS_H_
#define _MIGRATION_CHEMOTAXIS_H_
// =============================================================================
#include "core/container/math_array.h"
#include <string>
// =============================================================================
namespace bdm {
class BiologicalCell;
namespace migration {
// -----------------------------------------------------------------------------
/// Active chemotactic migration for one cell over all configured substances.
/// Updates `active_displacement` and `has_migrated` when a move is chosen.
/// Returns true when no further substances should be processed this step
/// (`can_migrate/accumulate_path` is false and a move was applied).
/// `step_length` is speed * dt (optionally capped); chemotaxis weights are
/// directional bias only (sign/attract-repel), not step size.
bool RunChemotaxis(BiologicalCell* cell,
                   const std::string& phenotype_name,
                   const double step_length,
                   Double3& active_displacement,
                   bool& has_migrated);
// =============================================================================
}  // namespace migration
}  // namespace bdm
// =============================================================================
#endif  // _MIGRATION_CHEMOTAXIS_H_
