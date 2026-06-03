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

#ifndef _INTRACELLULAR_DDR_PATHWAY_H_
#define _INTRACELLULAR_DDR_PATHWAY_H_
// =============================================================================
#include "global.h"
#include <string>
// =============================================================================
namespace bdm {
class BiologicalCell;
// -----------------------------------------------------------------------------
/// Intracellular DDR signalling nodes (normalized activity in [0, 1]).
enum class DdrMolecule {
  kAtm,
  kAtr,
  kChk1,
  kChk2,
  kP53,
  kP21,
  kCdc25,
  kCdk
};
// -----------------------------------------------------------------------------
/// Update ATM/ATR–CHK–p53–p21–Cdc25–CDK dynamics from `dna_damage_`.
void UpdateDdrPathway(BiologicalCell* cell);
// -----------------------------------------------------------------------------
/// True when molecular G1/S checkpoint is active (p21/CDK gate).
bool IsMolecularG1SCheckpointBlocked(const BiologicalCell* cell);
// -----------------------------------------------------------------------------
/// True when molecular G2/M checkpoint is active (CHK1/CHK2/CDK gate).
bool IsMolecularG2MCheckpointBlocked(const BiologicalCell* cell);
// -----------------------------------------------------------------------------
bool IsDdrPathwayEnabled(const ::Parameters& params,
                         const std::string& phenotype_name);
// =============================================================================
}  // namespace bdm
// =============================================================================
#endif  // _INTRACELLULAR_DDR_PATHWAY_H_
