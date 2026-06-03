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

#include "intracellular/ddr_pathway.h"
#include "model_parameters.h"
#include "biological_cell.h"
#include <algorithm>
#include <cmath>
// =============================================================================
namespace bdm {
namespace {
// -----------------------------------------------------------------------------
inline
double Clamp01(const double value)
{
  return std::clamp(value, 0.0, 1.0);
}
// -----------------------------------------------------------------------------
inline
void ActivateNode(const double dt,
                  const double signal,
                  const double k_on,
                  const double k_off,
                  double& level)
{
  const double drive = k_on * std::max(0.0, signal) * (1.0 - level);
  const double decay = k_off * level;
  level = Clamp01(level + dt * (drive - decay));
}
// -----------------------------------------------------------------------------
inline
double ParamOr(const ::Parameters& params,
               const std::string& name,
               const double fallback)
{
  return params.have_parameter<double>(name)
    ? params.get<double>(name) : fallback;
}
// -----------------------------------------------------------------------------
inline
bool ParamEnabled(const ::Parameters& params, const std::string& name)
{
  return params.have_parameter<bool>(name) && params.get<bool>(name);
}
// =============================================================================
}  // namespace
// -----------------------------------------------------------------------------
bool IsDdrPathwayEnabled(const ::Parameters& params,
                         const std::string& phenotype_name)
{
  const std::string flag = phenotype_name + "/intracellular/ddr/enabled";
  if (params.have_parameter<bool>(flag))
    return params.get<bool>(flag);
  return true;
}
// -----------------------------------------------------------------------------
void UpdateDdrPathway(BiologicalCell* cell)
{
  if (!cell || !cell->GetPhenotype() || !cell->params()) return;
  const std::string& CP_name =
    cell->params()->get<std::string>("phenotype_ID/"
                                     + std::to_string(cell->GetPhenotype()));
  if (!IsDdrPathwayEnabled(*cell->params(), CP_name)) return;
  const double dt = cell->params()->get<double>("time_step");
  const ::Parameters& params = *cell->params();
  const std::string base = CP_name + "/intracellular/ddr/";
  const double damage_signal = std::max(0.0, cell->GetDNADamage());
  // --- Sensor kinases (dsDNA breaks -> ATM; stalled forks -> ATR) ---
  ActivateNode(dt, damage_signal,
               ParamOr(params, base + "ATM/k_activation", 2.0),
               ParamOr(params, base + "ATM/k_deactivation", 0.5),
               cell->atm_active_);
  ActivateNode(dt, damage_signal,
               ParamOr(params, base + "ATR/k_activation", 1.5),
               ParamOr(params, base + "ATR/k_deactivation", 0.4),
               cell->atr_active_);
  // --- Effector kinases downstream of ATM/ATR ---
  const double chk1_signal =
    ParamOr(params, base + "CHK1/ATR_weight", 1.0) * cell->atr_active_
    + ParamOr(params, base + "CHK1/ATM_weight", 0.3) * cell->atm_active_;
  ActivateNode(dt, chk1_signal,
               ParamOr(params, base + "CHK1/k_activation", 1.2),
               ParamOr(params, base + "CHK1/k_deactivation", 0.3),
               cell->chk1_active_);
  ActivateNode(dt, cell->atm_active_,
               ParamOr(params, base + "CHK2/k_activation", 1.0),
               ParamOr(params, base + "CHK2/k_deactivation", 0.25),
               cell->chk2_active_);
  // --- p53 stabilisation (ATM/CHK2) and p21 CDK inhibitor ---
  const double p53_signal =
    cell->atm_active_
    + ParamOr(params, base + "p53/CHK2_weight", 0.8) * cell->chk2_active_;
  ActivateNode(dt, p53_signal,
               ParamOr(params, base + "p53/k_activation", 0.8),
               ParamOr(params, base + "p53/k_deactivation", 0.15),
               cell->p53_active_);
  ActivateNode(dt, cell->p53_active_,
               ParamOr(params, base + "p21/k_activation", 1.0),
               ParamOr(params, base + "p21/k_deactivation", 0.2),
               cell->p21_level_);
  // --- Cdc25 phosphatase (promotes CDK activation; inhibited by CHK1/CHK2) ---
  const double chk_inhibition =
    ParamOr(params, base + "Cdc25/CHK1_weight", 0.7) * cell->chk1_active_
    + ParamOr(params, base + "Cdc25/CHK2_weight", 0.8) * cell->chk2_active_;
  const double cdc25_baseline =
    ParamOr(params, base + "Cdc25/baseline_activity", 1.0);
  cell->cdc25_active_ =
    Clamp01(cdc25_baseline * (1.0 - std::min(1.0, chk_inhibition)));
  // --- CDK activity (Cdc25-dependent, p21-inhibited) ---
  const double cdk_baseline =
    ParamOr(params, base + "CDK/baseline_activity", 1.0);
  cell->cdk_activity_ = Clamp01(
    cdk_baseline * cell->cdc25_active_ * (1.0 - cell->p21_level_));
}
// -----------------------------------------------------------------------------
bool IsMolecularG1SCheckpointBlocked(const BiologicalCell* cell)
{
  if (!cell || !cell->GetPhenotype() || !cell->params()) return false;
  const std::string& CP_name =
    cell->params()->get<std::string>("phenotype_ID/"
                                     + std::to_string(cell->GetPhenotype()));
  if (!IsDdrPathwayEnabled(*cell->params(), CP_name)) return false;
  const ::Parameters& params = *cell->params();
  const std::string base = CP_name + "/intracellular/ddr/checkpoint/G1S/";
  const double p21_thr =
    ParamOr(params, base + "p21_threshold", 0.35);
  if (cell->GetP21Level() >= p21_thr) return true;
  const double cdk_thr =
    ParamOr(params, base + "CDK_min_activity", 0.45);
  if (cell->GetCdkActivity() < cdk_thr) return true;
  if (ParamEnabled(params, base + "use_p53")
      && cell->GetP53Active() >= ParamOr(params, base + "p53_threshold", 0.5))
    return true;
  return false;
}
// -----------------------------------------------------------------------------
bool IsMolecularG2MCheckpointBlocked(const BiologicalCell* cell)
{
  if (!cell || !cell->GetPhenotype() || !cell->params()) return false;
  const std::string& CP_name =
    cell->params()->get<std::string>("phenotype_ID/"
                                     + std::to_string(cell->GetPhenotype()));
  if (!IsDdrPathwayEnabled(*cell->params(), CP_name)) return false;
  const ::Parameters& params = *cell->params();
  const std::string base = CP_name + "/intracellular/ddr/checkpoint/G2M/";
  const double chk1_thr =
    ParamOr(params, base + "CHK1_threshold", 0.25);
  if (cell->GetChk1Active() >= chk1_thr) return true;
  const double chk2_thr =
    ParamOr(params, base + "CHK2_threshold", 0.20);
  if (cell->GetChk2Active() >= chk2_thr) return true;
  const double atm_thr =
    ParamOr(params, base + "ATM_threshold", 0.30);
  if (cell->GetAtmActive() >= atm_thr) return true;
  const double cdk_thr =
    ParamOr(params, base + "CDK_min_activity", 0.40);
  if (cell->GetCdkActivity() < cdk_thr) return true;
  return false;
}
// =============================================================================
}  // namespace bdm
