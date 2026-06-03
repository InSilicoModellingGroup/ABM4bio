// =============================================================================
// Inline implementation for BooleanRegulatoryModel.
// =============================================================================
#ifndef _BOOLEAN_GRN_INLINE_H_
#define _BOOLEAN_GRN_INLINE_H_
// =============================================================================
#include "./boolean_grn.h"
// =============================================================================
namespace bdm {
namespace regulatory {
// =============================================================================
namespace {
inline double GrnGetNumeric(const RegulatoryParameters& p,
                            const std::string& key,
                            double fallback) {
  auto it = p.numeric.find(key);
  return (it == p.numeric.end()) ? fallback : it->second;
}

inline bool GrnGetBool(const RegulatoryParameters& p,
                       const std::string& key,
                       bool fallback) {
  auto it = p.boolean.find(key);
  return (it == p.boolean.end()) ? fallback : it->second;
}

inline double Sigmoid(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

inline double AsBinary(double x) {
  return x > 0.5 ? 1.0 : 0.0;
}

inline void EvaluateCompactDdrEcmNodes(const RegulatoryParameters& params,
                                       const RegulatoryInput& in,
                                       std::map<std::string, double>* nodes) {
  if (!nodes) return;
  const double ros_thr = GrnGetNumeric(params, "ros_high_threshold", 0.35);
  const double dmg_thr = GrnGetNumeric(params, "dna_damage_high_threshold", 0.35);
  const double hyp_thr = GrnGetNumeric(params, "hypoxia_threshold", 0.15);
  const double nut_thr = GrnGetNumeric(params, "nutrient_low_threshold", 0.25);
  const double ecm_stiff_thr = GrnGetNumeric(params, "ecm_stiffness_high_threshold", 1.2);
  const double crowd_thr = GrnGetNumeric(params, "crowding_high_threshold", 0.6);
  const double tgfb_thr = GrnGetNumeric(params, "tgfb_high_threshold", 0.4);
  const double adh_low_thr = GrnGetNumeric(params, "adhesion_low_threshold", 0.2);

  const bool ros_high = (in.intracellular_ros > ros_thr) || (in.extracellular_rons > ros_thr);
  const bool damage_high = (in.dna_damage > dmg_thr);
  const bool hypoxia_high = (in.oxygen < hyp_thr);
  const bool nutrient_low = (in.nutrient < nut_thr);
  const bool crowding_high = (in.crowding > crowd_thr);
  const bool ecm_stiff_high = (in.ecm_stiffness > ecm_stiff_thr);
  const bool adhesion_low = (in.adhesion_signal < adh_low_thr);
  const bool tgfb_high = (in.tgfb > tgfb_thr);

  const bool nrf2 = ros_high && !damage_high;
  const bool atm_atr = damage_high;
  const bool chk1_chk2 = atm_atr;
  const bool p53 = damage_high;
  const bool p21 = p53;
  const bool repair = (!damage_high) || nrf2;
  const bool caspase3 = p53 && damage_high && !repair;
  const bool hif1a = hypoxia_high;
  const bool integrin_fak = !adhesion_low;
  const bool yap_taz = ecm_stiff_high && integrin_fak && !crowding_high;
  const bool emt = tgfb_high || hypoxia_high;
  const bool mmp = emt;
  const bool quiescence = crowding_high || hypoxia_high || nutrient_low;
  const bool proliferation = !quiescence && !p21;
  const bool migration = emt || yap_taz;

  (*nodes)["ROS_high"] = ros_high ? 1.0 : 0.0;
  (*nodes)["Damage_high"] = damage_high ? 1.0 : 0.0;
  (*nodes)["NRF2"] = nrf2 ? 1.0 : 0.0;
  (*nodes)["ATM_ATR"] = atm_atr ? 1.0 : 0.0;
  (*nodes)["CHK1_CHK2"] = chk1_chk2 ? 1.0 : 0.0;
  (*nodes)["p53"] = p53 ? 1.0 : 0.0;
  (*nodes)["p21"] = p21 ? 1.0 : 0.0;
  (*nodes)["Repair"] = repair ? 1.0 : 0.0;
  (*nodes)["Caspase3"] = caspase3 ? 1.0 : 0.0;
  (*nodes)["HIF1A"] = hif1a ? 1.0 : 0.0;
  (*nodes)["Integrin_FAK"] = integrin_fak ? 1.0 : 0.0;
  (*nodes)["YAP_TAZ"] = yap_taz ? 1.0 : 0.0;
  (*nodes)["EMT"] = emt ? 1.0 : 0.0;
  (*nodes)["MMP"] = mmp ? 1.0 : 0.0;
  (*nodes)["Quiescence"] = quiescence ? 1.0 : 0.0;
  (*nodes)["Proliferation"] = proliferation ? 1.0 : 0.0;
  (*nodes)["Migration"] = migration ? 1.0 : 0.0;
}

inline void NodesToOutput(const RegulatoryParameters& params,
                          const std::map<std::string, double>& nodes,
                          RegulatoryOutput* out) {
  if (!out) return;
  const double apop_base = GrnGetNumeric(params, "apoptosis_hazard_base", 0.0);
  const double nec_base = GrnGetNumeric(params, "necrosis_hazard_base", 0.0);
  const double qui_base = GrnGetNumeric(params, "quiescence_hazard_base", 0.0);

  const double p21 = AsBinary(nodes.at("p21"));
  const double chk = AsBinary(nodes.at("CHK1_CHK2"));
  const double casp = AsBinary(nodes.at("Caspase3"));
  const double hif = AsBinary(nodes.at("HIF1A"));
  const double nrf2 = AsBinary(nodes.at("NRF2"));
  const double mmp = AsBinary(nodes.at("MMP"));
  const double qui = AsBinary(nodes.at("Quiescence"));
  const double prol = AsBinary(nodes.at("Proliferation"));
  const double mig = AsBinary(nodes.at("Migration"));

  out->proliferation_signal = std::clamp(0.1 + 0.9 * prol, 0.0, 1.0);
  out->can_enter_S = (p21 < 0.5) ? 1.0 : 0.0;
  out->can_enter_M = (chk < 0.5) ? 1.0 : 0.0;
  out->apoptosis_hazard = std::max(apop_base, casp * GrnGetNumeric(params, "apoptosis_hazard_scale", 0.15));
  out->necrosis_hazard = std::max(nec_base, hif * GrnGetNumeric(params, "necrosis_hazard_scale", 0.05));
  out->quiescence_hazard = std::max(qui_base, qui * GrnGetNumeric(params, "quiescence_hazard_scale", 0.08));
  out->repair_capacity = std::clamp(0.2 + 0.8 * nrf2, 0.0, 1.0);
  out->antioxidant_capacity = std::clamp(0.2 + 0.8 * nrf2, 0.0, 1.0);
  out->migration_modifier = std::clamp(0.4 + 0.6 * mig, 0.0, 1.5);
  out->ecm_degradation_rate = std::clamp(mmp, 0.0, 1.0);
  out->ecm_deposition_rate = std::clamp((1.0 - mmp) * 0.1, 0.0, 1.0);
  out->secretion_modifier = std::clamp(1.0 + 0.2 * hif, 0.0, 2.0);
}
}  // namespace

inline void BooleanRegulatoryModel::Initialize(const RegulatoryParameters& params) {
  params_ = params;
}

inline void BooleanRegulatoryModel::Update(const RegulatoryInput& input,
                                           RegulatoryState* state,
                                           RegulatoryOutput* output,
                                           double dt) {
  if (!state || !output) return;
  state->time_since_last_update += dt;
  EvaluateCompactDdrEcmNodes(params_, input, &state->node_activity);

  if (GrnGetBool(params_, "enable_soft_boolean", false)) {
    for (auto& kv : state->node_activity) {
      kv.second = Sigmoid((kv.second - 0.5) * 6.0);
    }
  }
  NodesToOutput(params_, state->node_activity, output);
}
// =============================================================================
}  // namespace regulatory
}  // namespace bdm
// =============================================================================
#endif  // _BOOLEAN_GRN_INLINE_H_
