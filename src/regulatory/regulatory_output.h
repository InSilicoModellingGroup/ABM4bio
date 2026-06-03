// =============================================================================
// Regulatory outputs consumed by Mechanism 11 / 12 controllers.
// =============================================================================
#ifndef _REGULATORY_OUTPUT_H_
#define _REGULATORY_OUTPUT_H_
// =============================================================================
namespace bdm {
namespace regulatory {
// =============================================================================
struct RegulatoryOutput {
  double proliferation_signal = 1.0;
  double can_enter_S = 1.0;
  double can_enter_M = 1.0;
  double apoptosis_hazard = 0.0;
  double necrosis_hazard = 0.0;
  double quiescence_hazard = 0.0;
  double repair_capacity = 1.0;
  double antioxidant_capacity = 1.0;
  double migration_modifier = 1.0;
  double ecm_degradation_rate = 0.0;
  double ecm_deposition_rate = 0.0;
  double secretion_modifier = 1.0;
};
// =============================================================================
}  // namespace regulatory
}  // namespace bdm
// =============================================================================
#endif  // _REGULATORY_OUTPUT_H_
