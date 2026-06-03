// =============================================================================
// Regulatory input state passed from BiologicalCell + local microenvironment.
// =============================================================================
#ifndef _REGULATORY_INPUT_H_
#define _REGULATORY_INPUT_H_
// =============================================================================
namespace bdm {
namespace regulatory {
// =============================================================================
struct RegulatoryInput {
  double oxygen = 1.0;
  double nutrient = 1.0;
  double extracellular_rons = 0.0;
  double intracellular_ros = 0.0;
  double dna_damage = 0.0;
  double ecm_density = 1.0;
  double ecm_stiffness = 1.0;
  double adhesion_signal = 1.0;
  double crowding = 0.0;
  double growth_factor = 0.0;
  double tgfb = 0.0;
  double inflammatory_signal = 0.0;
  int phase = 0;
  int phenotype_id = 0;
};
// =============================================================================
}  // namespace regulatory
}  // namespace bdm
// =============================================================================
#endif  // _REGULATORY_INPUT_H_
