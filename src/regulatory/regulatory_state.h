// =============================================================================
// Per-cell regulatory internal state.
// =============================================================================
#ifndef _REGULATORY_STATE_H_
#define _REGULATORY_STATE_H_
// =============================================================================
#include <map>
#include <string>
// =============================================================================
namespace bdm {
namespace regulatory {
// =============================================================================
struct RegulatoryState {
  std::map<std::string, double> node_activity;
  double time_since_last_update = 0.0;
  int last_phase = -999;
};
// =============================================================================
}  // namespace regulatory
}  // namespace bdm
// =============================================================================
#endif  // _REGULATORY_STATE_H_
