// =============================================================================
// Inline implementation for StochasticBooleanRegulatoryModel.
// =============================================================================
#ifndef _STOCHASTIC_BOOLEAN_GRN_INLINE_H_
#define _STOCHASTIC_BOOLEAN_GRN_INLINE_H_
// =============================================================================
#include "./stochastic_boolean_grn.h"
#include "./boolean_grn-inline.h"
// =============================================================================
namespace bdm {
namespace regulatory {
// =============================================================================
namespace {
inline double SbgGetNumeric(const RegulatoryParameters& p,
                            const std::string& key,
                            double fallback) {
  auto it = p.numeric.find(key);
  return (it == p.numeric.end()) ? fallback : it->second;
}

inline double NodeTransition(double current_state,
                             bool desired_on,
                             double k_on,
                             double k_off,
                             double dt,
                             std::mt19937* gen) {
  if (!gen) return desired_on ? 1.0 : 0.0;
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  const bool current_on = current_state > 0.5;
  if (desired_on && !current_on) {
    const double p = 1.0 - std::exp(-std::max(0.0, k_on) * dt);
    if (uni(*gen) < p) return 1.0;
  }
  if (!desired_on && current_on) {
    const double p = 1.0 - std::exp(-std::max(0.0, k_off) * dt);
    if (uni(*gen) < p) return 0.0;
  }
  return current_on ? 1.0 : 0.0;
}
}  // namespace

inline void StochasticBooleanRegulatoryModel::Initialize(const RegulatoryParameters& params) {
  params_ = params;
}

inline void StochasticBooleanRegulatoryModel::Update(const RegulatoryInput& input,
                                                     RegulatoryState* state,
                                                     RegulatoryOutput* output,
                                                     double dt) {
  if (!state || !output) return;
  state->time_since_last_update += dt;

  std::map<std::string, double> desired;
  EvaluateCompactDdrEcmNodes(params_, input, &desired);

  thread_local std::mt19937 gen(std::random_device{}());
  for (const auto& kv : desired) {
    const std::string on_key = std::string("rate_on_") + kv.first;
    const std::string off_key = std::string("rate_off_") + kv.first;
    const double k_on = SbgGetNumeric(params_, on_key, 1.0);
    const double k_off = SbgGetNumeric(params_, off_key, 0.5);
    const double cur = state->node_activity.count(kv.first) ? state->node_activity[kv.first] : 0.0;
    state->node_activity[kv.first] = NodeTransition(cur, kv.second > 0.5, k_on, k_off, dt, &gen);
  }

  NodesToOutput(params_, state->node_activity, output);
}
// =============================================================================
}  // namespace regulatory
}  // namespace bdm
// =============================================================================
#endif  // _STOCHASTIC_BOOLEAN_GRN_INLINE_H_
