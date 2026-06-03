// =============================================================================
// Regulatory model interface and concrete backends.
// =============================================================================
#ifndef _REGULATORY_MODEL_H_
#define _REGULATORY_MODEL_H_
// =============================================================================
#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <string>
#include "./regulatory_input.h"
#include "./regulatory_output.h"
#include "./regulatory_state.h"
// =============================================================================
namespace bdm {
namespace regulatory {
// =============================================================================
struct RegulatoryParameters {
  std::map<std::string, double> numeric;
  std::map<std::string, bool> boolean;
  std::string backend = "none";
  int update_interval = 1;
};
// =============================================================================
class RegulatoryModel {
 public:
  virtual ~RegulatoryModel() = default;
  virtual void Initialize(const RegulatoryParameters& params) = 0;
  virtual void Update(const RegulatoryInput& input,
                      RegulatoryState* state,
                      RegulatoryOutput* output,
                      double dt) = 0;
  virtual std::string GetName() const = 0;
};
// =============================================================================
class NullRegulatoryModel : public RegulatoryModel {
 public:
  void Initialize(const RegulatoryParameters& params) override { (void) params; }
  void Update(const RegulatoryInput& input,
              RegulatoryState* state,
              RegulatoryOutput* output,
              double dt) override {
    (void) input;
    (void) dt;
    if (state) {
      state->time_since_last_update += dt;
      state->node_activity.clear();
    }
    if (output) {
      *output = RegulatoryOutput();
    }
  }
  std::string GetName() const override { return "none"; }
};
// =============================================================================
class BooleanRegulatoryModel : public RegulatoryModel {
 public:
  void Initialize(const RegulatoryParameters& params) override;
  void Update(const RegulatoryInput& input,
              RegulatoryState* state,
              RegulatoryOutput* output,
              double dt) override;
  std::string GetName() const override { return "boolean"; }

 private:
  RegulatoryParameters params_;
};
// =============================================================================
class StochasticBooleanRegulatoryModel : public RegulatoryModel {
 public:
  void Initialize(const RegulatoryParameters& params) override;
  void Update(const RegulatoryInput& input,
              RegulatoryState* state,
              RegulatoryOutput* output,
              double dt) override;
  std::string GetName() const override { return "stochastic_boolean"; }

 private:
  RegulatoryParameters params_;
};
// =============================================================================
}  // namespace regulatory
}  // namespace bdm
// =============================================================================
#endif  // _REGULATORY_MODEL_H_
