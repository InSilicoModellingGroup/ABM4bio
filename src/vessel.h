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

// =============================================================================
#ifndef _VESSEL_H_
#define _VESSEL_H_
// =============================================================================
#include "./global.h"
#include "./biochemical.h"
#include "neuroscience/new_agent_event/split_neurite_element_event.h"
#include "neuroscience/new_agent_event/neurite_branching_event.h"
#include "./obstacles.h"
// =============================================================================
namespace bdm {
// =============================================================================
class Vessel : public bdm::neuroscience::NeuriteElement {
BDM_AGENT_HEADER(Vessel, bdm::neuroscience::NeuriteElement, 1);
//
public:
  Vessel() : bdm::neuroscience::NeuriteElement() {
    age_ = 1;
    can_grow_ = can_branch_ = can_sprout_ = false;
    params_ = 0; // nullify pointer...
  }
  //
  void Initialize(const bdm::NewAgentEvent& event) override {
    Base::Initialize(event);
    if (auto* mother = dynamic_cast<Vessel*>(event.existing_agent))
      {
        if (event.GetUid() == bdm::neuroscience::SplitNeuriteElementEvent::kUid)
          {
            age_ = 1; // ...age cannot be inherited
            can_grow_   = mother->GetCanGrow();
            can_branch_ = mother->GetCanBranch();
            can_sprout_ = mother->GetCanSprout();
            vessel_ID_ = mother->GetVesselID();
            params_ = mother->params_; // copy parameters pointer...
          }
        else if (event.GetUid() == bdm::neuroscience::NeuriteBranchingEvent::kUid)
          {
            age_ = 1; // ...age cannot be inherited
            can_grow_   = mother->GetCanGrow();
            can_branch_ = mother->GetCanBranch();
            can_sprout_ = mother->GetCanSprout();
            vessel_ID_ = mother->GetVesselID();
            params_ = mother->params_; // copy parameters pointer...
          }
        else
          ABORT_("an exception is caught");
      }
  }
  //
  void SetAge(unsigned int a =1) const { age_ = a; }
  int  GetAge() const { return age_; }
  void IncrementAge() { age_++; }
  //
  void SetCanGrow(bool grows) { can_grow_ = grows; }
  bool GetCanGrow() const { return can_grow_; }
  //
  void SetCanBranch(bool branches) { can_branch_ = branches; }
  bool GetCanBranch() { return can_branch_; }
  //
  void SetCanSprout(bool sprouts) { can_sprout_ = sprouts; }
  bool GetCanSprout() { return can_sprout_; }
  //
  void SetVesselID(int ID) { vessel_ID_ = ID; }
  int  GetVesselID() { return vessel_ID_; }
  //
  void SetParametersPointer(Parameters* p) { params_ = p; }
  Parameters* params() const { return params_; }
  //
  void RunBiochemics();
  bool CheckPositionValidity();
  bool CheckGrowth();
  bool CheckBranching();
  bool CheckSprouting();
  void ScanAge() const;
  //
//
private:
  //
  void CheckSproutingValidity(double& extend_rate, bdm::Double3& direction);
  //
//
private:
  // vessel age (non-fractional time)
  mutable
  int age_ = 1;
  // flags to designate (individual) vessel behaviour
  bool can_grow_, can_branch_, can_sprout_;
  // ID of the vessel to control branching, sprouting, bifurcations
  int vessel_ID_;
  // pointer to all simulation parameters
  mutable
  Parameters* params_ = 0;
};
// =============================================================================
} // ...end of namespace
// =============================================================================
#endif // _VESSEL_H_
// =============================================================================
