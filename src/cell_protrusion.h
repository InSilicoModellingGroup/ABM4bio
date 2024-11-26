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
#ifndef _CELL_PROTRUSION_H_
#define _CELL_PROTRUSION_H_
// =============================================================================
#include "./global.h"
#include "./biochemical.h"
#include "neuroscience/new_agent_event/split_neurite_element_event.h"
#include "neuroscience/new_agent_event/neurite_branching_event.h"
#include "./obstacles.h"
// =============================================================================
namespace bdm {
// =============================================================================
// class forward declarations:
class BiologicalCell;
// =============================================================================
class CellProtrusion : public bdm::neuroscience::NeuriteElement {
BDM_AGENT_HEADER(CellProtrusion, bdm::neuroscience::NeuriteElement, 1);
//
public:
  CellProtrusion() : bdm::neuroscience::NeuriteElement() {
    cell_ = 0;
    age_ = 1;
    can_branch_ = can_sprout_ = true;
    params_ = 0; // nullify pointer...
  }
  //
  void Initialize(const bdm::NewAgentEvent& event) override {
    Base::Initialize(event);
    if (auto* mother = dynamic_cast<CellProtrusion*>(event.existing_agent))
      {
        if (event.GetUid() == bdm::neuroscience::SplitNeuriteElementEvent::kUid)
          {
            cell_ = mother->GetCell();
            SetAge(); // ...age cannot be inherited
            can_branch_ = mother->GetCanBranch();
            can_sprout_ = mother->GetCanSprout();
            params_ = mother->params_; // copy parameters pointer...
          }
        else if (event.GetUid() == bdm::neuroscience::NeuriteBranchingEvent::kUid)
          {
            cell_ = mother->GetCell();
            SetAge(); // ...age cannot be inherited
            can_branch_ = mother->GetCanBranch();
            can_sprout_ = mother->GetCanSprout();
            params_ = mother->params_; // copy parameters pointer...
          }
        else
          ABORT_("an exception is caught");
      }
  }
  //
  void SetCell(const BiologicalCell* c) const { cell_ = c; }
  const BiologicalCell* GetCell() const { return cell_; }
  //
  void SetAge(unsigned int a =1) { age_ = a; }
  int  GetAge() const { return age_; }
  void IncrementAge() { age_++; }
  //
  int GetTimeRepeats() const;
  //
  void SetCanBranch(bool branches) { can_branch_ = branches; }
  bool GetCanBranch() { return can_branch_; }
  //
  void SetCanSprout(bool sprouts) { can_sprout_ = sprouts; }
  bool GetCanSprout() { return can_sprout_; }
  //
  void SetParametersPointer(Parameters* p) { params_ = p; }
  Parameters* params() const { return params_; }
  //
  bool CheckPositionValidity();
  bool CheckTimeWindow() const;
  void CheckOutGrowth();
  void CheckBranching();
  void CheckState();
  void Set2Delete() const;
//
private:
  // pointer to the cell that this protrusion is associated with
  mutable
  const BiologicalCell* cell_ = 0;
  // cell protrusion age (non-fractional time)
  int age_ = 1;
  // flags to designate (individual) protrusion behaviour
  mutable
  bool can_branch_, can_sprout_;
  // pointer to all simulation parameters
  mutable
  Parameters* params_ = 0;
};
// =============================================================================
} // ...end of namespace
// =============================================================================
#endif // _CELL_PROTRUSION_H_
// =============================================================================
