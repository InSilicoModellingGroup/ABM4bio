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
#ifndef _BIOLOGY_H_
#define _BIOLOGY_H_
// =============================================================================
#include "core/behavior/behavior.h"
// =============================================================================
namespace bdm {
// =============================================================================
struct Biology4BiologicalCell_10 : public bdm::Behavior {
BDM_BEHAVIOR_HEADER(Biology4BiologicalCell_10, bdm::Behavior, 1);
//
public:
  Biology4BiologicalCell_10() { AlwaysCopyToNew(); }
  virtual ~Biology4BiologicalCell_10() {}
  //
  void Initialize(const bdm::NewAgentEvent& event) override {
    Base::Initialize(event);
  }
  //
  void Run(bdm::Agent* a) override;
  //
};
// -----------------------------------------------------------------------------
struct Biology4BiologicalCell_11 : public bdm::Behavior {
BDM_BEHAVIOR_HEADER(Biology4BiologicalCell_11, bdm::Behavior, 1);
//
public:
  Biology4BiologicalCell_11() { AlwaysCopyToNew(); }
  virtual ~Biology4BiologicalCell_11() {}
  //
  void Initialize(const bdm::NewAgentEvent& event) override {
    Base::Initialize(event);
  }
  //
  void Run(bdm::Agent* a) override;
  //
};
// -----------------------------------------------------------------------------
struct Biology4Vessel : public bdm::Behavior {
BDM_BEHAVIOR_HEADER(Biology4Vessel, bdm::Behavior, 1);
//
public:
  Biology4Vessel() { AlwaysCopyToNew(); }
  virtual ~Biology4Vessel() {}
  //
  void Initialize(const bdm::NewAgentEvent& event) override {
    Base::Initialize(event);
  }
  //
  void Run(bdm::Agent* a) override;
  //
};
// -----------------------------------------------------------------------------
struct Biology4CellProtrusion : public bdm::Behavior {
BDM_BEHAVIOR_HEADER(Biology4CellProtrusion, bdm::Behavior, 1);
//
public:
  Biology4CellProtrusion() { AlwaysCopyToNew(); }
  virtual ~Biology4CellProtrusion() {}
  //
  void Initialize(const bdm::NewAgentEvent& event) override {
    Base::Initialize(event);
  }
  //
  void Run(bdm::Agent* a) override;
  //
};
// =============================================================================
} // ...end of namespace
// =============================================================================
#endif // _BIOLOGY_H_
// =============================================================================
