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
#ifndef _BIOLOGY4VESSEL_INLINE_H_
#define _BIOLOGY4VESSEL_INLINE_H_
// =============================================================================
#include "./vessel.h"
// =============================================================================
inline
void bdm::Biology4Vessel::Run(bdm::Agent* a)
{
  if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
    {
      // firstly, we should check if vessel is inside the simulation domain
      if (!vessel->CheckPositionValidity())
        {
          vessel->RemoveBehavior(this);
          return;
        }
      // simply update the vessel age
      vessel->IncrementAge();
      // vessel produces/consumes substances
      vessel->RunBiochemics();
      // check if vessel can grow
      vessel->CheckGrowth();
      // othewise, check if vessel has sprouted or branched out
      if (vessel->CheckSprouting()) return;
      if (vessel->CheckBranching()) return;
      // ...end of mechanisms list for vessel behavior
    }
  else
    ABORT_("an exception is caught");
}
// =============================================================================
#endif // _BIOLOGY4VESSEL_INLINE_H_
// =============================================================================
