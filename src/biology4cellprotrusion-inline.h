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
#ifndef _BIOLOGY4CELLPROTRUSION_INLINE_H_
#define _BIOLOGY4CELLPROTRUSION_INLINE_H_
// =============================================================================
#include "./cell_protrusion.h"
// =============================================================================
inline
void bdm::Biology4CellProtrusion::Run(bdm::Agent* a)
{
  if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
    {
      // check if corresponding cell (to this protrusion) exists in simulation
      auto* cell = protrusion->GetCell();
      if (nullptr==cell)
        {
          protrusion->RemoveBehavior(this);
          protrusion->RemoveFromSimulation();
          return;
        }
      // simply update the cell protrusion age
      protrusion->IncrementAge();
      // check if cell age is within appropriate time-window to allow development
      if (!protrusion->CheckTimeWindow()) return;
      //
      const int n_repeat = protrusion->GetTimeRepeats();
      for (int repeat=0; repeat<n_repeat; repeat++)
        {
          // firstly, we should check if protrusion is inside the simulation domain
          if (!protrusion->CheckPositionValidity())
            {
              protrusion->RemoveBehavior(this);
              return;
            }
          // check if cell protrusion outgrows or/and branches
          protrusion->CheckOutGrowth();
          protrusion->CheckBranching();
        }
      //
      protrusion->CheckState();
    }
  else
    ABORT_("an exception is caught");
}
// =============================================================================
#endif // _BIOLOGY4CELLPROTRUSION_INLINE_H_
// =============================================================================
