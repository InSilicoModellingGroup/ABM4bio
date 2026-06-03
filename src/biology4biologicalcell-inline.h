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
#ifndef _BIOLOGY4BIOLOGICALCELL_INLINE_H_
#define _BIOLOGY4BIOLOGICALCELL_INLINE_H_
// =============================================================================
#include "./biological_cell.h"
// =============================================================================
inline
void bdm::Biology4BiologicalCell_10::Run(bdm::Agent* a)
{
  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
      // firstly, we should check if cell is inside the simulation domain
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // we check for cell apoptosis
      if (cell->CheckApoptosis())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // simply update the cell age
      cell->IncrementAge();
      // cell produces/consumes substances
      cell->RunBiochemics();
      // intracellular ROS and damage dynamics
      cell->RunIntracellular();
      if (cell->CheckApoptosisByDamage())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      if (cell->CheckTransformation()) return;
      // now check if cell can migrate
      if (cell->CheckMigration())
        {
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
        }
      // then, check if cell can polarize
      cell->CheckPolarization();
      cell->CheckProtrusion();
      // check if cell can grow
      if (cell->CheckGrowth()) return;
      // check if cell can divide (summetrically or unsymmetrically)
      if (cell->CheckTransformationAndDivision()) return;
      if (cell->CheckAsymmetricDivision()) return;
      if (cell->CheckDivision()) return;
      // finally, we check for cell apoptosis due to aging
      if (cell->CheckApoptosisAging())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // ...end of mechanisms list for cell behavior
    }
  else
    ABORT_("an exception is caught");
}
// -----------------------------------------------------------------------------
inline
void bdm::Biology4BiologicalCell_11::Run(bdm::Agent* a)
{
  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
      if (bdm::BiologicalCell::Phase::Ap==cell->GetPhase()) {
        cell->IncrementAge();
        //
        if (cell->CheckAfterApoptosis())
          {
            cell->Set2DeleteProtrusions();
            cell->RemoveBehavior(this);
            cell->RemoveFromSimulation();
            return;
          }
        return;
      }
      // firstly, we should check if cell is inside the simulation domain
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // simply update the cell age
      cell->IncrementAge();
      // cell produces/consumes substances
      cell->RunBiochemics();
      // intracellular ROS and damage dynamics
      cell->RunIntracellular();
      if (cell->CheckApoptosisByDamage())
        {
          cell->SetAge();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      //
      // post-division quiescence is modeled as an early G1 arrest
      if (bdm::BiologicalCell::Phase::G1==cell->GetPhase())
        {
          if (cell->CheckQuiescenceAfterDivision())
            return;
        }
      // we check for cell apoptosis
      if (cell->CheckApoptosis())
        {
          cell->SetAge(); // reset age (time) counter
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      if (cell->CheckTransformation()) return;
      // check if cell can polarize
      if (cell->CheckPolarization())
        {
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
        }
      // now check if cell can migrate
      if (cell->CheckMigration())
        {
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
        }
      // mitosis-like stage: only Di can divide
      if (bdm::BiologicalCell::Phase::Di==cell->GetPhase())
        if (cell->CheckDivision() || cell->CheckAsymmetricDivision())
          {
            cell->SetAge(); // reset age (time) counter
            // after successful division, continue from early-cycle state
            cell->SetPhase(bdm::BiologicalCell::Phase::G1);
            return;
          }
      // biomass growth is restricted to G1
      if (bdm::BiologicalCell::Phase::G1==cell->GetPhase())
        if (cell->CheckGrowth())
          {
            return;
          }
      // finally, we check for cell apoptosis due to aging
      if (bdm::BiologicalCell::Phase::G1==cell->GetPhase()||
          bdm::BiologicalCell::Phase::Sy==cell->GetPhase()||
          bdm::BiologicalCell::Phase::G2==cell->GetPhase()||
          bdm::BiologicalCell::Phase::Di==cell->GetPhase())
        if (cell->CheckApoptosisAging())
          {
            cell->SetAge(); // reset age (time) counter
            cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
            return;
          }
      //
      // default phase progression: G1 -> Sy -> G2 -> Di
      if (bdm::BiologicalCell::Phase::I0==cell->GetPhase() ||
          bdm::BiologicalCell::Phase::Tr==cell->GetPhase())
        cell->SetPhase(bdm::BiologicalCell::Phase::G1);
      else if (bdm::BiologicalCell::Phase::G1==cell->GetPhase())
        cell->SetPhase(bdm::BiologicalCell::Phase::Sy);
      else if (bdm::BiologicalCell::Phase::Sy==cell->GetPhase())
        cell->SetPhase(bdm::BiologicalCell::Phase::G2);
      else if (bdm::BiologicalCell::Phase::G2==cell->GetPhase())
        cell->SetPhase(bdm::BiologicalCell::Phase::Di);
    }
  else
    ABORT_("an exception is caught");
}
// =============================================================================
#endif // _BIOLOGY4BIOLOGICALCELL_INLINE_H_
// =============================================================================
