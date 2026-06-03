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
      // ================================================================
      // STEP 1 — Ap-phase: delayed removal after committed apoptosis
      // ================================================================
      if (bdm::BiologicalCell::Phase::Ap == cell->GetPhase())
        {
          cell->IncrementAge();
          cell->IncrementPhaseAge();
          if (cell->CheckAfterApoptosis())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
            }
          return;
        }
      // ================================================================
      // STEP 2 — Validate position (domain boundaries)
      // ================================================================
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // ================================================================
      // STEP 3 — Increment global age and per-phase timer
      // ================================================================
      cell->IncrementAge();
      cell->IncrementPhaseAge();
      // ================================================================
      // STEP 4+5 — Secretion/uptake and intracellular dynamics
      // ================================================================
      cell->RunBiochemics();
      cell->RunIntracellular();
      // ================================================================
      // STEP 6 — ECM interaction: adhesion sensing, remodelling, anoikis
      // ================================================================
      if (cell->RunECMInteraction())
        {
          // anoikis triggered by low ECM adhesion
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // ================================================================
      // STEP 7 — RONS/damage-triggered apoptosis
      // ================================================================
      if (cell->CheckApoptosisByDamage())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // ================================================================
      // STEP 8 — Post-division quiescence (early G1 arrest)
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
        if (cell->CheckQuiescenceAfterDivision())
          return;
      // ================================================================
      // STEP 9 — Chemical-threshold apoptosis (O2/nutrient-driven)
      // ================================================================
      if (cell->CheckApoptosis())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // ================================================================
      // STEP 10 — Phenotype transformation (e.g. cancer -> necrotic)
      // ================================================================
      if (cell->CheckTransformation()) return;
      // ================================================================
      // STEP 11 — G0 quiescence: crowding/nutrient-driven reversible arrest
      // ================================================================
      {
        const std::string& CP_name =
          cell->params()->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
        //
        if (cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_threshold"))
          {
            const double influence_ratio =
              cell->params()->have_parameter<double>(CP_name+"/can_divide/influence_ratio")
              ? cell->params()->get<double>(CP_name+"/can_divide/influence_ratio") : 2.0;
            const double occ = cell->ComputeLocalOccupancyRatio(cell->GetPosition(), influence_ratio);
            const double crowd_entry = cell->params()->get<double>(CP_name+"/quiescence/crowding_threshold");
            const double crowd_exit  =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_exit")
              ? cell->params()->get<double>(CP_name+"/quiescence/crowding_exit")
              : crowd_entry * 0.8;
            //
            if (!cell->IsQuiescent() && occ >= crowd_entry)
              {
                cell->SetQuiescent(true);
                cell->IncrementArrestTime();
              }
            else if (cell->IsQuiescent())
              {
                if (occ < crowd_exit)
                  {
                    cell->SetQuiescent(false);
                    cell->ResetArrestTime();
                  }
                else
                  cell->IncrementArrestTime();
              }
          }
        //
        if (cell->IsQuiescent())
          {
            // Necrosis may occur during prolonged quiescence under severe hypoxia
            if (cell->CheckNecrosis()) return;
            // Quiescent cells may still migrate to escape crowded/hypoxic areas
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
            // Skip growth, division, and phase progression while quiescent
            return;
          }
      }
      // ================================================================
      // STEP 12 — Phase-cycle progression with dwell times + checkpoints
      //  Biological order: I0/Tr -> G1 -> Sy -> G2 -> Di -> (divide)
      //  G1/S gate: damage, hypoxia, sparse ECM, crowding
      //  G2/M gate: damage, hypoxia
      // ================================================================
      {
        const std::string& CP_name =
          cell->params()->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
        // Configurable per-phase dwell times (number of time steps).
        // Defaults correspond to typical human cancer cell cycle at dt = 0.1 h:
        //   G1 ≈ 8 h → 80 steps, S ≈ 7 h → 70, G2 ≈ 4 h → 40, M ≈ 1.5 h → 15
        const int G1_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G1")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/G1") : 80;
        const int Sy_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Sy")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/Sy") : 70;
        const int G2_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G2")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/G2") : 40;
        //
        // I0/Tr -> G1 (entry into cell cycle)
        if (bdm::BiologicalCell::Phase::I0 == cell->GetPhase() ||
            bdm::BiologicalCell::Phase::Tr == cell->GetPhase())
          {
            cell->SetPhase(bdm::BiologicalCell::Phase::G1);
            cell->ResetPhaseAge();
            cell->ResetArrestTime();
          }
        // G1 -> Sy: dwell time met AND G1/S checkpoint cleared
        else if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G1_dwell)
              {
                if (!cell->EvaluateG1SCheckpoint())
                  {
                    cell->SetPhase(bdm::BiologicalCell::Phase::Sy);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                  }
                else
                  cell->IncrementArrestTime(); // G1 arrest at checkpoint
              }
          }
        // Sy -> G2: S-phase dwell complete
        else if (bdm::BiologicalCell::Phase::Sy == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= Sy_dwell)
              {
                cell->SetPhase(bdm::BiologicalCell::Phase::G2);
                cell->ResetPhaseAge();
                cell->ResetArrestTime();
              }
          }
        // G2 -> Di: dwell time met AND G2/M checkpoint cleared
        else if (bdm::BiologicalCell::Phase::G2 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G2_dwell)
              {
                if (!cell->EvaluateG2MCheckpoint())
                  {
                    cell->SetPhase(bdm::BiologicalCell::Phase::Di);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                  }
                else
                  cell->IncrementArrestTime(); // G2/M arrest
              }
          }
        // Di remains Di until division is attempted (Step 15)
      }
      // ================================================================
      // STEP 13 — Polarization and protrusions
      // ================================================================
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
      cell->CheckProtrusion();
      // ================================================================
      // STEP 14 — Migration (ECM haptotaxis + chemotaxis, any viable phase)
      // ================================================================
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
      // ================================================================
      // STEP 15 — Biomass growth: restricted to G1, O2/nutrient gated
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
        if (cell->CheckGrowth())
          return;
      // ================================================================
      // STEP 16 — Division: only from Di/M after G2/M checkpoint cleared
      // ================================================================
      if (bdm::BiologicalCell::Phase::Di == cell->GetPhase())
        {
          if (cell->CheckDivision() || cell->CheckAsymmetricDivision())
            {
              cell->SetAge();
              cell->ResetPhaseAge();
              cell->ResetArrestTime();
              cell->SetQuiescent(false);
              cell->SetPhase(bdm::BiologicalCell::Phase::G1);
              return;
            }
        }
      // ================================================================
      // STEP 17 — Apoptosis from aging (active cycling phases only)
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::Sy == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::G2 == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::Di == cell->GetPhase())
        if (cell->CheckApoptosisAging())
          {
            cell->SetAge();
            cell->ResetPhaseAge();
            cell->ResetArrestTime();
            cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
            return;
          }
      // ...end of Mechanism 11
    }
  else
    ABORT_("an exception is caught");
}
// -----------------------------------------------------------------------------
inline
void bdm::Biology4BiologicalCell_12::Run(bdm::Agent* a)
{
  // Mechanism 12: Mechanism 11 base + CAP/RONS treatment-response logic.
  // Adds: strict G2/M arrest with max-arrest-time-based apoptosis escalation.
  // All other steps are identical to Mechanism 11.
  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
      // ================================================================
      // STEP 1 — Ap-phase: delayed removal after committed apoptosis
      // ================================================================
      if (bdm::BiologicalCell::Phase::Ap == cell->GetPhase())
        {
          cell->IncrementAge();
          cell->IncrementPhaseAge();
          if (cell->CheckAfterApoptosis())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
            }
          return;
        }
      // ================================================================
      // STEP 2 — Validate position (domain boundaries)
      // ================================================================
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // ================================================================
      // STEP 3 — Increment global age and per-phase timer
      // ================================================================
      cell->IncrementAge();
      cell->IncrementPhaseAge();
      // ================================================================
      // STEP 4+5 — Secretion/uptake and intracellular CAP/RONS dynamics
      // RunIntracellular() handles H2O2/NO2 uptake, ROS accumulation,
      // antioxidant buffering, and DNA damage/repair — all CAP-specific
      // parameters are read from the CSV and applied per cell line.
      // ================================================================
      cell->RunBiochemics();
      cell->RunIntracellular();
      // ================================================================
      // STEP 6 — ECM interaction: adhesion sensing, remodelling, anoikis
      // ================================================================
      if (cell->RunECMInteraction())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // ================================================================
      // STEP 7 — RONS/damage-triggered apoptosis (high damage path)
      // ================================================================
      if (cell->CheckApoptosisByDamage())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // ================================================================
      // STEP 8 — Post-division quiescence (early G1 arrest)
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
        if (cell->CheckQuiescenceAfterDivision())
          return;
      // ================================================================
      // STEP 9 — Chemical-threshold apoptosis (O2/nutrient-driven)
      // ================================================================
      if (cell->CheckApoptosis())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // ================================================================
      // STEP 10 — Phenotype transformation (e.g. cancer -> necrotic)
      // ================================================================
      if (cell->CheckTransformation()) return;
      // ================================================================
      // STEP 11 — G0 quiescence + necrosis (same as Mechanism 11)
      // ================================================================
      {
        const std::string& CP_name =
          cell->params()->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
        //
        if (cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_threshold"))
          {
            const double influence_ratio =
              cell->params()->have_parameter<double>(CP_name+"/can_divide/influence_ratio")
              ? cell->params()->get<double>(CP_name+"/can_divide/influence_ratio") : 2.0;
            const double occ = cell->ComputeLocalOccupancyRatio(cell->GetPosition(), influence_ratio);
            const double crowd_entry = cell->params()->get<double>(CP_name+"/quiescence/crowding_threshold");
            const double crowd_exit  =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_exit")
              ? cell->params()->get<double>(CP_name+"/quiescence/crowding_exit")
              : crowd_entry * 0.8;
            //
            if (!cell->IsQuiescent() && occ >= crowd_entry)
              {
                cell->SetQuiescent(true);
                cell->IncrementArrestTime();
              }
            else if (cell->IsQuiescent())
              {
                if (occ < crowd_exit)
                  {
                    cell->SetQuiescent(false);
                    cell->ResetArrestTime();
                  }
                else
                  cell->IncrementArrestTime();
              }
          }
        //
        if (cell->IsQuiescent())
          {
            if (cell->CheckNecrosis()) return;
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
            return;
          }
      }
      // ================================================================
      // STEP 12 — Phase-cycle progression (Mechanism 12 version)
      //  Adds max_arrest_time gate: prolonged checkpoint arrest -> Ap
      //  This models CHK1/p53 activation leading to apoptosis after
      //  irreparable DNA damage (documented in CAP/CCA studies).
      // ================================================================
      {
        const std::string& CP_name =
          cell->params()->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
        const int G1_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G1")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/G1") : 80;
        const int Sy_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Sy")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/Sy") : 70;
        const int G2_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G2")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/G2") : 40;
        // Maximum time steps a cell may remain checkpoint-arrested before
        // the damage is considered irreparable and apoptosis is committed.
        const int max_arrest = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/max_arrest_time")
                             ? cell->params()->get<int>(CP_name+"/phase_dwell/max_arrest_time") : 999999;
        //
        // I0/Tr -> G1
        if (bdm::BiologicalCell::Phase::I0 == cell->GetPhase() ||
            bdm::BiologicalCell::Phase::Tr == cell->GetPhase())
          {
            cell->SetPhase(bdm::BiologicalCell::Phase::G1);
            cell->ResetPhaseAge();
            cell->ResetArrestTime();
          }
        // G1 -> Sy: dwell + G1/S checkpoint (+ max-arrest escalation)
        else if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G1_dwell)
              {
                if (!cell->EvaluateG1SCheckpoint())
                  {
                    cell->SetPhase(bdm::BiologicalCell::Phase::Sy);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                  }
                else
                  {
                    cell->IncrementArrestTime();
                    // Prolonged G1 arrest with unrepaired damage -> apoptosis
                    if (cell->GetArrestTime() > max_arrest)
                      {
                        cell->SetAge();
                        cell->ResetPhaseAge();
                        cell->ResetArrestTime();
                        cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
                        return;
                      }
                  }
              }
          }
        // Sy -> G2
        else if (bdm::BiologicalCell::Phase::Sy == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= Sy_dwell)
              {
                cell->SetPhase(bdm::BiologicalCell::Phase::G2);
                cell->ResetPhaseAge();
                cell->ResetArrestTime();
              }
          }
        // G2 -> Di: dwell + strict G2/M checkpoint (+ max-arrest escalation)
        else if (bdm::BiologicalCell::Phase::G2 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G2_dwell)
              {
                if (!cell->EvaluateG2MCheckpoint())
                  {
                    cell->SetPhase(bdm::BiologicalCell::Phase::Di);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                  }
                else
                  {
                    cell->IncrementArrestTime();
                    // Prolonged G2/M arrest -> apoptosis (irreparable damage)
                    if (cell->GetArrestTime() > max_arrest)
                      {
                        cell->SetAge();
                        cell->ResetPhaseAge();
                        cell->ResetArrestTime();
                        cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
                        return;
                      }
                  }
              }
          }
        // Di remains Di until division is attempted (Step 15)
      }
      // ================================================================
      // STEP 13 — Polarization and protrusions
      // ================================================================
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
      cell->CheckProtrusion();
      // ================================================================
      // STEP 14 — Migration (ECM haptotaxis + chemotaxis, any viable phase)
      // ================================================================
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
      // ================================================================
      // STEP 15 — Biomass growth: G1, nutrients sufficient
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
        if (cell->CheckGrowth())
          return;
      // ================================================================
      // STEP 16 — Division: Di only, after strict G2/M checkpoint cleared
      // ================================================================
      if (bdm::BiologicalCell::Phase::Di == cell->GetPhase())
        {
          if (cell->CheckDivision() || cell->CheckAsymmetricDivision())
            {
              cell->SetAge();
              cell->ResetPhaseAge();
              cell->ResetArrestTime();
              cell->SetQuiescent(false);
              cell->SetPhase(bdm::BiologicalCell::Phase::G1);
              return;
            }
        }
      // ================================================================
      // STEP 17 — Apoptosis from aging
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::Sy == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::G2 == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::Di == cell->GetPhase())
        if (cell->CheckApoptosisAging())
          {
            cell->SetAge();
            cell->ResetPhaseAge();
            cell->ResetArrestTime();
            cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
            return;
          }
      // ...end of Mechanism 12
    }
  else
    ABORT_("an exception is caught");
}
// =============================================================================
#endif // _BIOLOGY4BIOLOGICALCELL_INLINE_H_
// =============================================================================
