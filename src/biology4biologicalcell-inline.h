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
  // ================================================================
  // Mechanism 11 — Biologically strict generic cell–ECM and
  // microenvironment controller for untreated/control tumour biology.
  //
  // Timestep order:
  //  1. Ap-phase: increment timers, delayed removal, return.
  //  2. Validate position (domain boundaries).
  //  3. Increment global age and per-phase timer.
  //  4. Sample microenvironment (O2, nutrient, ECM, crowding).
  //  5. RunBiochemics: secretion and uptake.
  //  6. RunIntracellular: ROS/damage/DDR dynamics.
  //  7. Apoptosis/necrosis/quiescence hazard evaluation.
  //  8. Central cell-cycle checkpoint controller.
  //  9. Post-division G1 quiescence check.
  // 10. Phenotype transformation.
  // 11. Phase-cycle progression (G1/S, intra-S, G2/M checkpoints; Di dwell).
  // 12. Polarisation (ECM-aware; durotaxis stub via stiffness field).
  // 13. Protrusion / focal adhesion update.
  // 14. Migration (chemotaxis + haptotaxis; ECM density barrier).
  // 15. Biomass growth (G1 + Sy; microenvironment modulation).
  // 16. Division (Di only; Di dwell + checkpoint + volume threshold).
  // 17. Optional ECM remodelling (already called in step 7b via RunECMInteraction).
  // 18. Apoptosis from aging (cycling phases only).
  // ================================================================
  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
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
      // STEP 4 — Ap/Nec state handling
      //   Ap: delayed removal via CheckAfterApoptosis()
      //   Nec: phenotype 0 is treated as terminal/non-cycling in Mechanism 11
      // ================================================================
      if (bdm::BiologicalCell::Phase::Ap == cell->GetPhase())
        {
          if (cell->CheckAfterApoptosis())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
            }
          return;
        }
      if (0 == cell->GetPhenotype())
        {
          return;
        }
      // ================================================================
      // STEP 5 — Sample microenvironment once per timestep
      //   Biological context: the cell continuously senses its local niche
      //   via integrins (ECM), RTKs (growth factors/nutrients), and ROS sensors.
      //   All grid lookups are done here to avoid repeated queries below.
      // ================================================================
      const bdm::BiologicalCell::MicroenvironmentState env =
        cell->SampleMicroenvironment();
      // ================================================================
      // STEP 6 — Secretion/uptake via RunBiochemics
      // ================================================================
      cell->RunBiochemics();
      // ================================================================
      // STEP 7 — Intracellular dynamics via RunIntracellular
      //   Updates: ROS accumulation, antioxidant buffering, DNA damage,
      //   and the DDR cascade (ATM/ATR → CHK1/CHK2 → p53 → p21/Cdc25/CDK).
      // ================================================================
      cell->RunIntracellular();
      // ================================================================
      // STEP 7b — Intracellular regulation module (optional GRN backend)
      //   Mechanism consumes only RegulatoryOutput; backend remains modular.
      // ================================================================
      cell->UpdateRegulatoryModel(env);
      const bdm::regulatory::RegulatoryOutput& grn = cell->GetRegulatoryOutput();
      // ================================================================
      // STEP 8 — Apoptosis / necrosis / quiescence hazard evaluation
      //   Uses rate-to-probability: P(event) = 1 - exp(-hazard * dt)
      //   Biological context: stochastic threshold-crossing events driven
      //   by accumulated microenvironmental and intracellular stress.
      // ================================================================
      const double dt = cell->params()->get<double>("time_step");
      const std::string& CP_name =
        cell->params()->get<std::string>(
          "phenotype_ID/"+std::to_string(cell->GetPhenotype()));
      auto* rg = bdm::Simulation::GetActive()->GetRandom();
      // GRN-driven apoptosis hazard (optional layer)
      if (cell->IsRegulatoryModelActive() && cell->GetCanApoptose()
          && grn.apoptosis_hazard > 0.0)
        {
          const double P_grn_ap = 1.0 - std::exp(-std::max(0.0, grn.apoptosis_hazard) * dt);
          if (rg->Uniform(0.0, 1.0) < P_grn_ap)
            {
              cell->SetAge();
              cell->ResetPhaseAge();
              cell->ResetArrestTime();
              cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
              return;
            }
        }
      //
      // 7a: DNA-damage / intracellular ROS → apoptosis (hard deterministic stop)
      //    Committed when p53 activity exceeds the apoptosis threshold or
      //    aggregate damage exceeds the damage/threshold parameter.
      //    Represents p53-mediated BAX/BCL-2 balance committing to apoptosis.
      if (cell->CheckApoptosisByDamage())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // 7b: ECM adhesion / anoikis and ECM remodelling
      //    Low adhesion → anoikis via BIM activation (integrin/FAK/Src loss).
      //    Anoikis resistance can be parameterised (e.g. BCL-2 overexpression).
      //    Also performs ECM degradation (MMP/mmp_activity) and deposition
      //    (ecm_deposition_rate / caf_ecm_deposition_rate).
      if (cell->RunECMInteraction())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // 7c: Substance-threshold apoptosis (legacy; O2 or nutrient < threshold)
      if (cell->CheckApoptosis())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // 7d: Non-quiescent necrosis hazard (severe O2/nutrient deprivation)
      //    P(necrosis) = 1 - exp(-hazard_rate * dt); active even outside G0.
      //    Biologically: acute energy collapse → lysosomal membrane rupture
      //    → HMGB1/DAMPs release → inflammatory necrosis (oncosis).
      //    Parameters: can_necrose/O2_hazard_rate, can_necrose/nutrient_hazard_rate
      //    (distinct from the quiescence-arrest necrosis path in CheckNecrosis).
      if (cell->GetCanApoptose() && cell->GetPhenotype())
        {
          bool necr_triggered = false;
          // Severe hypoxia hazard (O2 below necrosis threshold)
          if (!necr_triggered &&
              cell->params()->have_parameter<double>(CP_name+"/can_necrose/O2_hazard_rate"))
            {
              const double o2_nec =
                cell->params()->have_parameter<double>(CP_name+"/can_necrose/O2_threshold")
                ? cell->params()->get<double>(CP_name+"/can_necrose/O2_threshold") : 0.0;
              if (env.local_O2 <= o2_nec)
                {
                  const double k_nec =
                    cell->params()->get<double>(CP_name+"/can_necrose/O2_hazard_rate");
                  const double P_nec = 1.0 - std::exp(-std::max(0.0, k_nec) * dt);
                  if (rg->Uniform(0.0, 1.0) < P_nec) necr_triggered = true;
                }
            }
          // Severe nutrient depletion hazard
          if (!necr_triggered &&
              cell->params()->have_parameter<double>(CP_name+"/can_necrose/nutrient_hazard_rate"))
            {
              const double nut_nec =
                cell->params()->have_parameter<double>(CP_name+"/can_necrose/nutrient_threshold")
                ? cell->params()->get<double>(CP_name+"/can_necrose/nutrient_threshold") : 0.0;
              if (env.local_nutrient <= nut_nec)
                {
                  const double k_nec =
                    cell->params()->get<double>(CP_name+"/can_necrose/nutrient_hazard_rate");
                  const double P_nec = 1.0 - std::exp(-std::max(0.0, k_nec) * dt);
                  if (rg->Uniform(0.0, 1.0) < P_nec) necr_triggered = true;
                }
            }
          if (necr_triggered)
            {
              // Transform to necrotic phenotype (ID 0) via CheckNecrosis().
              // If the O2-threshold path is not configured in CheckNecrosis(),
              // force necrosis directly.
              if (cell->CheckNecrosis()) return;
            }
        }
      // 7e: G0/quiescence entry hazard (reversible; crowding, O2, nutrient)
      //    Rate-to-probability: P(G0 entry) = 1 - exp(-k_entry * dt).
      //    Exit: P(G0 → G1) = 1 - exp(-k_exit * dt) when conditions improve.
      //    Biologically: contact inhibition → p27/Kip1 upregulation;
      //    hypoxia → HIF-1α → p21/BNIP3-mediated quiescence.
      //    If entry_hazard_rate is absent, falls back to deterministic
      //    threshold comparison for backward compatibility.
      {
        const bool has_prob_quiescence =
          cell->params()->have_parameter<double>(CP_name+"/quiescence/entry_hazard_rate");
        if (cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_threshold")
            || has_prob_quiescence
            || cell->params()->have_parameter<double>(CP_name+"/quiescence/O2_threshold")
            || cell->params()->have_parameter<double>(CP_name+"/quiescence/nutrient_threshold"))
          {
            const double crowd_entry =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_threshold")
              ? cell->params()->get<double>(CP_name+"/quiescence/crowding_threshold") : 1.0e99;
            const double crowd_exit =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_exit")
              ? cell->params()->get<double>(CP_name+"/quiescence/crowding_exit")
              : crowd_entry * 0.8;
            const double o2_quies =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/O2_threshold")
              ? cell->params()->get<double>(CP_name+"/quiescence/O2_threshold") : -1.0;
            const double nut_quies =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/nutrient_threshold")
              ? cell->params()->get<double>(CP_name+"/quiescence/nutrient_threshold") : -1.0;
            // Microenvironmental stress: is the niche unfavourable for cycling?
            const bool stressed_env =
              (env.local_crowding >= crowd_entry)
              || (o2_quies  >= 0.0 && env.local_O2       < o2_quies)
              || (nut_quies >= 0.0 && env.local_nutrient  < nut_quies)
              || (cell->IsRegulatoryModelActive() && grn.quiescence_hazard > 0.05);
            //
            if (!cell->IsQuiescent() && stressed_env)
              {
                bool enter_quies = false;
                if (has_prob_quiescence)
                  {
                    // Probabilistic entry: P = 1 - exp(-k_entry * dt)
                    const double k_e =
                      cell->params()->get<double>(CP_name+"/quiescence/entry_hazard_rate");
                    enter_quies =
                      (rg->Uniform(0.0, 1.0) < 1.0 - std::exp(-std::max(0.0, k_e) * dt));
                  }
                else
                  enter_quies = true; // deterministic (backward-compatible default)
                if (enter_quies)
                  {
                    cell->SetQuiescent(true);
                    cell->IncrementArrestTime();
                  }
              }
            else if (cell->IsQuiescent())
              {
                const bool conditions_improved =
                  (!stressed_env) || (env.local_crowding < crowd_exit);
                if (conditions_improved)
                  {
                    // G0 → G1 re-entry hazard: P = 1 - exp(-k_exit * dt)
                    bool exit_quies = false;
                    if (cell->params()->have_parameter<double>(
                          CP_name+"/quiescence/exit_hazard_rate"))
                      {
                        const double k_x =
                          cell->params()->get<double>(CP_name+"/quiescence/exit_hazard_rate");
                        exit_quies =
                          (rg->Uniform(0.0, 1.0) < 1.0 - std::exp(-std::max(0.0, k_x) * dt));
                      }
                    else
                      exit_quies = true; // deterministic (backward-compatible)
                    if (exit_quies)
                      {
                        cell->SetQuiescent(false);
                        cell->ResetArrestTime();
                      }
                    else
                      cell->IncrementArrestTime();
                  }
                else
                  cell->IncrementArrestTime();
              }
          }
        //
        // Quiescent (G0) cells: check necrosis, allow migration, then return.
        // Growth, division, and phase progression are suspended in G0.
        if (cell->IsQuiescent())
          {
            // Prolonged G0 under severe hypoxia may lead to necrosis
            // (CheckNecrosis uses O2_threshold + min_arrest_time gate)
            if (cell->CheckNecrosis()) return;
            // G0 cells may still migrate to escape crowded/hypoxic zones.
            // Biologically: amoeboid/mesenchymal switching under metabolic stress;
            // driven by RhoA/ROCK and FAK-independent protrusions.
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
            return; // skip cycling steps while quiescent
          }
      }
      // ================================================================
      // STEP 9 — Central cell-cycle checkpoint controller
      //   Consolidates all gate decisions using the pre-sampled env state.
      //   Returns CellCycleCheckpointState encoding each gate outcome.
      //   Molecular annotations: ATM/ATR–CHK–CDC25–CDK–RB/E2F pathway.
      // ================================================================
      const bdm::BiologicalCell::CellCycleCheckpointState ckpt =
        cell->EvaluateCellCycleCheckpoints(env);
      //
      // Checkpoint-driven forced fate decisions (extreme / irreversible states)
      if (ckpt.should_enter_Ap && cell->GetCanApoptose())
        {
          // Strict hazard form for checkpoint-driven apoptosis in Mechanism 11:
          //   P = 1 - exp(-k_checkpoint_apoptosis * dt)
          // Preferred parameter: can_apoptose/checkpoint_hazard_rate [1/time]
          // Backward compatibility: can_apoptose/checkpoint_probability.
          bool commit_apoptosis = false;
          if (cell->params()->have_parameter<double>(CP_name+"/can_apoptose/checkpoint_hazard_rate"))
            {
              const double k = std::max(0.0,
                cell->params()->get<double>(CP_name+"/can_apoptose/checkpoint_hazard_rate"));
              const double P = 1.0 - std::exp(-k * dt);
              commit_apoptosis = (rg->Uniform(0.0, 1.0) <= P);
            }
          else if (cell->params()->have_parameter<double>(CP_name+"/can_apoptose/checkpoint_probability"))
            {
              const double p_raw =
                cell->params()->get<double>(CP_name+"/can_apoptose/checkpoint_probability");
              const double p = std::max(0.0, std::min(1.0, p_raw));
              if (p >= 1.0)
                {
                  commit_apoptosis = true;
                }
              else
                {
                  const double k = (dt > 0.0 && p > 0.0) ? -std::log1p(-p) / dt : 0.0;
                  const double P = 1.0 - std::exp(-k * dt);
                  commit_apoptosis = (rg->Uniform(0.0, 1.0) <= P);
                }
            }
          else
            {
              // Legacy deterministic behavior when no checkpoint probability is set.
              commit_apoptosis = true;
            }
          if (commit_apoptosis)
            {
              cell->SetAge();
              cell->ResetPhaseAge();
              cell->ResetArrestTime();
              cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
              return;
            }
        }
      if (ckpt.should_enter_Nec)
        {
          // Attempt necrotic transformation; CheckNecrosis() handles
          // the phenotype switch and behavior reset.
          if (cell->CheckNecrosis()) return;
        }
      // ================================================================
      // STEP 10 — Post-division G1 quiescence (early arrest after cell split)
      //   Newly divided daughter cells wait a short period before committing
      //   to another round of G1→S. Prevents immediate re-division.
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
        if (cell->CheckQuiescenceAfterDivision())
          return;
      // ================================================================
      // STEP 11 — Phenotype transformation (e.g. cancer subtype switching)
      //   Checked before phase cycling to allow new phenotype parameters
      //   to govern subsequent steps in this same timestep.
      // ================================================================
      if (cell->CheckTransformation()) return;
      // ================================================================
      // STEP 12 — Cell-cycle phase progression with dwell times + checkpoints
      //   Biological order: I0/Tr → G1 → Sy/S → G2 → Di/M → (divide)
      //
      //   G1 dwell (default 80 steps ≈ 8h at dt=0.1h):
      //     CDK4/6–CyclinD accumulation; RB phosphorylation; E2F activation
      //   G1/S gate: ATM/ATR–CHK1/CHK2–CDC25A–CDK2/CyclinE–RB/E2F
      //     Also gated by: O2, nutrient, ECM density, adhesion, crowding
      //
      //   Sy dwell (default 70 steps ≈ 7h at dt=0.1h):
      //     DNA replication; histone synthesis; PCNA/RPA complex
      //   Intra-S gate: ATR–CHK1–CDC25A ubiquitination → CDK2 inhibition
      //     Blocks S→G2 when replication forks stall or DSBs arise
      //
      //   G2 dwell (default 40 steps ≈ 4h at dt=0.1h):
      //     Mitotic entry preparation; CDK1/CyclinB accumulation
      //   G2/M gate: ATM/ATR–CHK1/CHK2–CDC25C–CDK1/CyclinB inhibition
      //
      //   Di/M dwell (default 0 steps = attempt division immediately):
      //     Spindle assembly; chromosome alignment; APC/C–Cdc20 (SAC)
      //     Set Di_dwell > 0 to enforce M-phase duration.
      // ================================================================
      {
        const int G1_dwell =
          cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G1")
          ? cell->params()->get<int>(CP_name+"/phase_dwell/G1") : 80;
        const int Sy_dwell =
          cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Sy")
          ? cell->params()->get<int>(CP_name+"/phase_dwell/Sy") : 70;
        const int G2_dwell =
          cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G2")
          ? cell->params()->get<int>(CP_name+"/phase_dwell/G2") : 40;
        // Di/M dwell: enforces M-phase duration before division attempt.
        // Default 0 → attempt division every step while in Di (backward-compatible).
        const int Di_dwell =
          cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Di")
          ? cell->params()->get<int>(CP_name+"/phase_dwell/Di") : 0;
        //
        // I0/Tr → G1 (enter cell cycle from initial state or post-transformation)
        if (bdm::BiologicalCell::Phase::I0 == cell->GetPhase() ||
            bdm::BiologicalCell::Phase::Tr == cell->GetPhase())
          {
            cell->SetPhase(bdm::BiologicalCell::Phase::G1);
            cell->ResetPhaseAge();
            cell->ResetArrestTime();
          }
        // G1 → Sy/S: dwell time met AND G1/S checkpoint cleared
        //   Checkpoint result from central controller (STEP 8).
        else if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G1_dwell)
              {
                if (ckpt.can_enter_S)
                  {
                    cell->SetPhase(bdm::BiologicalCell::Phase::Sy);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                  }
                else
                  cell->IncrementArrestTime(); // G1 arrest at G1/S checkpoint
              }
          }
        // Sy/S → G2: S-phase dwell complete AND intra-S checkpoint cleared
        //   Intra-S checkpoint: ATR–CHK1–CDC25A axis. Blocks S→G2 when
        //   DNA damage arises during replication or replication forks stall.
        else if (bdm::BiologicalCell::Phase::Sy == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= Sy_dwell)
              {
                if (!cell->EvaluateIntraSCheckpoint())
                  {
                    cell->SetPhase(bdm::BiologicalCell::Phase::G2);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                  }
                else
                  cell->IncrementArrestTime(); // intra-S arrest
              }
          }
        // G2 → Di/M: dwell time met AND G2/M checkpoint cleared
        //   Checkpoint result from central controller (STEP 8).
        else if (bdm::BiologicalCell::Phase::G2 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G2_dwell)
              {
                if (ckpt.can_enter_M)
                  {
                    cell->SetPhase(bdm::BiologicalCell::Phase::Di);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                  }
                else
                  cell->IncrementArrestTime(); // G2/M arrest
              }
          }
        // Di/M: remains until Di_dwell met; division attempted in STEP 16.
      }
      // ================================================================
      // STEP 13 — Polarisation
      //   ECM fiber orientation (contact guidance): durotaxis is stubbed via
      //   the ecm_stiffness field (sampled in step 4). Stiffer ECM regions
      //   activate focal adhesions → FAK/Src → RhoA/ROCK → actomyosin tension
      //   → YAP/TAZ-mediated polarisation bias.
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
      // ================================================================
      // STEP 14 — Protrusion / focal adhesion update
      //   Filopodial extension driven by Rac1/Cdc42; retraction by RhoA/ROCK.
      //   Protrusion formation is biochemically gated (existing logic).
      // ================================================================
      cell->CheckProtrusion();
      // ================================================================
      // STEP 15 — Migration
      //   Direction optionally combines:
      //     · random motility
      //     · chemotaxis (soluble gradient fields)
      //     · haptotaxis (ECM density gradient; handled by chemotaxis
      //         subsystem when can_migrate/chemotaxis/ECM is configured)
      //     · durotaxis (ECM_stiffness gradient; use chemotaxis/ECM_stiffness)
      //     · contact guidance (ECM_fiber_orientation; stub field)
      //   ECM density barrier: very dense ECM (desmoplastic stroma) blocks
      //   migration. Enable with can_migrate/ECM_barrier parameter.
      //   Crowding penalty is applied inside RunChemotaxis (use_crowding flag).
      // ================================================================
      {
        // ECM density barrier check (optional, inactive by default)
        // Biologically: crosslinked collagen gel (e.g. PDAC stroma, BM)
        // physically blocks cell motility. MMP activity can degrade ECM
        // to re-open migratory paths (ecm_degradation_rate / mmp_activity).
        bool ecm_blocks_migration = false;
        if (cell->params()->have_parameter<double>(CP_name+"/can_migrate/ECM_barrier"))
          {
            const double ecm_bar =
              cell->params()->get<double>(CP_name+"/can_migrate/ECM_barrier");
            if (ecm_bar > 0.0 && env.ecm_density >= ecm_bar)
              ecm_blocks_migration = true;
          }
        if (!ecm_blocks_migration)
          {
            const double mig_gate =
              cell->IsRegulatoryModelActive()
              ? std::clamp(grn.migration_modifier, 0.0, 1.0) : 1.0;
            if (rg->Uniform(0.0, 1.0) <= mig_gate && cell->CheckMigration())
              {
                if (!cell->CheckPositionValidity())
                  {
                    cell->Set2DeleteProtrusions();
                    cell->RemoveBehavior(this);
                    cell->RemoveFromSimulation();
                    return;
                  }
              }
          }
      }
      // ================================================================
      // STEP 16 — Biomass growth (G1 and Sy/S phases)
      //   Growth is phase-gated:
      //     G1:  primary biomass accumulation (protein synthesis, organelles)
      //     Sy:  continued growth during DNA replication
      //   Checkpoint gate: ckpt.can_grow must be true.
      //   Microenvironment modulation: ComputeGrowthModulation(env) ∈ [0,1]
      //     factors in O2 (mTOR/HIF-1α), nutrient (PI3K/Akt/mTOR),
      //     ECM stiffness (YAP/TAZ), crowding (contact inhibition), and
      //     intracellular ROS (p38/JNK growth arrest).
      //   If no modulation parameters are defined, f_env = 1.0 and
      //   CheckGrowth() is called directly (backward-compatible).
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::Sy == cell->GetPhase())
        {
          if (ckpt.can_grow)
            {
              const double f_env = cell->ComputeGrowthModulation(env);
              const double grn_growth =
                cell->IsRegulatoryModelActive()
                ? std::clamp(grn.proliferation_signal, 0.0, 1.0) : 1.0;
              const double f_total = std::clamp(f_env * grn_growth, 0.0, 1.0);
              // Microenvironment modulation gate: skip growth if environment
              // is too poor. When f_env = 1.0 (defaults), gate always passes.
              if (f_total > 0.0 && rg->Uniform(0.0, 1.0) <= f_total)
                {
                  if (cell->CheckGrowth())
                    return; // cell grew; exit to avoid double-acting this step
                }
            }
        }
      // ================================================================
      // STEP 17 — Cell division (Di/M phase only)
      //   Division requires ALL of:
      //     1. Phase == Di (completed G2 → Di transition)
      //     2. Di dwell time met (M-phase duration enforced)
      //     3. Checkpoint permits (ckpt.can_divide from STEP 8)
      //     4. Volume ≥ division_volume_threshold (biomass checkpoint,
      //          if configured; Wee1/CDK1 size sensor)
      //   After division, daughter cell resets to G1, age = 1, phase_age = 0.
      //   Damage / RONS inheritance controlled by daughter_partition_factor.
      //   Biologically: spindle assembly checkpoint (SAC/APC-Cdc20) ensures
      //   correct chromosome alignment before anaphase onset.
      // ================================================================
      if (bdm::BiologicalCell::Phase::Di == cell->GetPhase())
        {
          const int Di_dwell =
            cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Di")
            ? cell->params()->get<int>(CP_name+"/phase_dwell/Di") : 0;
          if (cell->GetPhaseAge() >= Di_dwell && ckpt.can_divide)
            {
              // Optional biomass (volume) threshold for division
              // Biologically: minimum cell mass required for CDK1 activation
              bool volume_ok = true;
              if (cell->params()->have_parameter<double>(
                    CP_name+"/can_divide/division_volume_threshold"))
                {
                  const double r = 0.5 * cell->GetDiameter();
                  const double cell_vol = (4.0/3.0) * bdm::Math::kPi * r * r * r;
                  const double vol_thr =
                    cell->params()->get<double>(
                      CP_name+"/can_divide/division_volume_threshold");
                  if (cell_vol < vol_thr) volume_ok = false;
                }
              if (volume_ok)
                {
                  // Try symmetric division first, then asymmetric division.
                  // CheckDivision / CheckAsymmetricDivision also honour:
                  //   · can_divide/probability, /max, /diameter_cutoff, /max_occupancy
                  //   · intracellular/damage/division_block (hard stop for damage)
                  //   · substance-threshold gates (if configured)
                  // After successful division:
                  //   · Daughter cells are created by BioDynaMo's Divide()
                  //   · Their Initialize() sets phase = I0 (→ G1 next step),
                  //     age = 1, phase_age = 0, arrest_time = 0, quiescent = false
                  //   · Damage / RONS partitioned via daughter_partition_factor
                  if (cell->CheckDivision() || cell->CheckAsymmetricDivision())
                    {
                      // Mother cell transitions: reset to G1 with fresh timers
                      cell->SetAge();
                      cell->ResetPhaseAge();
                      cell->ResetArrestTime();
                      cell->SetQuiescent(false);
                      cell->SetPhase(bdm::BiologicalCell::Phase::G1);
                      return;
                    }
                }
            }
        }
      // ================================================================
      // STEP 18 — ECM remodelling note
      //   ECM degradation (MMP/mmp_activity) and deposition
      //   (ecm_deposition_rate / caf_ecm_deposition_rate) are executed in
      //   STEP 7b via RunECMInteraction(). No additional code needed here.
      //   Extension points for future improvements:
      //     · Fibronectin assembly (autocrine loop)
      //     · Collagen crosslinking (LOX-mediated stiffening)
      //     · Hyaluronan synthesis / CD44-HA signalling
      // ================================================================
      // ================================================================
      // STEP 19 — Apoptosis from aging (active cycling phases only)
      //   Represents cumulative cellular damage and telomere shortening.
      //   Quiescent cells are excluded (they returned early at STEP 7e).
      //   Extension points: senescence (separate flag, not yet implemented),
      //   mitotic catastrophe (separate fate from checkpoint override).
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
      // ================================================================
      // STEP 20 — Statistics
      //   Phase counts (G0/G1/Sy/G2/Di/Ap/Nec), migration events, division
      //   events, and mean microenvironmental fields are exported externally
      //   by save_stats() in ABM4bio.h at each statistics_interval step.
      //   All required information is accessible via cell state fields:
      //     IsQuiescent(), GetPhase(), GetPhenotype(), GetAge(), etc.
      //   Mean O2/nutrient/crowding/ECM are accumulated there per phenotype.
      // ================================================================
      // ...end of Mechanism 11
    }
  else
    ABORT_("an exception is caught");
}
// -----------------------------------------------------------------------------
inline
void bdm::Biology4BiologicalCell_12::Run(bdm::Agent* a)
{
  // ================================================================
  // MECHANISM 12 — CAP/PAM-Induced Cholangiocarcinoma (CCA) Treatment Response
  //
  // Scientific basis:
  //   CAP/PAM exposure → extracellular RONS (H2O2 + NO2- cocktail)
  //   → intracellular ROS/RNS stress → oxidative DNA damage:
  //     • 8-oxoguanine (8-oxoG) via OH· (Fenton/Haber-Weiss from H2O2)
  //     • DNA DSBs detected by γH2AX foci
  //   → DDR activation:
  //     • ATM (DSB sensor) + ATR (ssDNA/replication-fork-stall sensor)
  //     • CHK1 phosphorylation (Ser345; ATR-driven): reported in EGI-1 and HuCCT1
  //     • p53 phosphorylation/accumulation (Ser15; ATM/CHK2-driven)
  //   → Cell-cycle arrest (appears before apoptosis — DDR first, death later):
  //     • G1/S block: p53/p21 → CDK2/CyclinE inhibition
  //     • G2/M block (primary): CHK1/CHK2 → CDC25C inhibition → CDK1/CyclinB inactive
  //     • Intra-S slowing: ATR–CHK1–CDC25A during active replication
  //   → Repair window: NHEJ/HR (DSBs), BER/NER (8-oxoG) during arrest
  //   → If unresolved: apoptosis commitment (delayed, dose-dependent)
  //     • PARP cleavage (caspase-3 substrate) and caspase-3 activation reported
  //     • Annexin V/7-AAD: early/late apoptosis by flow cytometry
  //     • EGI-1 ~72h apoptosis; HuCCT1 ~48h (encoded via phenotype parameters)
  //   → In vivo (CCA xenografts): 8-oxoG IHC + cleaved caspase-3 IHC positive
  //   → Primary hepatocytes: lower CHK1/p53 activation; no PARP cleavage observed
  //     (encoded as high repair_capacity + low k_damage phenotype params, not names)
  //
  // Timestep order (16 biological steps):
  //  1.  Validate domain/position.
  //  2.  Increment age, phase_age, time_since_cap.
  //  3.  Ap/Nec phase handling: delayed removal; return.
  //  4.  Sample local microenvironment (RONS, O2, nutrient, ECM, crowding).
  //  5.  RunBiochemics: secretion and uptake.
  //  6.  RunIntracellular: basic ROS/damage/DDR dynamics.
  //  7.  UpdateCAPIntracellular: RONS-specific layer (RNS, dose integral,
  //      marker proxies: 8-oxoG, γH2AX, pCHK1, p53, PARP, caspase-3,
  //      commitment state, membrane permeability, repair capacity).
  //  8.  Evaluate CAP checkpoint controller (EvaluateCAPCheckpointState).
  //  9.  ECM interaction (anoikis/remodelling).
  // 10.  Hazard-based apoptosis from damage/commitment (CheckApoptosisByDamage).
  // 11.  Hazard-based necrosis (severe RONS/hypoxia).
  // 12.  G0/quiescence entry and exit (crowding/O2/nutrient-driven).
  // 13.  Post-division G1 quiescence check.
  // 14.  Phenotype transformation.
  // 15.  Cell-cycle phase progression with CAP checkpoint gating:
  //       G1/S block (p53/p21), intra-S slowing (ATR–CHK1),
  //       G2/M block (CHK1/CHK2–CDC25C), max-arrest-time → apoptosis.
  //       Division NEVER occurs with unresolved γH2AX / G2/M block.
  // 16.  Polarisation and protrusion (only if viable, not strongly arrested).
  // 17.  Migration (haptotaxis + chemotaxis; ECM barrier check).
  // 18.  Biomass growth (G1 + Sy; microenvironment modulation; damage gate).
  // 19.  Division (Di/M only; strict checkpoint clearance; volume threshold).
  // 20.  Aging-based apoptosis (cumulative damage / telomere shortening).
  // ================================================================
  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
      // ================================================================
      // STEP 1 — Validate position (domain boundaries)
      // ================================================================
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // ================================================================
      // STEP 2 — Increment global age and per-phase timer
      // ================================================================
      cell->IncrementAge();
      cell->IncrementPhaseAge();
      // ================================================================
      // STEP 3 — Ap/Nec phase handling: delayed removal
      //   Ap phase (apoptosis committed): increment timers, delayed removal.
      //   The delayed removal models:
      //     - phagocytic clearance window (Annexin V+ cells persist briefly)
      //     - Western blot / flow cytometry detection window for markers
      //   Necrotic phenotype (ID 0): terminal; handled separately (no cycling).
      // ================================================================
      if (bdm::BiologicalCell::Phase::Ap == cell->GetPhase())
        {
          if (cell->CheckAfterApoptosis())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
            }
          return;
        }
      if (0 == cell->GetPhenotype())
        return; // necrotic/debris: no cycling
      // ================================================================
      // STEP 4 — Sample microenvironment once per timestep
      //   Fields: O2, nutrient, ECM density, ECM stiffness, ECM adhesion,
      //           local crowding, local RONS (ros_internal_ proxy).
      //   All absent fields return safe defaults → backward-compatible.
      // ================================================================
      const bdm::BiologicalCell::MicroenvironmentState env =
        cell->SampleMicroenvironment();
      // ================================================================
      // STEP 5 — Secretion/uptake via RunBiochemics
      //   Handles all substance secretion/consumption including H2O2 and NO2_.
      //   Biochemical signalling (EGF, VEGF, HIF-1α targets) also processed here.
      // ================================================================
      cell->RunBiochemics();
      // ================================================================
      // STEP 6 — Intracellular dynamics via RunIntracellular
      //   Updates: ROS accumulation from H2O2/NO2 uptake (aquaporin-mediated),
      //   antioxidant buffering (GSH/thioredoxin/catalase), DNA damage/repair ODE,
      //   and the DDR cascade:
      //     ATM/ATR → CHK1/CHK2 → p53 → p21/Cdc25/CDK (via UpdateDdrPathway).
      // ================================================================
      cell->RunIntracellular();
      // ================================================================
      // STEP 7 — CAP/RONS-specific intracellular update
      //   Updates: rns_internal_ (NO2-/peroxynitrite tracking),
      //            cap_dose_integral_ (AUC of extracellular RONS),
      //            time_since_cap_ (steps since first treatment exposure),
      //            Marker proxies:
      //              oxidative_damage_8oxoG_proxy_ (8-oxoguanine; BER substrate)
      //              dsb_damage_gammaH2AX_proxy_ (γH2AX foci; DSB marker)
      //              parp_cleavage_proxy_ (execution; activated after caspase-3)
      //              caspase3_activation_proxy_ (executioner; post-commitment)
      //            membrane_permeability_ (lipid peroxidation effects)
      //            repair_capacity_ (oxidative suppression of NER/NHEJ enzymes)
      //            apoptosis_commitment_state_ (progressive BAX/BCL-2 tipping)
      // ================================================================
      cell->UpdateCAPIntracellular();
      // ================================================================
      // STEP 7b — Intracellular regulation module (optional GRN backend)
      // ================================================================
      cell->UpdateRegulatoryModel(env);
      const bdm::regulatory::RegulatoryOutput& grn = cell->GetRegulatoryOutput();
      // ================================================================
      // STEP 8 — Evaluate CAP checkpoint controller
      //   Returns CAPCheckpointState encoding:
      //     can_enter_S:     G1/S gate cleared
      //     can_enter_M:     G2/M gate cleared (critical for CAP response)
      //     intra_s_blocked: Sy→G2 blocked by ATR/CHK1 (S-phase damage)
      //     repair_allowed:  NHEJ/HR/BER permitted (O2/nutrient sufficient)
      //     must_enter_Ap:   commit to apoptosis (commitment threshold or max arrest)
      //     must_enter_Nec:  commit to necrosis (extreme RONS/energy collapse)
      // ================================================================
      const bdm::BiologicalCell::CAPCheckpointState ckpt =
        cell->EvaluateCAPCheckpointState(env);
      //
      // Retrieve commonly used parameters
      const double dt = cell->params()->get<double>("time_step");
      const std::string& CP_name =
        cell->params()->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
      auto* rg = bdm::Simulation::GetActive()->GetRandom();
      // Optional GRN-driven apoptosis hazard
      if (cell->IsRegulatoryModelActive() && cell->GetCanApoptose()
          && grn.apoptosis_hazard > 0.0)
        {
          const double P_grn_ap = 1.0 - std::exp(-std::max(0.0, grn.apoptosis_hazard) * dt);
          if (rg->Uniform(0.0, 1.0) < P_grn_ap)
            {
              cell->SetAge();
              cell->ResetPhaseAge();
              cell->ResetArrestTime();
              cell->SetCapArrestPhase(0);
              cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
              return;
            }
        }
      // ================================================================
      // STEP 9 — ECM interaction (anoikis / ECM remodelling)
      //   Low ECM adhesion → anoikis (BIM/BAD via integrin/FAK loss).
      //   ECM degradation (MMP/protease) and deposition (CAF) also performed.
      //   Uses hazard-rate form for Mechanism 12 (same as Mechanism 11).
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
      // STEP 10 — Damage/commitment-triggered apoptosis
      //   CheckApoptosisByDamage() applies the CAP-specific hazard form:
      //     P(apoptosis) = 1 - exp(-k_damage * dt)
      //   Also incorporates p53 activation threshold from DDR pathway.
      //   Biologically: p53/BAX–BCL-2 imbalance → cytochrome C →
      //   caspase-9 → caspase-3 → PARP cleavage → DNA ladder → Annexin V+.
      //   Immediate apoptosis commitment if CAP checkpoint requires it.
      // ================================================================
      if (ckpt.must_enter_Ap && cell->GetCanApoptose())
        {
          // Stochastic apoptosis commitment: P = 1 - exp(-k_ap * dt)
          // Preferred: cap/apoptosis/commitment_hazard_rate
          // Fallback: intracellular/damage/hazard_rate or deterministic
          double k_ap = 0.0;
          if (cell->params()->have_parameter<double>(CP_name+"/cap/apoptosis/commitment_hazard_rate"))
            {
              k_ap = std::max(0.0,
                cell->params()->get<double>(CP_name+"/cap/apoptosis/commitment_hazard_rate"));
            }
          else if (cell->params()->have_parameter<double>(CP_name+"/intracellular/damage/hazard_rate"))
            {
              k_ap = std::max(0.0,
                cell->params()->get<double>(CP_name+"/intracellular/damage/hazard_rate"));
            }
          else
            {
              k_ap = (dt > 0.0) ? 1.0 / dt : 1.0; // deterministic: always commit
            }
          const double P_ap = 1.0 - std::exp(-k_ap * dt);
          if (rg->Uniform(0.0, 1.0) < P_ap)
            {
              cell->SetAge();
              cell->ResetPhaseAge();
              cell->ResetArrestTime();
              cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
              return;
            }
        }
      // Also apply the standard damage-hazard apoptosis check
      if (cell->CheckApoptosisByDamage())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // Chemical-threshold apoptosis (O2/nutrient-driven; legacy mechanism)
      if (cell->CheckApoptosis())
        {
          cell->SetAge();
          cell->ResetPhaseAge();
          cell->ResetArrestTime();
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // ================================================================
      // STEP 11 — Necrosis hazard
      //   Extreme ROS/RNS + hypoxia → energy collapse → membrane rupture.
      //   Peroxynitrite (ONOO-) destroys mitochondrial complex I and III.
      //   Also handles O2-threshold + prolonged arrest path (original necrosis).
      //   Biologically distinct from apoptosis (passive, uncontrolled).
      // ================================================================
      if (ckpt.must_enter_Nec)
        {
          if (cell->CheckNecrosis()) return;
        }
      // Non-quiescent stochastic necrosis hazard (severe O2/nutrient deprivation)
      if (cell->GetCanApoptose() && cell->GetPhenotype())
        {
          bool necr_triggered = false;
          if (!necr_triggered &&
              cell->params()->have_parameter<double>(CP_name+"/can_necrose/O2_hazard_rate"))
            {
              const double o2_nec =
                cell->params()->have_parameter<double>(CP_name+"/can_necrose/O2_threshold")
                ? cell->params()->get<double>(CP_name+"/can_necrose/O2_threshold") : 0.0;
              if (env.local_O2 <= o2_nec)
                {
                  const double k_nec =
                    cell->params()->get<double>(CP_name+"/can_necrose/O2_hazard_rate");
                  if (rg->Uniform(0.0, 1.0) < 1.0 - std::exp(-std::max(0.0, k_nec) * dt))
                    necr_triggered = true;
                }
            }
          if (!necr_triggered &&
              cell->params()->have_parameter<double>(CP_name+"/can_necrose/nutrient_hazard_rate"))
            {
              const double nut_nec =
                cell->params()->have_parameter<double>(CP_name+"/can_necrose/nutrient_threshold")
                ? cell->params()->get<double>(CP_name+"/can_necrose/nutrient_threshold") : 0.0;
              if (env.local_nutrient <= nut_nec)
                {
                  const double k_nec =
                    cell->params()->get<double>(CP_name+"/can_necrose/nutrient_hazard_rate");
                  if (rg->Uniform(0.0, 1.0) < 1.0 - std::exp(-std::max(0.0, k_nec) * dt))
                    necr_triggered = true;
                }
            }
          if (necr_triggered)
            {
              if (cell->CheckNecrosis()) return;
            }
        }
      // ================================================================
      // STEP 12 — G0/quiescence entry and exit
      //   Reversible quiescence driven by crowding, hypoxia, or nutrient depletion.
      //   G0 cells: skip cycling steps; may still migrate; check necrosis.
      //   CAP-treated G0 cells still accumulate damage and run UpdateCAPIntracellular.
      // ================================================================
      {
        const bool has_prob_quiescence =
          cell->params()->have_parameter<double>(CP_name+"/quiescence/entry_hazard_rate");
        if (cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_threshold")
            || has_prob_quiescence
            || cell->params()->have_parameter<double>(CP_name+"/quiescence/O2_threshold")
            || cell->params()->have_parameter<double>(CP_name+"/quiescence/nutrient_threshold"))
          {
            const double crowd_entry =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_threshold")
              ? cell->params()->get<double>(CP_name+"/quiescence/crowding_threshold") : 1.0e99;
            const double crowd_exit =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/crowding_exit")
              ? cell->params()->get<double>(CP_name+"/quiescence/crowding_exit")
              : crowd_entry * 0.8;
            const double o2_quies =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/O2_threshold")
              ? cell->params()->get<double>(CP_name+"/quiescence/O2_threshold") : -1.0;
            const double nut_quies =
              cell->params()->have_parameter<double>(CP_name+"/quiescence/nutrient_threshold")
              ? cell->params()->get<double>(CP_name+"/quiescence/nutrient_threshold") : -1.0;
            const bool stressed_env =
              (env.local_crowding >= crowd_entry)
              || (o2_quies  >= 0.0 && env.local_O2      < o2_quies)
              || (nut_quies >= 0.0 && env.local_nutrient < nut_quies)
              || (cell->IsRegulatoryModelActive() && grn.quiescence_hazard > 0.05);
            if (!cell->IsQuiescent() && stressed_env)
              {
                bool enter_quies = false;
                if (has_prob_quiescence)
                  {
                    const double k_e =
                      cell->params()->get<double>(CP_name+"/quiescence/entry_hazard_rate");
                    enter_quies =
                      (rg->Uniform(0.0, 1.0) < 1.0 - std::exp(-std::max(0.0, k_e) * dt));
                  }
                else
                  enter_quies = true;
                if (enter_quies) {
                  cell->SetQuiescent(true);
                  cell->IncrementArrestTime();
                }
              }
            else if (cell->IsQuiescent())
              {
                const bool conditions_improved =
                  (!stressed_env) || (env.local_crowding < crowd_exit);
                if (conditions_improved)
                  {
                    bool exit_quies = false;
                    if (cell->params()->have_parameter<double>(CP_name+"/quiescence/exit_hazard_rate"))
                      {
                        const double k_x =
                          cell->params()->get<double>(CP_name+"/quiescence/exit_hazard_rate");
                        exit_quies =
                          (rg->Uniform(0.0, 1.0) < 1.0 - std::exp(-std::max(0.0, k_x) * dt));
                      }
                    else
                      exit_quies = true;
                    if (exit_quies) {
                      cell->SetQuiescent(false);
                      cell->ResetArrestTime();
                    } else
                      cell->IncrementArrestTime();
                  }
                else
                  cell->IncrementArrestTime();
              }
          }
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
      // STEP 13 — Post-division G1 quiescence (daughter re-entry delay)
      //   Prevents immediate re-division of daughter cells.
      //   Also applies in CAP context: daughters may inherit partial damage
      //   and must wait before committing to another round of replication.
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
        if (cell->CheckQuiescenceAfterDivision())
          return;
      // ================================================================
      // STEP 14 — Phenotype transformation
      //   Cancer subtype switching; also handles necrotic transformation.
      //   Mechanism 12 is preserved after transformation:
      //   CheckNecrosis() and CheckTransformation() use the new phenotype's
      //   mechanism_order to assign the correct behavior class, including
      //   Biology4BiologicalCell_12 if the new phenotype also uses Mech 12.
      // ================================================================
      if (cell->CheckTransformation()) return;
      // ================================================================
      // STEP 15 — Cell-cycle phase progression with CAP checkpoints
      //
      //   Dwell times define the baseline phase duration before checkpoint
      //   evaluation. Checkpoints then gate each transition.
      //
      //   CAP-specific additions vs Mechanism 11:
      //   (a) G1/S block uses CAPCheckpointState (damage + DDR + microenv)
      //   (b) Intra-S: damage during Sy slows Sy→G2 (ATR–CHK1–CDC25A)
      //   (c) G2/M block is ABSOLUTE for γH2AX-positive cells (no division
      //       with unresolved DSBs — mitotic catastrophe prevention)
      //   (d) max_arrest_time gate: prolonged arrest → apoptosis commitment
      //       This encodes the delayed-apoptosis biology seen in EGI-1/HuCCT1.
      //
      //   Parameter cap_arrest_phase_ is updated to track which checkpoint
      //   is responsible (for statistics export).
      // ================================================================
      {
        const int G1_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G1")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/G1") : 80;
        const int Sy_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Sy")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/Sy") : 70;
        const int G2_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/G2")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/G2") : 40;
        const int Di_dwell = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Di")
                           ? cell->params()->get<int>(CP_name+"/phase_dwell/Di") : 0;
        const int max_arrest = cell->params()->have_parameter<int>(CP_name+"/cap/repair/max_repair_time")
                             ? cell->params()->get<int>(CP_name+"/cap/repair/max_repair_time")
                             : (cell->params()->have_parameter<int>(CP_name+"/phase_dwell/max_arrest_time")
                                ? cell->params()->get<int>(CP_name+"/phase_dwell/max_arrest_time") : 999999);
        //
        // I0/Tr → G1 (enter cell cycle from initial state or post-transformation)
        if (bdm::BiologicalCell::Phase::I0 == cell->GetPhase() ||
            bdm::BiologicalCell::Phase::Tr == cell->GetPhase())
          {
            cell->SetPhase(bdm::BiologicalCell::Phase::G1);
            cell->ResetPhaseAge();
            cell->ResetArrestTime();
            cell->SetCapArrestPhase(0);
          }
        // G1 → Sy/S: dwell time met AND G1/S checkpoint cleared
        //   CAP: CHK1 phosphorylation → CDC25A degradation → CDK2/CyclinE inactive
        //        p53 → p21 → RB hypophosphorylated → E2F blocked → S-phase genes off
        else if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G1_dwell)
              {
                if (ckpt.can_enter_S)
                  {
                    if (cell->GetCapArrestPhase() == 1)
                      cell->IncrementCapRecoveredCount();
                    cell->SetPhase(bdm::BiologicalCell::Phase::Sy);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                    cell->SetCapArrestPhase(0);
                  }
                else
                  {
                    cell->IncrementArrestTime();
                    cell->SetCapArrestPhase(1); // arrested at G1/S
                    // Prolonged G1 arrest with unrepaired damage → apoptosis commitment
                    if (cell->GetArrestTime() > max_arrest)
                      {
                        cell->SetAge();
                        cell->ResetPhaseAge();
                        cell->ResetArrestTime();
                        cell->SetCapArrestPhase(0);
                        cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
                        return;
                      }
                  }
              }
          }
        // Sy/S → G2: S-phase dwell complete AND intra-S checkpoint cleared
        //   CAP: ATR–CHK1–CDC25A degradation → CDK2 inhibition → S-phase arrest
        //   Moderate damage during replication → slowed Sy→G2 progression.
        else if (bdm::BiologicalCell::Phase::Sy == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= Sy_dwell)
              {
                if (!ckpt.intra_s_blocked && !cell->EvaluateIntraSCheckpoint())
                  {
                    if (cell->GetCapArrestPhase() == 2)
                      cell->IncrementCapRecoveredCount();
                    cell->SetPhase(bdm::BiologicalCell::Phase::G2);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                    cell->SetCapArrestPhase(0);
                  }
                else
                  {
                    cell->IncrementArrestTime();
                    cell->SetCapArrestPhase(2); // arrested at intra-S
                    // Prolonged intra-S arrest → apoptosis
                    if (cell->GetArrestTime() > max_arrest)
                      {
                        cell->SetAge();
                        cell->ResetPhaseAge();
                        cell->ResetArrestTime();
                        cell->SetCapArrestPhase(0);
                        cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
                        return;
                      }
                  }
              }
          }
        // G2 → Di/M: dwell time met AND G2/M checkpoint cleared
        //   ABSOLUTE: γH2AX-positive cells (dsb_damage_gammaH2AX_proxy_ > threshold)
        //   must NEVER enter mitosis. Division with unresolved DSBs causes
        //   mitotic catastrophe → aneuploidy → cell death.
        //   CHK1/CHK2 → CDC25C inhibition keeps CDK1/CyclinB inactive.
        //   This is the most prominently reported CAP/CCA arrest in the literature.
        else if (bdm::BiologicalCell::Phase::G2 == cell->GetPhase())
          {
            if (cell->GetPhaseAge() >= G2_dwell)
              {
                if (ckpt.can_enter_M)
                  {
                    if (cell->GetCapArrestPhase() == 3)
                      cell->IncrementCapRecoveredCount();
                    cell->SetPhase(bdm::BiologicalCell::Phase::Di);
                    cell->ResetPhaseAge();
                    cell->ResetArrestTime();
                    cell->SetCapArrestPhase(0);
                  }
                else
                  {
                    cell->IncrementArrestTime();
                    cell->SetCapArrestPhase(3); // arrested at G2/M
                    // Prolonged G2/M arrest → apoptosis (irreparable DSBs)
                    if (cell->GetArrestTime() > max_arrest)
                      {
                        cell->SetAge();
                        cell->ResetPhaseAge();
                        cell->ResetArrestTime();
                        cell->SetCapArrestPhase(0);
                        cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
                        return;
                      }
                  }
              }
          }
        // Di/M: remains until Di_dwell met; division attempted in STEP 19.
        // Di_dwell = 0 (default) → attempt division every step while in Di.
      }
      // ================================================================
      // STEP 16 — Polarisation (ECM-aware; durotaxis via stiffness field)
      //   CAP-treated cells: arrested cells may depolarise due to cytoskeletal
      //   stress from ROS/RNS (p38/JNK → cofilin activation → F-actin depolymerisation).
      //   Polarisation skip when strongly arrested (arrest_time > skip_polarise_threshold)
      //   is optional; default: always run.
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
      // ================================================================
      // STEP 17 — Protrusion / focal adhesion update
      //   CAP-treated cells: high ROS suppresses lamellipodia (Rac1 oxidation).
      //   Protrusion formation gated by existing biochemical logic.
      // ================================================================
      cell->CheckProtrusion();
      // ================================================================
      // STEP 18 — Migration
      //   ECM density barrier (desmoplastic stroma) checked first.
      //   Chemotaxis + haptotaxis (ECM density gradient) combined.
      //   CAP: arrested cells (cap_arrest_phase_ > 0) may still migrate
      //   (amoeboid switching under stress), but arrest does not block migration here.
      // ================================================================
      {
        bool ecm_blocks_migration = false;
        if (cell->params()->have_parameter<double>(CP_name+"/can_migrate/ECM_barrier"))
          {
            const double ecm_bar =
              cell->params()->get<double>(CP_name+"/can_migrate/ECM_barrier");
            if (ecm_bar > 0.0 && env.ecm_density >= ecm_bar)
              ecm_blocks_migration = true;
          }
        if (!ecm_blocks_migration)
          {
            const double mig_gate =
              cell->IsRegulatoryModelActive()
              ? std::clamp(grn.migration_modifier, 0.0, 1.0) : 1.0;
            if (rg->Uniform(0.0, 1.0) <= mig_gate && cell->CheckMigration())
              {
                if (!cell->CheckPositionValidity())
                  {
                    cell->Set2DeleteProtrusions();
                    cell->RemoveBehavior(this);
                    cell->RemoveFromSimulation();
                    return;
                  }
              }
          }
      }
      // ================================================================
      // STEP 19 — Biomass growth (G1 and Sy phases)
      //   Growth gated by: microenvironment (O2, nutrient, ECM, crowding, ROS)
      //   and checkpoint state (no growth during strong arrest at G1/S or G2/M).
      //   CAP damage gate: high dna_damage_ or ros_internal_ suppresses growth
      //   (p38/JNK → mTOR inhibition → reduced protein synthesis).
      //   Damage gate uses cap/growth/damage_growth_block parameter.
      // ================================================================
      if (bdm::BiologicalCell::Phase::G1 == cell->GetPhase() ||
          bdm::BiologicalCell::Phase::Sy == cell->GetPhase())
        {
          // Growth blocked when checkpoint-arrested at G1/S (cap_arrest_phase_ == 1)
          // Allow growth during G2 arrest window (cells continue protein synthesis).
          const bool arrested_no_grow = (cell->GetCapArrestPhase() == 1);
          if (!arrested_no_grow)
            {
              // CAP damage gate: high damage suppresses growth (mTOR inhibition)
              bool damage_blocks_growth = false;
              if (cell->params()->have_parameter<double>(CP_name+"/intracellular/damage/growth_block"))
                {
                  const double grow_block_thr =
                    cell->params()->get<double>(CP_name+"/intracellular/damage/growth_block");
                  if (cell->GetDNADamage() > grow_block_thr)
                    damage_blocks_growth = true;
                }
              if (cell->params()->have_parameter<double>(CP_name+"/cap/growth/damage_growth_block"))
                {
                  const double grow_block_thr =
                    cell->params()->get<double>(CP_name+"/cap/growth/damage_growth_block");
                  if (cell->GetDNADamage() > grow_block_thr)
                    damage_blocks_growth = true;
                }
              if (!damage_blocks_growth)
                {
                  const double f_env = cell->ComputeGrowthModulation(env);
                  const double grn_growth =
                    cell->IsRegulatoryModelActive()
                    ? std::clamp(grn.proliferation_signal, 0.0, 1.0) : 1.0;
                  const double f_total = std::clamp(f_env * grn_growth, 0.0, 1.0);
                  if (f_total > 0.0 && rg->Uniform(0.0, 1.0) <= f_total)
                    {
                      if (cell->CheckGrowth())
                        return;
                    }
                }
            }
        }
      // ================================================================
      // STEP 20 — Cell division (Di/M phase only)
      //   Division requires ALL of:
      //     1. Phase == Di (completed G2→Di transition: G2/M checkpoint cleared)
      //     2. Di dwell time met (M-phase duration)
      //     3. G2/M checkpoint still clear (dsb_damage_gammaH2AX_proxy_ < block thr)
      //     4. Volume ≥ division_volume_threshold (Wee1/CDK1 size sensor)
      //
      //   Critical CAP constraint: γH2AX proxy check is re-evaluated here
      //   to catch any damage that arose during M-phase itself (Di dwell).
      //   Spindle assembly checkpoint (SAC/APC-Cdc20) gated by existing logic.
      // ================================================================
      if (bdm::BiologicalCell::Phase::Di == cell->GetPhase())
        {
          // Re-read Di_dwell here (was scoped inside step 15 block above)
          const int Di_dwell_div = cell->params()->have_parameter<int>(CP_name+"/phase_dwell/Di")
                               ? cell->params()->get<int>(CP_name+"/phase_dwell/Di") : 0;
          if (cell->GetPhaseAge() >= Di_dwell_div)
            {
              // Re-check G2/M status: no division with active γH2AX / damage
              const double gh2ax_div_thr =
                cell->params()->have_parameter<double>(CP_name+"/cap/checkpoint/G2M/gammaH2AX_block_threshold")
                ? cell->params()->get<double>(CP_name+"/cap/checkpoint/G2M/gammaH2AX_block_threshold") : 0.5;
              const bool division_blocked_by_damage =
                (cell->GetDSBDamageGammaH2AXProxy() > gh2ax_div_thr)
                || bdm::IsMolecularG2MCheckpointBlocked(cell);
              if (!division_blocked_by_damage)
                {
                  // Optional biomass threshold
                  bool volume_ok = true;
                  if (cell->params()->have_parameter<double>(CP_name+"/can_divide/division_volume_threshold"))
                    {
                      const double r = 0.5 * cell->GetDiameter();
                      const double cell_vol = (4.0/3.0) * bdm::Math::kPi * r * r * r;
                      const double vol_thr =
                        cell->params()->get<double>(CP_name+"/can_divide/division_volume_threshold");
                      if (cell_vol < vol_thr) volume_ok = false;
                    }
                  if (volume_ok)
                    {
                      if (cell->CheckDivision() || cell->CheckAsymmetricDivision())
                        {
                          cell->SetAge();
                          cell->ResetPhaseAge();
                          cell->ResetArrestTime();
                          cell->SetQuiescent(false);
                          cell->SetCapArrestPhase(0);
                          cell->SetPhase(bdm::BiologicalCell::Phase::G1);
                          return;
                        }
                    }
                }
            }
        }
      // ================================================================
      // (Di_dwell re-declared locally in STEP 20 to avoid scope issue)
      // ================================================================
      // STEP 21 — Apoptosis from aging (active cycling phases only)
      //   Represents cumulative cellular damage and telomere shortening.
      //   Quiescent cells are excluded (they returned early at STEP 12).
      //   CAP cells with prior damage accumulation may age-apoptose sooner.
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
            cell->SetCapArrestPhase(0);
            cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
            return;
          }
      // ================================================================
      // STEP 22 — Statistics
      //   Phase counts (G0/G1/Sy/G2/Di/Ap/Nec), CAP-specific marker means,
      //   checkpoint arrest counts (G1/S, intra-S, G2/M), dose integral means,
      //   and recovered/viable counts are exported externally by save_stats()
      //   in ABM4bio.h at each statistics_interval step.
      // ================================================================
      // ...end of Mechanism 12
    }
  else
    ABORT_("an exception is caught");
}
// =============================================================================
#endif // _BIOLOGY4BIOLOGICALCELL_INLINE_H_
// =============================================================================
