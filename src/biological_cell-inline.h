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
#ifndef _BIOLOGICAL_CELL_INLINE_H_
#define _BIOLOGICAL_CELL_INLINE_H_
// =============================================================================
#include "intracellular/ddr_pathway.h"
inline
void bdm::BiologicalCell::RunBiochemics()
{
  // by design only viable (non-necrotic) cells could secrete biochemicals
  if (!this->GetPhenotype()) return;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0,
               deltaCOORD = maxCOORD - meanCOORD;
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const double dt = this->params()->get<double>("time_step");
  const int mo = this->params()->get<int>(CP_name+"/mechanism_order");
  //
  // cell spatial coordinates
  const bdm::Double3 xyz = this->GetPosition();
  // exit function if cell resides at the edge of the simulation domain,
  // but first identify the mode / dimension of the simulation domain
  if (this->params()->get<bool>("simulation_domain_is_polar"))
    {
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        {
          if ( sqrt(pow2(xyz[0]-meanCOORD)
                   +pow2(xyz[1]-meanCOORD)) > deltaCOORD-tol ) return;
        }
      else
        {
          if ( sqrt(pow2(xyz[0]-meanCOORD)
                   +pow2(xyz[1]-meanCOORD)
                   +pow2(xyz[2]-meanCOORD)) > deltaCOORD-tol ) return;
        }
    }
  else
    {
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol ) return;
        }
      else
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol ||
               xyz[2] < minCOORD+tol || xyz[2] > maxCOORD-tol ) return;
        }
    }
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, xyz, tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      Biochemical BC_id =
        static_cast<Biochemical>(dg->GetContinuumId()); // biochemical ID
      //
      // skip following calculations for radiation!!!
      if ( Biochemical::RAD == BC_id ) continue;
      //
      const double concentration = GetInterpolatedValue(dg, xyz, this->params());
      //
      // --- Michaelis-Menten kinetics model (optional, per-substance) ---
      // When michaelis_menten_model = true, the standard net_balance pathway
      // is REPLACED by a concentration-dependent rate:
      //   R = Vmax * C / (Km + C)
      // The sign of Vmax determines the direction:
      //   Vmax < 0 → consumption,  Vmax > 0 → production.
      // This is a self-contained model: no net_balance, saturation, or
      // dependency parameters are needed.
      if (this->params()->have_parameter<bool>(CP_name+"/"+BC_name+"/secretion/michaelis_menten_model") &&
          this->params()->get<bool>(CP_name+"/"+BC_name+"/secretion/michaelis_menten_model"))
        {
          const double Vmax = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Vmax");
          const double Km   = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Km");
          if (concentration > 0.0 && Km > 0.0)
            {
              double mm_rate = Vmax * concentration / (Km + concentration);
              // apply stochastic variability if defined
              if (this->params()->have_parameter<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/std") &&
                  this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/std") > 0.0)
                {
                  const double mm_std = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/std");
                  mm_rate *= rg->Uniform(1.0 - mm_std, 1.0 + mm_std);
                }
              // apply the rate to the grid
              if (mm_rate > 0.0)
                { // production
                  dg->ChangeConcentrationBy(xyz, mm_rate);
                }
              else
                { // consumption
                  if (concentration + mm_rate > 0.0) dg->ChangeConcentrationBy(xyz, mm_rate);
                  else                               dg->ChangeConcentrationBy(xyz, -concentration);
                }
            }
          // done with MM for this substance — skip the net_balance pathway
          continue;
        }
      //
      // --- Standard net_balance pathway ---
      if (! this->params()->have_parameter<double>(CP_name+"/"+BC_name+"/secretion/net_balance"))
        continue;
      //
      // parameters that modulate biochemical cue secretion (production or consumption)
      const double BC_stdev = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std")<=0.0 ? 1.0 :
                              rg->Uniform(1.0-this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std"),
                                          1.0+this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std"));
      double net_balance = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance") * BC_stdev;
      //
      // skip subsequent calculations if net balance of this biochemical cue secretion is
      // equal to absolute zero!!!
      if (0.0 == net_balance) continue;
      //
      const double saturation = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/saturation");
      //
      if (! this->params()->get<bool>(CP_name+"/"+BC_name+"/secretion/dependency"))
        {
          if (net_balance > 0.0)
            {
              // increase concentration
              if ( saturation > 0.0 )
                {
                  if (concentration>saturation) ; // ... do nothing!
                  else                          dg->ChangeConcentrationBy(xyz, net_balance);
                }
              else
                { // no saturation effects are accounted in
                  dg->ChangeConcentrationBy(xyz, net_balance);
                }
              // increased concentration
            }
          else
            { // ... net_balance < 0.0
              // decrease concentration
              if ( concentration > 0.0 )
                {
                  if (concentration+net_balance>0.0) dg->ChangeConcentrationBy(xyz, net_balance);
                  else                               dg->ChangeConcentrationBy(xyz, -concentration);
                }
              // decreased concentration
            }
          // ...continue to next substance
          continue;
        }
      //
      // check for positive or negative feedback loop from other substances
      // iterate for all substances - simply ignore for itself ;)
      for ( std::vector<std::string>::const_iterator
            cj=substances.begin(); cj!=substances.end(); cj++ )
        {
          // skip if the same substance...
          if (ci==cj) continue;
          // access the BioDynaMo diffusion grid
          auto* dg_other = rm->GetDiffusionGrid(*cj);
          const std::string BC_other_name = dg_other->GetContinuumName(); // biochemical name
          //
          const double concentration_other = GetInterpolatedValue(dg_other, xyz, this->params()),
                       threshold_other = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/"+BC_other_name+"/threshold");
          // check if other substances regulate secretion of this substance...
          if ( ( threshold_other > 0.0 && concentration_other > +threshold_other ) ||
               ( threshold_other < 0.0 && concentration_other < -threshold_other ) )
            {
              //
              if (net_balance > 0.0)
                {
                  // increase concentration
                  if ( saturation > 0.0 )
                    {
                      if (concentration>saturation) ; // ... do nothing!
                      else                          dg->ChangeConcentrationBy(xyz, net_balance);
                    }
                  else
                    { // no saturation effects are accounted in
                      dg->ChangeConcentrationBy(xyz, net_balance);
                    }
                  // increased concentration
                }
              else
                { // ... net_balance < 0.0
                  // decrease concentration
                  if ( concentration > 0.0 )
                    {
                      if (concentration+net_balance>0.0) dg->ChangeConcentrationBy(xyz, net_balance);
                      else                               dg->ChangeConcentrationBy(xyz, -concentration);
                    }
                  // decreased concentration
                }
              //
            }
          //...end of other substances loop
        }
      //...end of substances loop
    }
  //...end of cell biochemics
}
// -----------------------------------------------------------------------------
inline
void bdm::BiologicalCell::RunIntracellular()
{
  // only viable cells maintain intracellular dynamics
  if (!this->GetPhenotype()) return;
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // timestep
  const double dt = this->params()->get<double>("time_step");
  // phenotype-specific namespace
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  // spatial location
  const bdm::Double3 xyz = this->GetPosition();
  // parameters (use safe defaults if not provided)
  const double k_uptake_H2O2 = this->params()->have_parameter<double>(CP_name+"/intracellular/uptake/H2O2")
                             ? this->params()->get<double>(CP_name+"/intracellular/uptake/H2O2") : 0.0;
  const double k_uptake_NO2  = this->params()->have_parameter<double>(CP_name+"/intracellular/uptake/NO2_")
                             ? this->params()->get<double>(CP_name+"/intracellular/uptake/NO2_")  : 0.0;
  antioxidant_capacity_ = this->params()->have_parameter<double>(CP_name+"/intracellular/antioxidant/capacity")
                        ? this->params()->get<double>(CP_name+"/intracellular/antioxidant/capacity") : 1.0;
  const double k_scavenge = this->params()->have_parameter<double>(CP_name+"/intracellular/antioxidant/k_scavenge")
                          ? this->params()->get<double>(CP_name+"/intracellular/antioxidant/k_scavenge") : 0.0;
  const double alpha_h2o2 = this->params()->have_parameter<double>(CP_name+"/intracellular/alpha/H2O2")
                          ? this->params()->get<double>(CP_name+"/intracellular/alpha/H2O2") : 1.0;
  const double alpha_no2  = this->params()->have_parameter<double>(CP_name+"/intracellular/alpha/NO2_")
                          ? this->params()->get<double>(CP_name+"/intracellular/alpha/NO2_")  : 0.0;
  const double k_induce   = this->params()->have_parameter<double>(CP_name+"/intracellular/damage/k_induction")
                          ? this->params()->get<double>(CP_name+"/intracellular/damage/k_induction") : 0.0;
  const double k_repair   = this->params()->have_parameter<double>(CP_name+"/intracellular/damage/k_repair")
                          ? this->params()->get<double>(CP_name+"/intracellular/damage/k_repair") : 0.0;
  // uptake from extracellular fields (reduce grid accordingly)
  double uptake_h2o2 = 0.0, uptake_no2 = 0.0;
  // only query grids if the substances are present in the input file
  const auto& substances = this->params()->get<std::vector<std::string>>("substances");
  // Optional ECM barrier: dense ECM attenuates RONS penetration to cells.
  // Relevant for desmoplastic tumours (PDAC/BTC) where collagen-rich stroma
  // restricts diffusion and uptake of reactive species. Parameterised as
  // exponential attenuation: penetration = exp(-k_barrier * ECM_density).
  // Only active when intracellular/ecm_barrier/k_barrier is defined.
  double ecm_penetration_factor = 1.0;
  if (this->params()->have_parameter<double>(CP_name+"/intracellular/ecm_barrier/k_barrier"))
    {
      const double k_bar = this->params()->get<double>(CP_name+"/intracellular/ecm_barrier/k_barrier");
      if (k_bar > 0.0 &&
          std::find(substances.begin(), substances.end(), "ECM") != substances.end())
        {
          if (auto* dg_ecm = rm->GetDiffusionGrid("ECM"))
            {
              const double ecm = GetInterpolatedValue(dg_ecm, xyz, this->params());
              ecm_penetration_factor = std::exp(-k_bar * std::max(0.0, ecm));
              if (ecm_penetration_factor < 0.0) ecm_penetration_factor = 0.0;
              if (ecm_penetration_factor > 1.0) ecm_penetration_factor = 1.0;
            }
        }
    }
  if (std::find(substances.begin(), substances.end(), "H2O2") != substances.end()) {
    if (auto* dg = rm->GetDiffusionGrid("H2O2")) {
      const double c = GetInterpolatedValue(dg, xyz, this->params());
      const double desired = dt * std::max(0.0, k_uptake_H2O2 * c * ecm_penetration_factor);
      const double removed = std::min(std::max(0.0, c), desired);
      if (removed > 0.0) dg->ChangeConcentrationBy(xyz, -removed);
      uptake_h2o2 = removed / dt;
    }
  }
  if (std::find(substances.begin(), substances.end(), "NO2_") != substances.end()) {
    if (auto* dg = rm->GetDiffusionGrid("NO2_")) {
      const double c = GetInterpolatedValue(dg, xyz, this->params());
      const double desired = dt * std::max(0.0, k_uptake_NO2 * c * ecm_penetration_factor);
      const double removed = std::min(std::max(0.0, c), desired);
      if (removed > 0.0) dg->ChangeConcentrationBy(xyz, -removed);
      uptake_no2 = removed / dt;
    }
  }
  // intracellular ROS balance (semi-implicit: explicit production, exponential decay)
  const double prod = alpha_h2o2*uptake_h2o2 + alpha_no2*uptake_no2;
  const double decay_rate = k_scavenge * antioxidant_capacity_; // per hour
  if (decay_rate > 1.0e-12)
    {
      // Exponential integrator: exact solution for linear decay with constant source
      // [ROS](t+dt) = prod/decay + ([ROS](t) - prod/decay) * exp(-decay*dt)
      const double steady_state = prod / decay_rate;
      ros_internal_ = steady_state + (ros_internal_ - steady_state) * std::exp(-decay_rate * dt);
    }
  else
    {
      // No scavenging: pure accumulation
      ros_internal_ += dt * prod;
    }
  if (ros_internal_ < 0.0) ros_internal_ = 0.0;
  // DNA damage accumulation with repair (implicit Euler for stability)
  // D_new = D_old + dt*(k_ind*ROS - k_rep*D_new)
  // Rearranged: D_new = (D_old + dt*k_ind*ROS) / (1 + dt*k_rep)
  dna_damage_ = (dna_damage_ + dt * k_induce * ros_internal_) / (1.0 + dt * k_repair);
  if (dna_damage_ < 0.0) dna_damage_ = 0.0;
  UpdateDdrPathway(this);
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::RunECMInteraction()
{
  // Only viable cells interact with ECM
  if (!this->GetPhenotype()) return false;
  //
  auto* rm  = bdm::Simulation::GetActive()->GetResourceManager();
  auto* rg  = bdm::Simulation::GetActive()->GetRandom();
  const double dt = this->params()->get<double>("time_step");
  //
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const int mo = this->params()->get<int>(CP_name+"/mechanism_order");
  const bdm::Double3 xyz = this->GetPosition();
  //
  // Check if ECM field is present in this simulation
  const auto& substances = this->params()->get<std::vector<std::string>>("substances");
  const bool has_ecm = (std::find(substances.begin(), substances.end(), "ECM") != substances.end());
  if (!has_ecm) return false;
  //
  auto* dg_ecm = rm->GetDiffusionGrid("ECM");
  if (!dg_ecm) return false;
  //
  const double ecm_conc = GetInterpolatedValue(dg_ecm, xyz, this->params());
  //
  // --- ECM adhesion signal (integrin-mediated) ---
  const double integrin_sensitivity =
    this->params()->have_parameter<double>(CP_name+"/ecm/integrin_sensitivity")
    ? this->params()->get<double>(CP_name+"/ecm/integrin_sensitivity") : 1.0;
  const double adhesion_signal = ecm_conc * integrin_sensitivity;
  //
  // --- Anoikis check: low adhesion => stochastic apoptosis ---
  if (this->GetCanApoptose() &&
      this->params()->have_parameter<double>(CP_name+"/ecm/anoikis_threshold"))
    {
      const double anoikis_thr = this->params()->get<double>(CP_name+"/ecm/anoikis_threshold");
      if (adhesion_signal < anoikis_thr)
        {
          if (11 == mo || 12 == mo)
            {
              // Strict hazard form (Mechanism 11 and 12):
              //   P(anoikis) = 1 - exp(-k_anoikis * dt)
              // Preferred parameter: ecm/anoikis_hazard_rate [1/time]
              // Backward compatibility: if legacy per-step probability is set,
              // convert it to an equivalent hazard rate.
              double k_anoikis = 0.0;
              if (this->params()->have_parameter<double>(CP_name+"/ecm/anoikis_hazard_rate"))
                {
                  k_anoikis = std::max(0.0, this->params()->get<double>(CP_name+"/ecm/anoikis_hazard_rate"));
                }
              else
                {
                  const double p_legacy_raw =
                    this->params()->have_parameter<double>(CP_name+"/ecm/anoikis_probability")
                    ? this->params()->get<double>(CP_name+"/ecm/anoikis_probability") : 0.01;
                  const double p_legacy = std::max(0.0, std::min(1.0, p_legacy_raw));
                  if (p_legacy >= 1.0)
                    return true;
                  if (dt > 0.0 && p_legacy > 0.0)
                    k_anoikis = -std::log1p(-p_legacy) / dt;
                }
              // Optional resistance factor in [0,1]: 1 => complete anoikis resistance.
              if (this->params()->have_parameter<double>(CP_name+"/ecm/anoikis_resistance"))
                {
                  const double r_raw = this->params()->get<double>(CP_name+"/ecm/anoikis_resistance");
                  const double r = std::max(0.0, std::min(1.0, r_raw));
                  k_anoikis *= (1.0 - r);
                }
              const double P_anoikis = 1.0 - std::exp(-std::max(0.0, k_anoikis) * dt);
              if (rg->Uniform(0.0, 1.0) < P_anoikis)
                return true;
            }
          else
            {
              const double anoikis_prob =
                this->params()->have_parameter<double>(CP_name+"/ecm/anoikis_probability")
                ? this->params()->get<double>(CP_name+"/ecm/anoikis_probability") : 0.01;
              // Legacy behavior for non-Mechanism-11 cells
              if (rg->Uniform(0.0, 1.0) < anoikis_prob * dt)
                return true; // anoikis: caller should trigger Ap phase
            }
        }
    }
  //
  // --- ECM degradation by cell (MMP/protease-like activity) ---
  // Supports two parameter names for backward compatibility:
  //   ecm_degradation_rate  (preferred, Mechanism 11 naming)
  //   k_degrade             (legacy)
  // mmp_activity is an optional multiplier (reflects MMP expression level;
  //   enzymes: MMP-2, MMP-9, MT1-MMP, cathepsins, uPA).
  // Focal-adhesion dynamics (FAK/Src/talin) regulate MMP secretion.
  {
    double k_deg = 0.0;
    if (this->params()->have_parameter<double>(CP_name+"/ecm/ecm_degradation_rate"))
      k_deg = this->params()->get<double>(CP_name+"/ecm/ecm_degradation_rate");
    else if (this->params()->have_parameter<double>(CP_name+"/ecm/k_degrade"))
      k_deg = this->params()->get<double>(CP_name+"/ecm/k_degrade");
    // Optional MMP activity multiplier (mmp_activity = 0 disables degradation)
    if (this->params()->have_parameter<double>(CP_name+"/ecm/mmp_activity"))
      k_deg *= this->params()->get<double>(CP_name+"/ecm/mmp_activity");
    if (k_deg > 0.0 && ecm_conc > 0.0) {
      const double deg_amount = dt * k_deg * ecm_conc;
      const double actual_deg = std::min(ecm_conc, deg_amount);
      if (actual_deg > 0.0) dg_ecm->ChangeConcentrationBy(xyz, -actual_deg);
    }
  }
  //
  // --- ECM deposition by cell (matrix synthesis / CAF remodelling) ---
  // Supports parameter names:
  //   ecm_deposition_rate    (preferred, Mechanism 11 naming)
  //   caf_ecm_deposition_rate (for CAF-like cells; fibronectin, type-I collagen)
  //   k_deposit              (legacy)
  // All three are optional; if multiple are present they are summed.
  // Biologically: CAFs deposit collagen/fibronectin/periostin to stiffen the
  // tumour stroma; cancer cells may also secrete fibronectin for autocrine
  // adhesion. Uses logistic deposition: rate decreases near ECM saturation.
  {
    double k_dep = 0.0;
    if (this->params()->have_parameter<double>(CP_name+"/ecm/ecm_deposition_rate"))
      k_dep += this->params()->get<double>(CP_name+"/ecm/ecm_deposition_rate");
    if (this->params()->have_parameter<double>(CP_name+"/ecm/caf_ecm_deposition_rate"))
      k_dep += this->params()->get<double>(CP_name+"/ecm/caf_ecm_deposition_rate");
    if (k_dep <= 0.0 && this->params()->have_parameter<double>(CP_name+"/ecm/k_deposit"))
      k_dep = this->params()->get<double>(CP_name+"/ecm/k_deposit");
    if (k_dep > 0.0) {
      const double max_ecm =
        this->params()->have_parameter<double>(CP_name+"/ecm/ecm_saturation")
        ? this->params()->get<double>(CP_name+"/ecm/ecm_saturation") : 2.0;
      if (ecm_conc < max_ecm) {
        const double dep_amount = dt * k_dep * (1.0 - ecm_conc / max_ecm);
        dg_ecm->ChangeConcentrationBy(xyz, dep_amount);
      }
    }
  }
  //
  return false; // no anoikis
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::EvaluateG1SCheckpoint()
{
  // G1/S checkpoint: returns true if the G1->S transition should be BLOCKED.
  // Conditions checked: DNA damage, hypoxia, sparse ECM, local crowding.
  if (!this->GetPhenotype()) return false;
  //
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // --- DNA damage checkpoint (molecular DDR or legacy aggregate damage) ---
  if (bdm::IsMolecularG1SCheckpointBlocked(this))
    return true;
  if (this->params()->have_parameter<double>(CP_name+"/checkpoint/G1S/damage_threshold")
      && (!this->params()->have_parameter<bool>(CP_name+"/intracellular/ddr/enabled")
          || !this->params()->get<bool>(CP_name+"/intracellular/ddr/enabled")))
    {
      const double thr = this->params()->get<double>(CP_name+"/checkpoint/G1S/damage_threshold");
      if (dna_damage_ > thr) return true;
    }
  //
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  const auto& substances = this->params()->get<std::vector<std::string>>("substances");
  //
  // --- Oxygen/nutrient sufficiency ---
  if (std::find(substances.begin(), substances.end(), "O2") != substances.end())
    {
      if (auto* dg = rm->GetDiffusionGrid("O2"))
        {
          const double o2 = GetInterpolatedValue(dg, this->GetPosition(), this->params());
          if (this->params()->have_parameter<double>(CP_name+"/checkpoint/G1S/O2_threshold"))
            {
              const double o2_thr = this->params()->get<double>(CP_name+"/checkpoint/G1S/O2_threshold");
              if (o2 < o2_thr) return true; // block: hypoxic
            }
        }
    }
  // --- Nutrient sufficiency (mTORC1/eIF4E/AMPK gate for G1/S entry) ---
  // Low nutrient levels suppress CDK2 via AMPK → p21/p27 upregulation and
  // decreased cyclin D1 synthesis. Applies to any substance named
  // "nutrient", "Nutrient", "glucose", or "Glucose" if present.
  // Only active when checkpoint/G1S/nutrient_min is configured.
  {
    for (const char* nname : {"nutrient", "Nutrient", "glucose", "Glucose"}) {
      const std::string ns(nname);
      if (std::find(substances.begin(), substances.end(), ns) != substances.end()) {
        if (auto* dg = rm->GetDiffusionGrid(ns)) {
          const double nut = GetInterpolatedValue(dg, this->GetPosition(), this->params());
          if (this->params()->have_parameter<double>(CP_name+"/checkpoint/G1S/nutrient_min")) {
            const double nut_min = this->params()->get<double>(CP_name+"/checkpoint/G1S/nutrient_min");
            if (nut < nut_min) return true; // block: nutrient-depleted for S-phase
          }
        }
        break; // only the first matching nutrient substance is used
      }
    }
  }
  //
  // --- ECM density (sparse ECM => poor survival conditions for replication) ---
  if (std::find(substances.begin(), substances.end(), "ECM") != substances.end())
    {
      if (auto* dg = rm->GetDiffusionGrid("ECM"))
        {
          const double ecm = GetInterpolatedValue(dg, this->GetPosition(), this->params());
          if (this->params()->have_parameter<double>(CP_name+"/checkpoint/G1S/ECM_threshold"))
            {
              const double ecm_thr = this->params()->get<double>(CP_name+"/checkpoint/G1S/ECM_threshold");
              if (ecm < ecm_thr) return true; // block: ECM too sparse for cycling
            }
        }
    }
  //
  // --- Local crowding (contact inhibition of proliferation) ---
  if (this->params()->have_parameter<double>(CP_name+"/checkpoint/G1S/crowding_threshold"))
    {
      const double influence_ratio =
        this->params()->have_parameter<double>(CP_name+"/can_divide/influence_ratio")
        ? this->params()->get<double>(CP_name+"/can_divide/influence_ratio") : 2.0;
      const double occupancy = ComputeLocalOccupancyRatio(this->GetPosition(), influence_ratio);
      const double crowd_thr = this->params()->get<double>(CP_name+"/checkpoint/G1S/crowding_threshold");
      if (occupancy >= crowd_thr) return true; // block: contact inhibition
    }
  //
  return false; // all checks passed: allow G1->S
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::EvaluateG2MCheckpoint()
{
  // G2/M checkpoint: returns true if the G2->Di (mitosis) transition should be BLOCKED.
  // This is the critical gate preventing damaged cells from entering mitosis.
  if (!this->GetPhenotype()) return false;
  //
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // --- DNA damage checkpoint (ATM/CHK1/CHK2–Cdc25–CDK or legacy damage threshold) ---
  if (bdm::IsMolecularG2MCheckpointBlocked(this))
    return true;
  if (this->params()->have_parameter<double>(CP_name+"/checkpoint/G2M/damage_threshold")
      && (!this->params()->have_parameter<bool>(CP_name+"/intracellular/ddr/enabled")
          || !this->params()->get<bool>(CP_name+"/intracellular/ddr/enabled")))
    {
      const double thr = this->params()->get<double>(CP_name+"/checkpoint/G2M/damage_threshold");
      if (dna_damage_ > thr) return true;
    }
  //
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  const auto& substances = this->params()->get<std::vector<std::string>>("substances");
  //
  // --- Oxygen sufficiency (cells need energy for mitosis) ---
  if (std::find(substances.begin(), substances.end(), "O2") != substances.end())
    {
      if (auto* dg = rm->GetDiffusionGrid("O2"))
        {
          const double o2 = GetInterpolatedValue(dg, this->GetPosition(), this->params());
          if (this->params()->have_parameter<double>(CP_name+"/checkpoint/G2M/O2_threshold"))
            {
              const double o2_thr = this->params()->get<double>(CP_name+"/checkpoint/G2M/O2_threshold");
              if (o2 < o2_thr) return true; // block: insufficient oxygen for mitosis
            }
        }
    }
  //
  return false; // all checks passed: allow G2->Di
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckNecrosis()
{
  // Necrosis pathway — two independent triggers:
  //   Path 1 (original): severe hypoxia (O2 < threshold) + prolonged arrest.
  //   Path 2 (new):      extreme intracellular ROS exceeding membrane / energy
  //                      tolerance (relevant for high-dose CAP or metabolic stress).
  // On any trigger, the cell transforms to the necrotic phenotype (ID 0).
  if (!this->GetCanApoptose()) return false;
  if (!this->GetPhenotype())   return false;
  //
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  bool triggered = false;
  //
  // --- Path 2: extreme intracellular ROS / energy collapse ---
  if (!triggered &&
      this->params()->have_parameter<double>(CP_name+"/can_necrose/ros_threshold"))
    {
      const double ros_thr = this->params()->get<double>(CP_name+"/can_necrose/ros_threshold");
      if (ros_internal_ > ros_thr)
        {
          const double prob_ros =
            this->params()->have_parameter<double>(CP_name+"/can_necrose/ros_probability")
            ? this->params()->get<double>(CP_name+"/can_necrose/ros_probability") : 0.5;
          if (rg->Uniform(0.0, 1.0) <= prob_ros)
            triggered = true;
        }
    }
  //
  // --- Path 1: severe hypoxia + prolonged arrest (original behaviour) ---
  if (!triggered &&
      this->params()->have_parameter<double>(CP_name+"/can_necrose/O2_threshold"))
    {
      auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
      const auto& substances = this->params()->get<std::vector<std::string>>("substances");
      if (std::find(substances.begin(), substances.end(), "O2") != substances.end())
        {
          auto* dg = rm->GetDiffusionGrid("O2");
          if (dg)
            {
              const double o2 = GetInterpolatedValue(dg, this->GetPosition(), this->params());
              const double necr_thr = this->params()->get<double>(CP_name+"/can_necrose/O2_threshold");
              if (o2 <= necr_thr)
                {
                  const int min_arrest =
                    this->params()->have_parameter<int>(CP_name+"/can_necrose/min_arrest_time")
                    ? this->params()->get<int>(CP_name+"/can_necrose/min_arrest_time") : 1;
                  if (arrest_time_ >= min_arrest)
                    {
                      const double prob =
                        this->params()->have_parameter<double>(CP_name+"/can_necrose/probability")
                        ? this->params()->get<double>(CP_name+"/can_necrose/probability") : 0.5;
                      if (rg->Uniform(0.0, 1.0) <= prob)
                        triggered = true;
                    }
                }
            }
        }
    }
  //
  if (!triggered) return false;
  //
  // --- Perform necrotic transformation to phenotype 0 ---
  const int new_phenotype = 0;
  this->SetPhenotype(new_phenotype);
  this->SetAge();
  phase_age_ = 0;
  arrest_time_ = 0;
  is_quiescent_ = false;
  this->IncrementNumberOfTrasformations();
  //
  const std::string CP_new_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(new_phenotype));
  //
  // Reset behavior to match the new (necrotic) phenotype's mechanism order
  {
    const bdm::InlineVector<bdm::Behavior*,2>& behavior = this->GetAllBehaviors();
    if (behavior.size() == 1)
      this->RemoveBehavior(behavior[0]);
  }
  const int mo = this->params()->get<int>(CP_new_name+"/mechanism_order");
  if      (10==mo) this->AddBehavior(new Biology4BiologicalCell_10());
  else if (11==mo) this->AddBehavior(new Biology4BiologicalCell_11());
  else if (12==mo) this->AddBehavior(new Biology4BiologicalCell_12());
  else             ABORT_("an exception is caught");
  //
  return true; // transformation completed
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckApoptosisByDamage()
{
  if (!this->GetCanApoptose()) return false;
  if (!this->GetPhenotype()) return false;
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  const double dt = this->params()->get<double>("time_step");
  // phenotype-specific namespace
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const int mo = this->params()->get<int>(CP_name+"/mechanism_order");
  const double thr = this->params()->have_parameter<double>(CP_name+"/intracellular/damage/threshold")
                   ? this->params()->get<double>(CP_name+"/intracellular/damage/threshold") : 1.0e+99;
  const bool ddr_enabled =
    this->params()->have_parameter<bool>(CP_name+"/intracellular/ddr/enabled")
    && this->params()->get<bool>(CP_name+"/intracellular/ddr/enabled");
  if (ddr_enabled
      && this->params()->have_parameter<double>(CP_name+"/intracellular/ddr/apoptosis/p53_threshold"))
    {
      if (p53_active_ <= this->params()->get<double>(CP_name+"/intracellular/ddr/apoptosis/p53_threshold"))
        return false;
    }
  else if (dna_damage_ <= thr)
    return false;
  if (11 == mo || 12 == mo)
    {
      // Strict hazard form for Mechanism 11 and 12:
      //   P(apoptosis) = 1 - exp(-k_damage * dt)
      // Preferred parameter: intracellular/damage/hazard_rate [1/time]
      // Legacy compatibility: convert intracellular/damage/probability to hazard.
      double k_apopt = 0.0;
      if (this->params()->have_parameter<double>(CP_name+"/intracellular/damage/hazard_rate"))
        {
          k_apopt = std::max(0.0,
                             this->params()->get<double>(CP_name+"/intracellular/damage/hazard_rate"));
        }
      else if (this->params()->have_parameter<double>(CP_name+"/intracellular/damage/probability"))
        {
          const double p_raw = this->params()->get<double>(CP_name+"/intracellular/damage/probability");
          const double p = std::max(0.0, std::min(1.0, p_raw));
          if (p >= 1.0) return true;
          if (dt > 0.0 && p > 0.0) k_apopt = -std::log1p(-p) / dt;
        }
      else
        {
          // Legacy deterministic behavior when no probability parameter exists.
          return true;
        }
      const double P_apopt = 1.0 - std::exp(-std::max(0.0, k_apopt) * dt);
      return rg->Uniform(0.0,1.0) <= P_apopt;
    }
  // Legacy behavior for non-Mechanism-11 cells.
  if (this->params()->have_parameter<double>(CP_name+"/intracellular/damage/probability"))
    {
      const double p = this->params()->get<double>(CP_name+"/intracellular/damage/probability");
      if (rg->Uniform(0.0,1.0) > p) return false;
    }
  return true;
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckPositionValidity()
{
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0,
               deltaCOORD = maxCOORD - meanCOORD;
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  // check if simulation domain is bounded or unbounded (for cell outflux)
  if ( this->params()->get<bool>("simulation_domain_is_bounded") )
    {
      bdm::Double3 xyz = this->GetPosition();
      // identify mode of simulation domain
      if      ( this->params()->get<bool>("simulation_domain_is_polar") &&
                this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          const double radius = sqrt(pow2(xyz[0]-meanCOORD)
                                    +pow2(xyz[1]-meanCOORD));
          const double phi = atan2(xyz[1], xyz[0]);
          bool updated_xyz = false;
          if ( radius > deltaCOORD-tol )
            {
              xyz[0] = meanCOORD + (deltaCOORD-tol)*cos(phi);
              xyz[1] = meanCOORD + (deltaCOORD-tol)*sin(phi);
              updated_xyz = true;
            }
          if        ( xyz[2] < meanCOORD-tol )
            {
              xyz[2] = meanCOORD - tol;
              updated_xyz = true;
            }
          else if ( xyz[2] > meanCOORD+tol )
            {
              xyz[2] = meanCOORD + tol;
              updated_xyz = true;
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else if ( this->params()->get<bool>("simulation_domain_is_polar") &&
              ! this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          const double radius = sqrt(pow2(xyz[0]-meanCOORD)
                                    +pow2(xyz[1]-meanCOORD)
                                    +pow2(xyz[2]-meanCOORD));
          const double phi = atan2(xyz[1], xyz[0]);
          const double theta = atan2(sqrt(pow2(xyz[0])+pow2(xyz[1])), xyz[2]);
          bool updated_xyz = false;
          if ( radius > deltaCOORD-tol )
            {
              xyz[0] = meanCOORD + (deltaCOORD-tol)*sin(theta)*cos(phi);
              xyz[1] = meanCOORD + (deltaCOORD-tol)*sin(theta)*sin(phi);
              xyz[2] = meanCOORD + (deltaCOORD-tol)*cos(theta);
              updated_xyz = true;
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else if ( ! this->params()->get<bool>("simulation_domain_is_polar") &&
                  this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          bool updated_xyz = false;
          for (size_t d=0; d<2; d++)
            {
              if      ( xyz[d] < minCOORD+tol )
                {
                  xyz[d] = minCOORD + tol;
                  updated_xyz = true;
                }
              else if ( xyz[d] > maxCOORD-tol )
                {
                  xyz[d] = maxCOORD - tol;
                  updated_xyz = true;
                }
            }
          if      ( xyz[2] < meanCOORD-tol )
            {
              xyz[2] = meanCOORD - tol;
              updated_xyz = true;
            }
          else if ( xyz[2] > meanCOORD+tol )
            {
              xyz[2] = meanCOORD + tol;
              updated_xyz = true;
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else if ( ! this->params()->get<bool>("simulation_domain_is_polar") &&
                ! this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          bool updated_xyz = false;
          for (size_t d=0; d<3; d++)
            {
              if      ( xyz[d] < minCOORD+tol )
                {
                  xyz[d] = minCOORD + tol;
                  updated_xyz = true;
                }
              else if ( xyz[d] > maxCOORD-tol )
                {
                  xyz[d] = maxCOORD - tol;
                  updated_xyz = true;
                }
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else
        ABORT_("an exception is caught");
    }
  else
    {
      bool reached_boundary = false;
      const bdm::Double3 xyz = this->GetPosition();
      // spherical coordinate: radius
      double radius = 0.0;
      // identify mode of simulation domain
      if (this->params()->get<bool>("simulation_domain_is_2D"))
      {
        radius = sqrt(pow2(xyz[0]-meanCOORD)
                     +pow2(xyz[1]-meanCOORD));
        if (this->params()->get<bool>("simulation_domain_is_polar"))
        {
          if ( radius > deltaCOORD-tol )
            reached_boundary = true;
        }
        else
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol )
            reached_boundary = true;
        }
      }
      else
      {
        radius = sqrt(pow2(xyz[0]-meanCOORD)
                     +pow2(xyz[1]-meanCOORD)
                     +pow2(xyz[2]-meanCOORD));
        if (this->params()->get<bool>("simulation_domain_is_polar"))
        {
          if ( radius > deltaCOORD-tol )
            reached_boundary = true;
        }
        else
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol ||
               xyz[2] < minCOORD+tol || xyz[2] > maxCOORD-tol )
            reached_boundary = true;
        }
      }
      //
      if (reached_boundary)
      {
        // spherical coordinates: inclination (theta), azimuth (phi)
        double theta = 0.0, phi = 0.0;
        // identify mode of simulation domain
        if (this->params()->get<bool>("simulation_domain_is_2D"))
          theta = 0.5*bdm::Math::kPi;
        else
          theta = atan2(sqrt(pow2(xyz[0]-meanCOORD)+pow2(xyz[1]-meanCOORD)),
                        (xyz[2]-meanCOORD));
        phi = atan2((xyz[1]-meanCOORD), (xyz[0]-meanCOORD));
        // reset the spherical coordinates for all escaping cells:
        // radius, inclination (theta) and azimuth (phi)
        std::vector<double>& escaped_cells_radius =
          this->params()->set<std::vector<double>>(CP_name+"/escaped_cells/radius");
        std::vector<double>& escaped_cells_theta  =
          this->params()->set<std::vector<double>>(CP_name+"/escaped_cells/theta");
        std::vector<double>& escaped_cells_phi    =
          this->params()->set<std::vector<double>>(CP_name+"/escaped_cells/phi");
        //
        std::vector<int>& escaped_cells_phase =
          this->params()->set<std::vector<int>>(CP_name+"/escaped_cells/phase");
        //
        escaped_cells_radius.push_back(radius);
        escaped_cells_theta.push_back(theta);
        escaped_cells_phi.push_back(phi);
        //
        escaped_cells_phase.push_back(this->GetPhase());
        // cell has not remained inside the simulation domain
        return false;
      }
      // in case of a 2D simulation domain, ensure you
      // "correct" and bring back any cells on-plane
      if (this->params()->get<bool>("simulation_domain_is_2D"))
      {
        bdm::Double3 xyz = this->GetPosition();
        if      ( xyz[2] < meanCOORD-tol )
        {
          xyz[2] = meanCOORD - tol;
          this->SetPosition(xyz);
        }
        else if ( xyz[2] > meanCOORD+tol )
        {
          xyz[2] = meanCOORD + tol;
          this->SetPosition(xyz);
        }
      }
    }
  // access the pointer to parameter of the simulation obstacles object
  const SimulationObstacles* obstacles =
    this->params()->get<SimulationObstacles*>("simulation_obstacles");
  // check if cell has reached any of the 'box' simulation obstacles
  for (size_t l=0; l<obstacles->box.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      // it true, then cell position needs to be corrected
      if ("box/inside"==obstacles->box[l].type)
        {
          // check if cell position is outside this obstacle
          if (! obstacles->box[l].is_inside(xyz))
            continue;
          // direction vector from point to center
          bdm::Double3 d_vector = xyz - obstacles->box[l].center;
          if (!normalize(d_vector, d_vector))
            ABORT_("could not normalize the direction vector");
          //
          const double dot_0 = (d_vector*obstacles->box[l].laxis_0),
                       dot_1 = (d_vector*obstacles->box[l].laxis_1),
                       dot_2 = (d_vector*obstacles->box[l].laxis_2);
          //
          std::map<double, int, std::greater<double>> sorted_map;
          sorted_map.insert( std::make_pair(fabs(dot_0), 0) );
          sorted_map.insert( std::make_pair(fabs(dot_1), 1) );
          sorted_map.insert( std::make_pair(fabs(dot_2), 2) );
          const int dir = sorted_map.begin()->second;
          // origin (center) point to box face
          bdm::Double3 origin;
          // outward unit normal vector to box face
          bdm::Double3 normal;
          if (0==dir)
            {
              if (dot_0>0.0)
                { origin = obstacles->box[l].center2face_0; normal = obstacles->box[l].laxis_0*(+1.0); }
              else
                { origin = obstacles->box[l].center2face_1; normal = obstacles->box[l].laxis_0*(-1.0); }
              normal = obstacles->box[l].laxis_0;
            }
          else if (1==dir)
            {
              if (dot_1>0.0)
                { origin = obstacles->box[l].center2face_2; normal = obstacles->box[l].laxis_1*(+1.0); }
              else
                { origin = obstacles->box[l].center2face_3; normal = obstacles->box[l].laxis_1*(-1.0); }
              normal = obstacles->box[l].laxis_1;
            }
          else
            {
              if (dot_2>0.0)
                { origin = obstacles->box[l].center2face_4; normal = obstacles->box[l].laxis_2*(+1.0); }
              else
                { origin = obstacles->box[l].center2face_5; normal = obstacles->box[l].laxis_2*(-1.0); }
            }
          //
          xyz = project_to_plane(normal, origin, xyz)
              + normal * (rg->Uniform(0.5,0.8)*this->GetDiameter());
          // enforce cell to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
        }
      else if ("box/outside"==obstacles->box[l].type)
        {
          // check if cell position is inside this obstacle
          if (obstacles->box[l].is_inside(xyz))
            continue;
          // calculate the outward unit normal vector and the projection vector of
          // the cell to each box face
          std::vector< std::pair<bdm::Double3,bdm::Double3> > data_vector;
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_0;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_0, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_1;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_1, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_2;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_2, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_3;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_3, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_4;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_4, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_5;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_5, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          //
          for (auto d : data_vector)
            {
              const bdm::Double3& normal = d.first;
              const bdm::Double3& proj   = d.second;
              //
              const bdm::Double3 xyz_proj = xyz - proj;
              if (xyz_proj*normal>0.0) continue;
              //
              xyz = proj
                  + normal * (rg->Uniform(0.5,0.8)*this->GetDiameter());
              // enforce cell to lie on the (box) obstacle surface
              this->SetPosition(xyz);
              //
              displace += xyz;
              this->UpdateTrail(L2norm(displace));
              // cell has remained inside the simulation domain
              return true;
            }
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'box' simulation obstacles loop
    }
  // check if cell has reached any of the 'sphere' simulation obstacles
  for (size_t l=0; l<obstacles->sphere.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      // it true, then cell position needs to be corrected
      if ("sphere/inside"==obstacles->sphere[l].type)
        {
          // check if cell position is outside this obstacle
          if (! obstacles->sphere[l].is_inside(xyz))
            continue;
          // direction vector from point to center
          bdm::Double3 d_vector = xyz - obstacles->sphere[l].center;
          if (!normalize(d_vector, d_vector))
            ABORT_("could not normalize the direction vector");
          //
          const double delta = obstacles->sphere[l].radius
                             + (rg->Uniform(0.4,0.8)*this->GetDiameter());
          xyz = obstacles->sphere[l].center + d_vector * delta;
          // enforce cell to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
        }
      else if ("sphere/outside"==obstacles->sphere[l].type)
        {
          // check if cell position is inside this obstacle
          if (obstacles->sphere[l].is_inside(xyz))
            continue;
          // direction vector from point to center
          bdm::Double3 d_vector = xyz - obstacles->sphere[l].center;
          if (!normalize(d_vector, d_vector))
            ABORT_("could not normalize the direction vector");
          //
          const double delta = obstacles->sphere[l].radius
                             - (rg->Uniform(0.4,0.8)*this->GetDiameter());
          xyz = obstacles->sphere[l].center + d_vector * delta;
          // enforce cell to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'sphere' simulation obstacles loop
    }
  // check if cell has reached any of the 'surface' simulation obstacles
  for (size_t l=0; l<obstacles->surface.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      //
      const double safe_distance = 0.75 * this->GetDiameter();
      //
      std::map<double, std::pair<unsigned int,bdm::Double3>> tri3_proj;
      //
      const unsigned int n_tri3 = obstacles->surface[l].triangle.size();
      for (unsigned int t=0; t<n_tri3; t++)
        {
          const ObstacleSTL::Triangle& tri3 = obstacles->surface[l].triangle[t];
          // origin (center) point to triangle
          const bdm::Double3& origin = tri3.center;
          // outward unit normal vector to triangle
          const bdm::Double3& normal = tri3.normal;
          // internal point wrt the user-defined surface
          const bdm::Double3& l0 = tri3.inside;
          // cell position intersection to the user-defined surface
          bdm::Double3 intx;
          if (! line_intersects_plane(l0, xyz, normal, origin, intx))
            continue;
          // skip following computations if projection point is outside triangle
          if (! is_inside_triangle(tri3.vertex_0, tri3.vertex_1, tri3.vertex_2, intx))
            continue;
          // cell position projection to the user-defined surface
          const bdm::Double3 proj = project_to_plane(normal, origin, xyz);
          //
          const bdm::Double3 xyz_proj = xyz - proj;
          const double distance = L2norm(xyz_proj);
          // skip following if cell is well outside the user-defined surface
          if ((xyz_proj*normal)>0.0 && distance>this->GetDiameter())
            continue;
          //
          auto data = std::make_pair(t, proj);
          tri3_proj.insert( std::make_pair(distance, data) );
          // ...end of triangles loop
        }
      //
      if (!tri3_proj.empty())
        {
          std::map<double, std::pair<unsigned int,bdm::Double3>>::const_iterator
            ci = tri3_proj.begin();
          // access the triangle first...
          const unsigned int t = ci->second.first;
          const ObstacleSTL::Triangle& tri3 = obstacles->surface[l].triangle[t];
          // outward unit normal vector to triangle
          bdm::Double3 normal = tri3.normal;
          if (!normalize(normal, normal))
            ABORT_("could not normalize the normal vector");
          // cell position projection to the user-defined surface
          const bdm::Double3& proj = ci->second.second;
          //
          xyz = proj + normal * safe_distance;
          // enforce cell to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
          // ...end of if-case
        }
      // ...end of 'surface' simulation obstacles loop
    }
  // check if cell has reached any of the 'scaffold' simulation obstacles
  for (size_t l=0; l<obstacles->scaffold.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      //
      const unsigned int n_segm = obstacles->scaffold[l].segment.size();
      for (unsigned int s=0; s<n_segm; s++)
        {
          const ObstacleScaffold::Segment& segm = obstacles->scaffold[l].segment[s];
          //
          const bdm::Double3 n0 = segm.vertex_0,
                             n1 = segm.vertex_1;
          // check if cell position is inside this obstacle
          if (! is_inside_segment(n0, n1, xyz))
            continue;
          // cell position projection to the user-defined segment
          const bdm::Double3 proj = project_to_line(n0, n1, xyz);
          const double distance = L2norm(bdm::Double3(xyz-proj))
                                - segm.radius;
          //
          if (distance>rg->Uniform(0.5,1.0)*this->GetDiameter()) continue;
          //
          bdm::Double3 normal = xyz - proj;
          if (!normalize(normal, normal))
            ABORT_("could not normalize the normal vector");
          //
          const double delta = (rg->Uniform(0.5,1.0)*this->GetDiameter());
          xyz = proj + normal * delta;
          // enforce cell to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
          // ...end of segments loop
        }
      // ...end of 'scaffold' simulation obstacles loop
    }
  // cell has remained inside the simulation domain
  return true;
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckApoptosisAging()
{
  if (!this->GetCanApoptose()) return false;
  //
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  const double dt = this->params()->get<double>("time_step");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const int mo = this->params()->get<int>(CP_name+"/mechanism_order");
  const int cell_maturity = this->params()->get<int>(CP_name+"/can_apoptose/time_window");
  //
  // Age-gated apoptosis only after maturity threshold.
  if (this->GetAge() <= cell_maturity) return false;

  if (11 == mo)
    {
      // Strict hazard form for Mechanism 11:
      //   P(apoptosis_age) = 1 - exp(-k_age * dt)
      // Preferred parameters:
      //   can_apoptose/hazard_rate
      //   can_apoptose/hazard_rate_increment_with_age
      // Legacy compatibility:
      //   can_apoptose/probability and probability_increment_with_age
      double k_age = 0.0;
      const bool has_hazard =
        this->params()->have_parameter<double>(CP_name+"/can_apoptose/hazard_rate") ||
        this->params()->have_parameter<double>(CP_name+"/can_apoptose/hazard_rate_increment_with_age");

      if (has_hazard)
        {
          const double k0 = this->params()->have_parameter<double>(CP_name+"/can_apoptose/hazard_rate")
                          ? this->params()->get<double>(CP_name+"/can_apoptose/hazard_rate") : 0.0;
          const double kinc = this->params()->have_parameter<double>(CP_name+"/can_apoptose/hazard_rate_increment_with_age")
                            ? this->params()->get<double>(CP_name+"/can_apoptose/hazard_rate_increment_with_age") : 0.0;
          k_age = std::max(0.0, k0 + kinc * this->GetAge());
        }
      else
        {
          const double p0 = this->params()->have_parameter<double>(CP_name+"/can_apoptose/probability")
                          ? this->params()->get<double>(CP_name+"/can_apoptose/probability") : 0.0;
          const double pinc = this->params()->have_parameter<double>(CP_name+"/can_apoptose/probability_increment_with_age")
                            ? this->params()->get<double>(CP_name+"/can_apoptose/probability_increment_with_age") : 0.0;
          const double p_raw = p0 + pinc * this->GetAge();
          const double p = std::max(0.0, std::min(1.0, p_raw));
          if (p >= 1.0) return true;
          if (dt > 0.0 && p > 0.0) k_age = -std::log1p(-p) / dt;
        }

      const double P_apopt_age = 1.0 - std::exp(-std::max(0.0, k_age) * dt);
      return rg->Uniform(0.0,1.0) <= P_apopt_age;
    }

  // Legacy behavior for non-Mechanism-11 cells.
  if (this->params()->get<double>(CP_name+"/can_apoptose/probability_increment_with_age")>0.0)
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_apoptose/probability")
                                +this->params()->get<double>(CP_name+"/can_apoptose/probability_increment_with_age")
                                *this->GetAge() )
        return false;
    }
  else
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_apoptose/probability"))
        return false;
    }
  return true;
  //...end of cell apoptosis
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckApoptosis()
{
  if (!this->GetCanApoptose()) return false;
  // by design only viable (non-necrotic) cells could apoptose
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const double dt = this->params()->get<double>("time_step");
  const int mo = this->params()->get<int>(CP_name+"/mechanism_order");
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string& BC_name = dg->GetContinuumName(); // biochemical name
      //
      if (! this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/threshold"))
        continue;
      //
      const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                   threshold = this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          // allow cell apoptosis controlled by a combination of two biochemical cues, therefore
          // cell survival is dependent from another substance as well
          if (this->params()->have_parameter<bool>(CP_name+"/can_apoptose/"+BC_name+"/dependency") &&
              this->params()->get<bool>(CP_name+"/can_apoptose/"+BC_name+"/dependency"))
            {
              // iterate for all OTHER substances
              for ( std::vector<std::string>::const_iterator
                    cj=substances.begin(); cj!=substances.end(); cj++ )
                {
                  if ( cj == ci ) continue;
                  // access the BioDynaMo diffusion grid
                  auto* dg_other = rm->GetDiffusionGrid(*cj);
                  const std::string& BC_other_name = dg_other->GetContinuumName(); // biochemical name
                  //
                  if (! this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/threshold"))
                    continue;
                  //
                  const double concentration_other = GetInterpolatedValue(dg_other, this->GetPosition(), this->params()),
                               threshold_other = this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/threshold");
                  //
                  if ( ( threshold_other > 0.0 && concentration_other > +threshold_other ) ||
                       ( threshold_other < 0.0 && concentration_other < -threshold_other ) )
                    {
                      if (11 == mo)
                        {
                          // Strict hazard form for Mechanism 11:
                          // dependency-specific apoptosis
                          //   P = 1 - exp(-k * dt)
                          double k_apopt = 0.0;
                          if (this->params()->have_parameter<double>(
                                CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/hazard_rate"))
                            {
                              k_apopt = std::max(0.0, this->params()->get<double>(
                                CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/hazard_rate"));
                            }
                          else if (this->params()->have_parameter<double>(
                                     CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                            {
                              const double p_raw = this->params()->get<double>(
                                CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/probability");
                              const double p = std::max(0.0, std::min(1.0, p_raw));
                              if (p >= 1.0) return true;
                              if (dt > 0.0 && p > 0.0) k_apopt = -std::log1p(-p) / dt;
                            }
                          else
                            {
                              // Legacy deterministic behavior if no probability/hazard exists.
                              return true;
                            }
                          const double P_apopt = 1.0 - std::exp(-std::max(0.0, k_apopt) * dt);
                          if (rg->Uniform(0.0,1.0) <= P_apopt)
                            return true;
                        }
                      else
                        {
                          if (! this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                            return true;
                          else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                            return true;
                        }
                    }
                  //...end of other substances loop
                }
            }
          if (11 == mo)
            {
              // Strict hazard form for Mechanism 11:
              // direct apoptosis from BC threshold crossing
              //   P = 1 - exp(-k * dt)
              double k_apopt = 0.0;
              if (this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/hazard_rate"))
                {
                  k_apopt = std::max(0.0,
                    this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/hazard_rate"));
                }
              else if (this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/probability"))
                {
                  const double p_raw = this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/probability");
                  const double p = std::max(0.0, std::min(1.0, p_raw));
                  if (p >= 1.0) return true;
                  if (dt > 0.0 && p > 0.0) k_apopt = -std::log1p(-p) / dt;
                }
              else
                {
                  // Legacy deterministic behavior if no probability/hazard exists.
                  return true;
                }
              const double P_apopt = 1.0 - std::exp(-std::max(0.0, k_apopt) * dt);
              if (rg->Uniform(0.0,1.0) <= P_apopt)
                return true;
            }
          // Legacy behavior for non-Mechanism-11 cells.
          else if (! this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/probability"))
            {
              // since cell has apoptosed, then it must be removed from simulation
              return true;
            }
          else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/probability"))
            {
              // since cell has apoptosed, then it must be removed from simulation
              return true;
            }
        }
      //...end of substances loop
    }
  // since cell has not been through apoptosis, then it can do other things
  return false;
  //...end of cell apoptosis
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckAfterApoptosis()
{
  if (!this->GetCanApoptose()) return false;
  // by design only viable (non-necrotic) cells could apoptose
  if (!this->GetPhenotype()) return false;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const int cell_time_window = this->params()->get<int>(CP_name+"/can_apoptose/time_window/to_delete");
  //
  // since cell has apoptosed (due to ageing), then it must be removed from simulation
  if (this->GetAge()>cell_time_window) return true;
  //
  return false;
  //...end of cell apoptosis inquiry
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckQuiescenceAfterDivision()
{
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could go through this phase
  if (!this->GetPhenotype()) return false;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const int n_div = this->GetNumberOfDivisions();
  if (n_div==0) return false;
  //
  if (! this->params()->have_parameter<int>(CP_name+"/can_divide/quiescence/time_window"))
    return false;
  const int cell_quiescence = this->params()->get<int>(CP_name+"/can_divide/quiescence/time_window");
  //
  if (this->GetAge()<cell_quiescence) return true;
  //
  return false;
  //...end of cell quiescence after division
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckMigration()
{
  if (!this->GetCanMigrate()) return false;
  // by design only viable (non-necrotic) cells could migrate
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  bool has_migrated = false;
  //
  this->passive_displacement_ = {0.0, 0.0, 0.0};
  //
  // check if cell migrates passively due to some convective field
  if ( this->params()->have_parameter<std::string>("convection/dynamic/from_file") )
    {
      const int SpaceDimension = this->params()->get<bool>("simulation_domain_is_2D")
                               ? 2 : 3;
      const double time_step = this->params()->get<double>("time_step");
      const double max_adhesion =
        this->params()->get<double>(CP_name+"/can_migrate/max_adhesion/displacement");
      //
      bdm::Double3 dvec = {0.0, 0.0, 0.0};
      // ensure cell is well within the simulation domain!
      if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
      // iterate for all components of the convection (vector) field
      for (int ispdm=0; ispdm<SpaceDimension; ispdm++)
        {
          const std::string name = "convection_" + std::to_string(ispdm);
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(name);
          // obtain the convection component
          const double velocity_comp = GetInterpolatedValue(dg, this->GetPosition(), this->params());
          // calculate corresponding displacement component
          const double displacement_comp = velocity_comp * time_step;
          //
          dvec[ispdm] += displacement_comp;
        }
      //
      const double d_magn = L2norm(dvec);
      // check if distance covered is above a minimum, else ignore
      if (d_magn > this->params()->get<double>("migration_tolerance"))
      // check if convection contribution is significant compared to the cell adhesion property
      if (d_magn > max_adhesion)
        {
          // update the (cell) displacement vector
          this->passive_displacement_ += dvec;
          // scale the (cell) displacement vector wrt the adhesion effect
          this->passive_displacement_ *= ((d_magn-max_adhesion)/d_magn);
          // update this flag
          has_migrated = true;
        }
    }
  //
  // Active migration: timestep-consistent speed * dt, persistence, and chemotaxis.
  bdm::migration::RunActiveMigration(
    this, CP_name, this->active_displacement_, has_migrated);
  //
  // check if cell has migrated, if so then revise its spatial coordinates and trail
  if ( has_migrated )
    {
      this->UpdatePosition(this->GetDisplacement());
      this->UpdateTrail(L2norm(this->GetDisplacement()));
    }
  //
  // check if cell has migrated, if so then revise the spatial coordinates
  // of the cell protrusions; however, check for current algorithmic limitations!!!
  if ( has_migrated && this->GetNumberOfProtrusions() )
  // current implementation assumes that cell migration cannot be applied when
  // a cell has produced protrusions that have created sprouts or/and branches
  {
    if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
      ABORT_("an internal error occurred");

    // iterate for all (existing) protrusions of this cell
    for (int p=0; p<this->GetNumberOfProtrusions(); p++)
    {
      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
      // only a cell with filodia can be displaced (together with its migrating cell)
      if ( protrusion->IsTerminal() )
        {
          const bdm::Double3 xyz = protrusion->GetPosition();
          protrusion->SetPosition(xyz+this->GetDisplacement());
        }
      else
        ABORT_("an exception is caught");
    }
    // completed the cell protrusion displacement task
  }
  //
  return has_migrated;
  //...end of cell migration
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckTransformation()
{
  if (!this->GetCanTransform()) return false;
  // by design only viable (non-necrotic) cells could transform
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_transform/probability"))
    return false;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
    // iterate for all substances
    for ( std::vector<std::string>::const_iterator
          ci=substances.begin(); ci!=substances.end(); ci++ )
      {
        // access the BioDynaMo diffusion grid
        auto* dg = rm->GetDiffusionGrid(*ci);
        const std::string BC_name = dg->GetContinuumName(); // biochemical name
        //
        if (! this->params()->have_parameter<int>(CP_name+"/can_transform/"+BC_name+"/new_phenotype"))
          continue;
        //
        const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                     threshold = this->params()->get<double>(CP_name+"/can_transform/"+BC_name+"/threshold");
        //
        if ( ( threshold > 0.0 && concentration > +threshold ) ||
             ( threshold < 0.0 && concentration < -threshold ) )
          {
            if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_transform/"+BC_name+"/probability"))
              // now check if cell age is within an appropriate time-window to allow its transformation
              if ( this->GetAge() >= this->params()->get<int>(CP_name+"/can_transform/"+BC_name+"/time_window_open" ) &&
                   this->GetAge() <= this->params()->get<int>(CP_name+"/can_transform/"+BC_name+"/time_window_close") )
                  {
                    const int new_phenotype = this->params()->get<int>(CP_name+"/can_transform/"+BC_name+"/new_phenotype");
                    // firstly, the cell transforms
                    this->SetPhenotype(new_phenotype);
                    // increment this index
                    this->IncrementNumberOfTrasformations();
                    // now reset the age of the cell
                    this->SetAge();
                    //
                    const std::string CP_new_name = // cell phenotype name
                      this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
                    // principal directions of the cell polarization matrix
                    double p0, p1, p2;
                    if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                      {
                        p0 = this->params()->get<double>(CP_new_name+"/principal/0");
                        p1 = this->params()->get<double>(CP_new_name+"/principal/1");
                        p2 = this->params()->get<double>(CP_new_name+"/principal/2");
                      }
                    //
                    this->SetCanApoptose(this->params()->get<bool>(CP_new_name+"/can_apoptose"));
                    this->SetCanGrow(this->params()->get<bool>(CP_new_name+"/can_grow"));
                    this->SetCanDivide(this->params()->get<bool>(CP_new_name+"/can_divide"));
                    this->SetCanMigrate(this->params()->get<bool>(CP_new_name+"/can_migrate"));
                    this->SetCanTransform(this->params()->get<bool>(CP_new_name+"/can_transform"));
                    this->SetCanProtrude(this->params()->get<bool>(CP_new_name+"/can_protrude"));
                    this->SetCanPolarize(this->params()->get<bool>(CP_new_name+"/can_polarize"));
                    // update adherence and density for the new phenotype
                    this->SetAdherence(this->params()->have_parameter<double>(CP_new_name+"/adherence")
                                     ? this->params()->get<double>(CP_new_name+"/adherence") : 0.0);
                    if (this->params()->have_parameter<double>(CP_new_name+"/density"))
                      this->SetDensity(this->params()->get<double>(CP_new_name+"/density"));
                    // reset the cell polarization matrix
                    if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                      this->SetPolarization(diag(p0, p1, p2));
                    // reset the cell protrusion phenotype
                    if ( this->GetNumberOfProtrusions() )
                      {
                        if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
                          ABORT_("an internal error occurred");
                        //
                        // iterate for all (existing) protrusions of this cell
                        for (int p=0; p<this->GetNumberOfProtrusions(); p++)
                          {
                            auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
                            // assign this cell (that is associated with) to the protrusion created
                            protrusion->SetCell(this);
                          }
                      }
                    // reset the cell behavior (mechanisms order) from old to new one
                    {
                      const bdm::InlineVector<bdm::Behavior*,2>& behavior = this->GetAllBehaviors();
                      if (behavior.size()!=1)
                        ABORT_("an internal error occurred");
                      //
                      this->RemoveBehavior(behavior[0]);
                    }
                    const int mo = this->params()->get<int>(CP_new_name+"/mechanism_order");
                    if      (10==mo) this->AddBehavior(new Biology4BiologicalCell_10());
                    else if (11==mo) this->AddBehavior(new Biology4BiologicalCell_11());
                    else if (12==mo) this->AddBehavior(new Biology4BiologicalCell_12());
                    else             ABORT_("an exception is caught");
                    // cell has transformed, then proceed to check if it can do other things
                    return true;
                  }
              //
            //
          }
        //...end of substances loop
      }
  // cell has not been through any transformation
  return false;
  //...end of cell transformation
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckPolarization()
{
  if (!this->GetCanPolarize()) return false;
  // by design only viable (non-necrotic) cells could polarize
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // current implementation assumes that cell polarization cannot be applied when
  // a cell has produced protrusions
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_polarize/probability"))
    return false;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  const double principal_min = this->params()->get<double>(CP_name+"/principal/min"),
               principal_max = this->params()->get<double>(CP_name+"/principal/max");
  const std::vector<int> perm =
    this->params()->get<std::vector<int>>(CP_name+"/principal/permutation");
  const bdm::Double3 principal = { this->params()->get<double>(CP_name+"/principal/"+std::to_string(perm[0])) ,
                              this->params()->get<double>(CP_name+"/principal/"+std::to_string(perm[1])) ,
                              this->params()->get<double>(CP_name+"/principal/"+std::to_string(perm[2])) };
  //
  if ( this->params()->get<bool>(CP_name+"/can_polarize/migration/dependency") )
    {
      const int pattern = this->params()->have_parameter<int>(CP_name+"/can_polarize/migration/pattern")
                        ? this->params()->get<int>(CP_name+"/can_polarize/migration/pattern") : 0;
      // check if cell displacement is considerable to allow the cell self-polarization
      if (L2norm(this->GetDisplacement()) > this->params()->get<double>("migration_tolerance"))
        {
          bdm::Double3 v = this->GetDisplacement(); // cell displacement (vector) field
          // eigenvectors
          bdm::Double3 n0 = this->params()->get<bool>("simulation_domain_is_2D")
                     ? bdm::Double3{v[0], v[1], 0.0}
                     : bdm::Double3{v[0], v[1], v[2]};
          bdm::Double3 n1 = this->params()->get<bool>("simulation_domain_is_2D")
                     ? bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), 0.0}
                     : bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), rg->Uniform(-1.,1.)};
          bdm::Double3 n2 = cross(n0, n1);
          n1 = cross(n2, n0);
          // normalize all vectors
          n0.Normalize();
          n1.Normalize();
          n2.Normalize();
          // eigenvalues
          double p[3];
          if (0==pattern)
            {
              p[0] = rg->Uniform(principal[0],  principal_max),
              p[1] = rg->Uniform(principal[1],  principal[0]),
              p[2] = rg->Uniform(principal_min, principal[2]);
            }
          else if (1==pattern)
            {
              p[0] = rg->Uniform(principal[0],  principal_max),
              p[1] = rg->Uniform(principal[2],  principal[1]),
              p[2] = rg->Uniform(principal_min, principal[2]);
            }
          else
            ABORT_("an exception is caught");
          // tensor products of eigenvectors, scaled by respective eigenvalues
          const bdm::Double3x3 n0Xn0_p0 = tensor(n0, n0, p[0]),
                               n1Xn1_p1 = tensor(n1, n1, p[1]),
                               n2Xn2_p2 = tensor(n2, n2, p[2]);
          // update the cell polarization matrix
          this->polarize_ = n0Xn0_p0 + n1Xn1_p1 + n2Xn2_p2;
          // cell has polarized, then proceed to check if it can do other things
          return true;
        }
    }
  //
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
    // iterate for all substances if cell can
    // re-orient / polarize
    for ( std::vector<std::string>::const_iterator
          ci=substances.begin(); ci!=substances.end(); ci++ )
      {
        // access the BioDynaMo diffusion grid
        auto* dg = rm->GetDiffusionGrid(*ci);
        const std::string BC_name = dg->GetContinuumName(); // biochemical name
        const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                     threshold = this->params()->get<double>(CP_name+"/can_polarize/"+BC_name+"/threshold");
        //
        if ( ( threshold > 0.0 && concentration > +threshold ) ||
             ( threshold < 0.0 && concentration < -threshold ) )
          {
            bdm::Double3 v; // substance gradient (vector) field
            dg->GetGradient(this->GetPosition(), &v);
            // eigenvectors
            bdm::Double3 n0 = this->params()->get<bool>("simulation_domain_is_2D")
                       ? bdm::Double3{v[0], v[1], 0.0}
                       : bdm::Double3{v[0], v[1], v[2]};
            bdm::Double3 n1 = this->params()->get<bool>("simulation_domain_is_2D")
                       ? bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), 0.0}
                       : bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), rg->Uniform(-1.,1.)};
            bdm::Double3 n2 = cross(n0, n1);
            n1 = cross(n2, n0);
            // normalize all direction vectors
            n0.Normalize();
            n1.Normalize();
            n2.Normalize();
            // eigenvalues
            const double p0 = rg->Uniform(principal[0],  principal_max),
                         p1 = rg->Uniform(principal[1],  p0),
                         p2 = rg->Uniform(principal_min, principal[2]);
            // tensor products of eigenvectors, scaled by respective eigenvalues
            const bdm::Double3x3 n0Xn0_p0 = tensor(n0, n0, p0),
                                 n1Xn1_p1 = tensor(n1, n1, p1),
                                 n2Xn2_p2 = tensor(n2, n2, p2);
            // update the cell polarization matrix
            this->polarize_ = n0Xn0_p0 + n1Xn1_p1 + n2Xn2_p2;
            // cell has polarized, then proceed to check if it can do other things
            return true;
          }
        //...end of substances loop
      }
  // cell has not been through any polarization
  return false;
  //...end of cell polarization
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckProtrusion()
{
  if (!this->GetCanProtrude()) return false;
  // by design only viable (non-necrotic) cells could develop protrusions
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  if (! this->params()->have_parameter<int>(CP_name+"/can_protrude/pattern"))
    return false;
  //
  // confirn if permissible number of protrusions have been createdon the cell
  const int max_protrusions = this->params()->get<int>(CP_name+"/can_protrude/max_protrusions");
  if (max_protrusions == this->GetNumberOfProtrusions())
    return false;
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_protrude/probability"))
    return false;
  // now check if cell age is within an appropriate time-window to allow protrusion generation
  if ( this->GetAge() < this->params()->get<int>(CP_name+"/can_protrude/time_window_open" ) ||
       this->GetAge() > this->params()->get<int>(CP_name+"/can_protrude/time_window_close") )
    return false;
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_protrude/diameter_cutoff");
  //
  if ( diameter < diameter_cutoff )
    return false;
  //
  bool stimulate = false;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                   threshold = this->params()->get<double>(CP_name+"/can_protrude/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          stimulate = true;
        }
      //...end of substances loop
    }
  // check if any biochemical stimulus leads towards cell protrusion generation
  if ( ! stimulate ) return false;
  //
  const double protrusion_tol =
    ! this->params()->have_parameter<double>(CP_name+"/can_protrude/tolerance")
    ? tol :      this->params()->get<double>(CP_name+"/can_protrude/tolerance");
  if (!check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), protrusion_tol))
    return false;
  //
  CellProtrusion c_p;
  //
  const int pattern = this->params()->get<int>(CP_name+"/can_protrude/pattern");
  // cell protrusion max diameter
  const double dia_max = this->params()->get<double>(CP_name+"/can_protrude/diameter/max");
  // identify mode of simulation domain
  if (this->params()->get<bool>("simulation_domain_is_2D"))
    {
      if (0==pattern)
        {
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              bdm::Double3 axis = {rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0), 0.0};
              // if protrusion is not valid the simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis)) {
                ip--;
                continue;
              }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                if (behavior.size()!=1)
                  ABORT_("an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (-1==pattern || +1==pattern)
        {
          //
          const std::string& chemo_substance = this->params()->get<std::string>(CP_name+"/can_protrude/chemotaxis");
          // iterate for all substances
          std::vector<std::string>::const_iterator
            ci = std::find(substances.begin(), substances.end(), chemo_substance);
          if (substances.end()==ci)
            ABORT_("an internal error occurred");
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(*ci);
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              const double radius = rg->Uniform(0.0,this->GetDiameter()),
                           phi    = rg->Uniform(0.0,2.0*bdm::Math::kPi);
              bdm::Double3 pnt = {radius*cos(phi), radius*sin(phi), 0.0};
              pnt += this->GetPosition();
              //
              bdm::Double3 axis;
              dg->GetGradient(pnt, &axis);
              // note these details...
              axis *= pattern;
              if (!normalize(axis, axis)) continue;
              // if protrusion is not valid then simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  // // // // // // //ip--;
                  continue;
                  //...well skip this for the time being and check again some other time,
                  // in order to avoid a perpetual loop at this point; however, please do
                  // fix this issue in the future!
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                    const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                    protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (+2==pattern)
        {
          //
          ABORT_("an internal error occurred");
          //
        }
      else
        ABORT_("an exception is caught");
    }
  else
    {
      if (0==pattern)
        {
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              bdm::Double3 axis = {rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0)};
              // if protrusion is not valid the simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  ip--;
                  continue;
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (-1==pattern || +1==pattern)
        {
          //
          const std::string& chemo_substance = this->params()->get<std::string>(CP_name+"/can_protrude/chemotaxis");
          // iterate for all substances
          std::vector<std::string>::const_iterator
            ci = std::find(substances.begin(), substances.end(), chemo_substance);
          ASSERT_(substances.end()!=ci,"an internal error occurred");
          //
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(*ci);
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              const double radius = rg->Uniform(0.0,this->GetDiameter()),
                           phi    = rg->Uniform(0.0,2.0*bdm::Math::kPi),
                           theta  = rg->Uniform(0.0,bdm::Math::kPi);
              bdm::Double3 pnt = {radius*cos(phi)*sin(theta), radius*sin(phi)*sin(theta), radius*cos(theta)};
              pnt += this->GetPosition();
              //
              bdm::Double3 axis;
              dg->GetGradient(pnt, &axis);
              // note these details...
              axis *= pattern;
              if (!normalize(axis, axis)) continue;
              // if protrusion is not valid then simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  // // // // // // //ip--;
                  continue;
                  //...well skip this for the time being and check again some other time,
                  // in order to avoid a perpetual loop at this point; however, please do
                  // fix this issue in the future!
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (+2==pattern)
        {
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              const double phi   = rg->Uniform()*bdm::Math::kPi*2.0,
                           theta = rg->Uniform()*bdm::Math::kPi;
              const bdm::Double3 axis = { cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta) };
              // if protrusion is not valid then simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  ip--;
                  continue;
                  //...well skip this for the time being and check again some other time,
                  // in order to avoid a perpetual loop at this point; however, please do
                  // fix this issue in the future!
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else
        ABORT_("an exception is caught");
    }
  //
  if ( this->GetNumberOfProtrusions() )
    {
      this->can_grow_      = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_grow");
      this->can_divide_    = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_divide");
      this->can_migrate_   = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_migrate");
      this->can_transform_ = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_transform");
      this->can_polarize_  = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_polarize");
    }
  // cell has developed protrusions (filodia or/and neurites)
  return true;
  //...end of cell protrusion
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckGrowth()
{
  if (!this->GetCanGrow()) return false;
  // by design only viable (non-necrotic) cells could grow
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  if ( this->params()->have_parameter<double>(CP_name+"/intracellular/damage/growth_block") )
    {
      const double block_g = this->params()->get<double>(CP_name+"/intracellular/damage/growth_block");
      if (dna_damage_ > block_g) return false;
    }
  //
  // Contact inhibition of growth: high local crowding suppresses biomass accumulation.
  // This is distinct from contact inhibition of proliferation (G1/S checkpoint).
  // Only active when the parameter can_grow/crowding_threshold is provided.
  if (this->params()->have_parameter<double>(CP_name+"/can_grow/crowding_threshold"))
    {
      const double influence_ratio =
        this->params()->have_parameter<double>(CP_name+"/can_divide/influence_ratio")
        ? this->params()->get<double>(CP_name+"/can_divide/influence_ratio") : 2.0;
      const double occ = ComputeLocalOccupancyRatio(this->GetPosition(), influence_ratio);
      if (occ >= this->params()->get<double>(CP_name+"/can_grow/crowding_threshold"))
        return false;
    }
  //
  const double diameter = this->GetDiameter(),
               diameter_min = this->params()->get<double>(CP_name+"/diameter/min"),
               diameter_max = this->params()->get<double>(CP_name+"/diameter/max");
  //
  // check if cell size is within acceptable (upper limit) before proceeding with growth
  if (diameter >= diameter_max)
    return false;
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_grow/probability"))
    return false;
  //
  double volume_rate = 0.0;
  {
    const double diameter_rate = this->params()->get<double>(CP_name+"/can_grow/diameter_rate");
    //
    volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
  }
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      // skip this substance if no growth threshold is defined for it
      if (! this->params()->have_parameter<double>(CP_name+"/can_grow/"+BC_name+"/threshold"))
        continue;
      //
      const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                   threshold = this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          const double diameter_rate =
            this->params()->have_parameter<double>(CP_name+"/can_grow/"+BC_name+"/diameter_rate")
            ? this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/diameter_rate")
            : this->params()->get<double>(CP_name+"/can_grow/diameter_rate");
          //
          // allow cell growth controlled by a combination of two biochemical cues, therefore
          // cell development is dependent from another substance as well
          if (this->params()->have_parameter<bool>(CP_name+"/can_grow/"+BC_name+"/dependency") &&
              this->params()->get<bool>(CP_name+"/can_grow/"+BC_name+"/dependency"))
            {
              // iterate for all OTHER substances
              for ( std::vector<std::string>::const_iterator
                    cj=substances.begin(); cj!=substances.end(); cj++ )
                {
                  if ( cj == ci ) continue;
                  // access the BioDynaMo diffusion grid
                  auto* dg_other = rm->GetDiffusionGrid(*cj);
                  const std::string& BC_other_name = dg_other->GetContinuumName(); // biochemical name
                  //
                  if (! this->params()->have_parameter<double>(CP_name+"/can_grow/"+BC_name+"/dependency/"+BC_other_name+"/threshold"))
                    continue;
                  //
                  const double concentration_other = GetInterpolatedValue(dg_other, this->GetPosition(), this->params()),
                               threshold_other = this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/dependency/"+BC_other_name+"/threshold");
                  //
                  if ( ( threshold_other > 0.0 && concentration_other > +threshold_other ) ||
                       ( threshold_other < 0.0 && concentration_other < -threshold_other ) )
                    {
                      if (! this->params()->have_parameter<double>(CP_name+"/can_grow/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                        {
                          volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
                        }
                      else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                        {
                          volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
                        }
                    }
                  //...end of other substances loop
                }
            // ...if no probability is provided by user, then cell simply grows!
            }
          else if (! this->params()->have_parameter<double>(CP_name+"/can_grow/"+BC_name+"/probability"))
            {
              volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
            }
          // ...otherwise, check the likelihood for cell growth ;)
          else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/probability"))
            {
              volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
            }
        }
      //...end of substances loop
    }
  //
  if ( ( diameter < diameter_max && volume_rate ) ||
       ( diameter < diameter_min && volume_rate > 0.0 ) )
    {
      this->ChangeVolume(volume_rate);
      // cell has grown, then proceed to check if it can do other things
      return true;
    }
  // cell has not been through any growth
  return false;
  //...end of cell growth
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckTransformationAndDivision()
{
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could divide after they transform
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  if ( this->params()->have_parameter<double>(CP_name+"/intracellular/damage/division_block") )
    {
      const double block = this->params()->get<double>(CP_name+"/intracellular/damage/division_block");
      if (dna_damage_ > block) return false;
    }
  //
  // cell cannot divide (not at least with current BioDynaMo implementation)
  // if it has developed protrusions (filopodia or/and neurites)
  // however, in division followed by transformation, cell protrusion phenotype
  // has been implemented to switch into new one (after transformation)
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  const int n_div = this->GetNumberOfDivisions();
  //
  // call cannot divide more times than it should
  if (n_div >= this->params()->get<int>(CP_name+"/can_divide/max"))
    return false;
  //
  // Compute effective division probability with CAP modulation
  double p0 = this->params()->get<double>(CP_name+"/can_divide/probability");
  double p_eff = p0;
  // Apply CAP-induced modulation if parameters are present and CAP is configured
  if (this->params()->have_parameter<double>(CP_name+"/can_divide/CAP_sensitivity"))
    {
      double beta_cap = this->params()->get<double>(CP_name+"/can_divide/CAP_sensitivity");
      if (beta_cap > 0.0 && this->params()->have_parameter<double>("CAP/duration_h"))
        {
          double t_cap_h = this->params()->get<double>("CAP/duration_h");
          double t_cap_s = t_cap_h * 3600.0; // convert hours to seconds
          double k_cap = 60.0; // default saturation time in seconds
          if (this->params()->have_parameter<double>(CP_name+"/can_divide/CAP_saturation_time"))
            k_cap = this->params()->get<double>(CP_name+"/can_divide/CAP_saturation_time");
          // Apply bounded modulation: p_eff = p0 * [1 + beta * t/(K + t)]
          double modulation = 1.0 + beta_cap * (t_cap_s / (k_cap + t_cap_s));
          if (modulation > 2.0) modulation = 2.0; // prevent excessive increase (max 2x)
          p_eff = p0 * modulation;
        }
    }
  // Apply age increment if configured
  double p_age_increment = this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age");
  if (p_age_increment > 0.0)
    p_eff = p_eff + p_age_increment * this->GetAge();
  // Final probability check
  if (rg->Uniform(0.0,1.0) > p_eff)
    return false;
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_divide/diameter_cutoff");
  const int cell_maturity = (n_div+1) * this->params()->get<int>(CP_name+"/can_divide/time_window");
  //
  if ( diameter < diameter_cutoff || this->GetAge() < cell_maturity )
    return false;
  //
  // produce the separation vector
  const bdm::Double3 axis =
    { rg->Uniform(-1.0,+1.0) ,
      rg->Uniform(-1.0,+1.0) ,
      (this->params()->get<bool>("simulation_domain_is_2D") ? 0.0 : rg->Uniform(-1.0,+1.0)) };
  //
  const double volume_ratio = rg->Uniform(0.9,1.1);
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances if cell can
  // transform and then divide (symmetrically)
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      if (! this->params()->have_parameter<int>(CP_name+"/can_transform_and_divide/"+BC_name+"/new_phenotype"))
        continue;
      //
      const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                   threshold = this->params()->get<double>(CP_name+"/can_transform_and_divide/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_transform_and_divide/"+BC_name+"/probability"))
            {
              const int new_phenotype = this->params()->get<int>(CP_name+"/can_transform_and_divide/"+BC_name+"/new_phenotype");
              // firstly, the cell transforms
              this->SetPhenotype(new_phenotype);
              // now reset the age of the cell
              this->SetAge();
              // increment this index
              this->IncrementNumberOfTrasformations();
              //
              const std::string CP_new_name = // cell phenotype name
                this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
              // principal directions of the cell polarization matrix
              double p0, p1, p2;
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                {
                  p0 = this->params()->get<double>(CP_new_name+"/principal/0");
                  p1 = this->params()->get<double>(CP_new_name+"/principal/1");
                  p2 = this->params()->get<double>(CP_new_name+"/principal/2");
                }
              //
              this->SetCanApoptose(this->params()->get<bool>(CP_new_name+"/can_apoptose"));
              this->SetCanGrow(this->params()->get<bool>(CP_new_name+"/can_grow"));
              this->SetCanDivide(this->params()->get<bool>(CP_new_name+"/can_divide"));
              this->SetCanMigrate(this->params()->get<bool>(CP_new_name+"/can_migrate"));
              this->SetCanTransform(this->params()->get<bool>(CP_new_name+"/can_transform"));
              this->SetCanProtrude(this->params()->get<bool>(CP_new_name+"/can_protrude"));
              this->SetCanPolarize(this->params()->get<bool>(CP_new_name+"/can_polarize"));
              // reset the cell polarization matrix
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                this->SetPolarization(diag(p0, p1, p2));
              // reset the cell protrusion phenotype
              if ( this->GetNumberOfProtrusions() )
                {
                  if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
                    ABORT_("an internal error occurred");
                  //
                  // iterate for all (existing) protrusions of this cell
                  for (int p=0; p<this->GetNumberOfProtrusions(); p++)
                    {
                      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
                      // assign this cell (that is associated with) to the protrusion created
                      protrusion->SetCell(this);
                    }
                }
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = this->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an exception is caught");
                //
                this->RemoveBehavior(behavior[0]);
              }
              const int mo = this->params()->get<int>(CP_new_name+"/mechanism_order");
              if      (10==mo) this->AddBehavior(new Biology4BiologicalCell_10());
              else if (11==mo) this->AddBehavior(new Biology4BiologicalCell_11());
              else if (12==mo) this->AddBehavior(new Biology4BiologicalCell_12());
              else             ABORT_("an exception is caught");
              // secondly, the cell divides
              this->Divide(volume_ratio, axis);
              // cell has transformed and divided, then proceed to check if it can do other things
              return true;
            }
        }
      //...end of substances loop
    }
  // cell has not been through any division
  return false;
  //...end of cell division
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckAsymmetricDivision()
{
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could divide and then transform
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // DNA damage division block gate (centralized checkpoint enforcement)
  if (this->params()->have_parameter<double>(CP_name+"/intracellular/damage/division_block"))
    {
      const double block = this->params()->get<double>(CP_name+"/intracellular/damage/division_block");
      if (dna_damage_ > block) return false;
    }
  //
  // cell cannot divide (not at least with current BioDynaMo implementation)
  // if it has developed protrusions (filopodia or/and neurites)
  // however, in division followed by transformation, cell protrusion phenotype
  // has been implemented to switch into new one (after transformation)
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  const int n_div = this->GetNumberOfDivisions();
  //
  // call cannot divide more times than it should
  if (n_div >= this->params()->get<int>(CP_name+"/can_divide/max"))
    return false;
  //
  // Compute effective division probability with CAP modulation
  double p0 = this->params()->get<double>(CP_name+"/can_divide/probability");
  double p_eff = p0;
  // Apply CAP-induced modulation if parameters are present and CAP is configured
  if (this->params()->have_parameter<double>(CP_name+"/can_divide/CAP_sensitivity"))
    {
      double beta_cap = this->params()->get<double>(CP_name+"/can_divide/CAP_sensitivity");
      if (beta_cap > 0.0 && this->params()->have_parameter<double>("CAP/duration_h"))
        {
          double t_cap_h = this->params()->get<double>("CAP/duration_h");
          double t_cap_s = t_cap_h * 3600.0; // convert hours to seconds
          double k_cap = 60.0; // default saturation time in seconds
          if (this->params()->have_parameter<double>(CP_name+"/can_divide/CAP_saturation_time"))
            k_cap = this->params()->get<double>(CP_name+"/can_divide/CAP_saturation_time");
          // Apply bounded modulation: p_eff = p0 * [1 + beta * t/(K + t)]
          double modulation = 1.0 + beta_cap * (t_cap_s / (k_cap + t_cap_s));
          if (modulation > 2.0) modulation = 2.0; // prevent excessive increase (max 2x)
          p_eff = p0 * modulation;
        }
    }
  // Apply age increment if configured
  double p_age_increment = this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age");
  if (p_age_increment > 0.0)
    p_eff = p_eff + p_age_increment * this->GetAge();
  // Final probability check
  if (rg->Uniform(0.0,1.0) > p_eff)
    return false;
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_divide/diameter_cutoff");
  const int cell_maturity = (n_div+1) * this->params()->get<int>(CP_name+"/can_divide/time_window");
  //
  if ( diameter < diameter_cutoff || this->GetAge() < cell_maturity )
    return false;
  //
  // produce the separation vector
  const bdm::Double3 axis =
    { rg->Uniform(-1.0,+1.0) ,
      rg->Uniform(-1.0,+1.0) ,
      (this->params()->get<bool>("simulation_domain_is_2D") ? 0.0 : rg->Uniform(-1.0,+1.0)) };
  //
  const double volume_ratio = rg->Uniform(0.9,1.1);
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances if cell can
  // divide and then transform
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      if (! this->params()->have_parameter<double>(CP_name+"/can_divide_and_transform/"+BC_name+"/threshold"))
        continue;
      //
      const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                   threshold = this->params()->get<double>(CP_name+"/can_divide_and_transform/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_divide_and_transform/"+BC_name+"/probability"))
            {
              const int new_phenotype = this->params()->get<int>(CP_name+"/can_divide_and_transform/"+BC_name+"/new_phenotype");
              // firstly, the cell divides
              this->Divide(volume_ratio, axis);
              // secondly, the cell transforms
              this->SetPhenotype(new_phenotype);
              // now reset the age of the cell
              this->SetAge();
              // increment this index
              this->IncrementNumberOfTrasformations();
              //
              const std::string CP_new_name = // cell phenotype name
                this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
              // principal directions of the cell polarization matrix
              double p0, p1, p2;
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                {
                  p0 = this->params()->get<double>(CP_new_name+"/principal/0");
                  p1 = this->params()->get<double>(CP_new_name+"/principal/1");
                  p2 = this->params()->get<double>(CP_new_name+"/principal/2");
                }
              //
              this->SetCanApoptose(this->params()->get<bool>(CP_new_name+"/can_apoptose"));
              this->SetCanGrow(this->params()->get<bool>(CP_new_name+"/can_grow"));
              this->SetCanDivide(this->params()->get<bool>(CP_new_name+"/can_divide"));
              this->SetCanMigrate(this->params()->get<bool>(CP_new_name+"/can_migrate"));
              this->SetCanTransform(this->params()->get<bool>(CP_new_name+"/can_transform"));
              this->SetCanProtrude(this->params()->get<bool>(CP_new_name+"/can_protrude"));
              this->SetCanPolarize(this->params()->get<bool>(CP_new_name+"/can_polarize"));
              // reset the cell polarization matrix
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                this->SetPolarization(diag(p0, p1, p2));
              // reset the cell protrusion phenotype
              if ( this->GetNumberOfProtrusions() )
                {
                  if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
                    ABORT_("an internal error occurred");
                  // iterate for all (existing) protrusions of this cell
                  for (int p=0; p<this->GetNumberOfProtrusions(); p++)
                    {
                      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
                      // assign this cell (that is associated with) to the protrusion created
                      protrusion->SetCell(this);
                    }
                }
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = this->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                this->RemoveBehavior(behavior[0]);
              }
              const int mo = this->params()->get<int>(CP_new_name+"/mechanism_order");
              if      (10==mo) this->AddBehavior(new Biology4BiologicalCell_10());
              else if (11==mo) this->AddBehavior(new Biology4BiologicalCell_11());
              else if (12==mo) this->AddBehavior(new Biology4BiologicalCell_12());
              else             ABORT_("an exception is caught");
              // cell has divided and transformed, then proceed to check if it can do other things
              return true;
            }
        }
      //...end of substances loop
    }
  // cell has not been through any division
  return false;
  //...end of cell division
}
// -----------------------------------------------------------------------------
inline
double bdm::BiologicalCell::ComputeLocalOccupancyRatio(
  const bdm::Double3& position, double influence_ratio) const
{
  if (influence_ratio <= 0.0) return 0.0;
  //
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  const bool is_2D = this->params()->get<bool>("simulation_domain_is_2D");
  const real_t R = this->GetDiameter() * static_cast<real_t>(influence_ratio);
  //
  if (is_2D)
    {
      const real_t A = bdm::Math::kPi * pow2(R);
      if (A <= 0.0) return 0.0;
      real_t area(0.0);
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* other = dynamic_cast<const BiologicalCell*>(a))
          if (other != this)
            {
              const real_t r = 0.5 * other->GetDiameter();
              const real_t d = L2norm(position - other->GetPosition());
              if      (d > (R+r))      ;
              else if (d <= abs(R-r))  area += bdm::Math::kPi * pow2(r);
              else
                area += pow2(r)*acos((pow2(d)+pow2(r)-pow2(R))/(2*d*r))
                      + pow2(R)*acos((pow2(d)-pow2(r)+pow2(R))/(2*d*R))
                      - 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
            }
      });
      return static_cast<double>(area / A);
    }
  else
    {
      const real_t V = bdm::Math::kPi * pow3(R) * (4.0/3.0);
      if (V <= 0.0) return 0.0;
      real_t volume(0.0);
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* other = dynamic_cast<const BiologicalCell*>(a))
          if (other != this)
            {
              const real_t r = 0.5 * other->GetDiameter();
              const real_t d = L2norm(position - other->GetPosition());
              if      (d > (R+r))      ;
              else if (d <= abs(R-r))  volume += bdm::Math::kPi * pow3(r) * (4.0/3.0);
              else
                volume += bdm::Math::kPi * pow2(R+r-d) / (12.0*d)
                        * (d*d + 2*d*r - 3*r*r + 2*d*R - 3*R*R + 6*r*R);
            }
      });
      return static_cast<double>(volume / V);
    }
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckDivision() {
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could divide
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // DNA damage division block gate (centralized checkpoint enforcement)
  if (this->params()->have_parameter<double>(CP_name+"/intracellular/damage/division_block"))
    {
      const double block = this->params()->get<double>(CP_name+"/intracellular/damage/division_block");
      if (dna_damage_ > block) return false;
    }
  //
  // cell cannot divide (not at least with current BioDynaMo implementation)
  // if it has developed protrusions (filopodia or/and neurites)
  // however, in division followed by transformation, cell protrusion phenotype
  // has been implemented to switch into new one (after transformation)
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  const int n_div = this->GetNumberOfDivisions();
  //
  // call cannot divide more times than it should
  if (n_div >= this->params()->get<int>(CP_name+"/can_divide/max"))
    return false;
  //
  // Compute effective division probability with CAP modulation
  double p0 = this->params()->get<double>(CP_name+"/can_divide/probability");
  double p_eff = p0;
  // Apply CAP-induced modulation if parameters are present and CAP is configured
  if (this->params()->have_parameter<double>(CP_name+"/can_divide/CAP_sensitivity"))
    {
      double beta_cap = this->params()->get<double>(CP_name+"/can_divide/CAP_sensitivity");
      if (beta_cap > 0.0 && this->params()->have_parameter<double>("CAP/duration_h"))
        {
          double t_cap_h = this->params()->get<double>("CAP/duration_h");
          double t_cap_s = t_cap_h * 3600.0; // convert hours to seconds
          double k_cap = 60.0; // default saturation time in seconds
          if (this->params()->have_parameter<double>(CP_name+"/can_divide/CAP_saturation_time"))
            k_cap = this->params()->get<double>(CP_name+"/can_divide/CAP_saturation_time");
          // Apply bounded modulation: p_eff = p0 * [1 + beta * t/(K + t)]
          double modulation = 1.0 + beta_cap * (t_cap_s / (k_cap + t_cap_s));
          if (modulation > 2.0) modulation = 2.0; // prevent excessive increase (max 2x)
          p_eff = p0 * modulation;
        }
    }
  // Apply age increment if configured
  double p_age_increment = this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age");
  if (p_age_increment > 0.0)
    p_eff = p_eff + p_age_increment * this->GetAge();
  // Final probability check
  if (rg->Uniform(0.0,1.0) > p_eff)
    return false;
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_divide/diameter_cutoff");
  const int cell_maturity = (n_div+1) * this->params()->get<int>(CP_name+"/can_divide/time_window");
  //
  if ( diameter < diameter_cutoff || this->GetAge() < cell_maturity )
    return false;
  //
  // local-crowding inhibition of division
  {
    const double influence_ratio = this->params()->get<double>(CP_name+"/can_divide/influence_ratio");
    if (influence_ratio > 0.0)
      {
        const double max_occupancy = this->params()->get<double>(CP_name+"/can_divide/max_occupancy");
        if (ComputeLocalOccupancyRatio(this->GetPosition(), influence_ratio) >= max_occupancy)
          return false;
      }
  }
  //
  // produce the separation vector
  const bdm::Double3 axis =
    { rg->Uniform(-1.0,+1.0) ,
      rg->Uniform(-1.0,+1.0) ,
      (this->params()->get<bool>("simulation_domain_is_2D") ? 0.0 : rg->Uniform(-1.0,+1.0)) };
  //
  const double volume_ratio = rg->Uniform(0.9,1.1);
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  // Check if any substance-based division thresholds are defined
  bool any_substance_threshold_defined = false;
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName();
      if (this->params()->have_parameter<double>(CP_name+"/can_divide/"+BC_name+"/threshold"))
        {
          any_substance_threshold_defined = true;
          break;
        }
    }
  //
  // If no substance thresholds are defined, allow division based on basic criteria
  if (! any_substance_threshold_defined)
    {
      this->Divide(volume_ratio, axis);
      return true;
    }
  //
  // Otherwise, check if at least one substance condition is met
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances if cell can
  // divide (symmetrically)
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      if (! this->params()->have_parameter<double>(CP_name+"/can_divide/"+BC_name+"/threshold"))
        continue;
      //
      const double concentration = GetInterpolatedValue(dg, this->GetPosition(), this->params()),
                   threshold = this->params()->get<double>(CP_name+"/can_divide/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          // since no symmetric (prior to cell transformation) or unsymmetric division
          // has occurred, then cell divides conventionally
          this->Divide(volume_ratio, axis);
          // cell has divided, then proceed to check if it can do other things
          return true;
        }
      //...end of substances loop
    }
  // cell has not been through any division
  return false;
  //...end of cell division
}
// -----------------------------------------------------------------------------
inline
void bdm::BiologicalCell::CheckAndFixDiameter()
{
  if (!this->GetPhenotype()) return;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const double diameter = this->GetDiameter(),
               diameter_min = this->params()->get<double>(CP_name+"/diameter/min"),
               diameter_max = this->params()->get<double>(CP_name+"/diameter/max");
  // check if cell diameter is within user-defined bounds
  if      (diameter < diameter_min)
    this->SetDiameter(diameter_min);
  else if (diameter > diameter_max)
    this->SetDiameter(diameter_max);
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckProtrusionAxis(bdm::Double3 axis)
{
  ASSERT_(0!=this->GetPhenotype(),"an internal error occurred");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // relative angle between protrusions
  const double rel_angle_min = this->params()->get<double>(CP_name+"/can_protrude/relative_angle/degrees/min"),
               rel_angle_max = this->params()->get<double>(CP_name+"/can_protrude/relative_angle/degrees/max");
  // make sure to normalize the axis vector first
  ASSERT_(normalize(axis, axis),"could not normalize the axis vector");
  //
  if (!this->protrusions_.empty())
    {
      for (unsigned int l=0; l<this->protrusions_.size(); l++)
        {
          const bdm::Double3& current_axis = this->protrusions_[l];
          const double angle = radians_to_degrees( acos(axis*current_axis) );
          // check if relative angle is within user-defined range
          if (angle<rel_angle_min || angle>rel_angle_max)
            return false;
        }
    }
  //
  this->protrusions_.push_back(axis);
  return true;
}
// -----------------------------------------------------------------------------
inline
void bdm::BiologicalCell::Set2DeleteProtrusions()
{
  // sanity check
  if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
    ABORT_("an internal error occurred");
  // iterate for all (existing) protrusions of this cell
  for (int p=0; p<this->GetNumberOfProtrusions(); p++)
    {
      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
      // assign this cell (that is associated with) to the protrusion created
      protrusion->Set2Delete();
    }
}
// -----------------------------------------------------------------------------
inline
double bdm::BiologicalCell::GetRegulatoryNodeActivity(
    const std::string& node_name) const
{
  auto it = regulatory_state_.node_activity.find(node_name);
  return (it == regulatory_state_.node_activity.end()) ? 0.0 : it->second;
}
// -----------------------------------------------------------------------------
inline
void bdm::BiologicalCell::ConfigureRegulatoryModel()
{
  regulatory_backend_id_ = 0;
  regulatory_model_type_ = "none";
  regulatory_update_interval_ = 1;
  regulatory_update_counter_ = 0;
  regulatory_model_phenotype_id_ = this->GetPhenotype();
  regulatory_state_.node_activity.clear();
  regulatory_state_.time_since_last_update = 0.0;
  regulatory_state_.last_phase = this->GetPhase();
  regulatory_output_ = bdm::regulatory::RegulatoryOutput();
  regulatory_parameters_ = bdm::regulatory::RegulatoryParameters();
  if (!params_ || !this->GetPhenotype()) return;
  // phenotype namespace
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  // model selection
  std::string model_type = "none";
  if (this->params()->have_parameter<std::string>(CP_name+"/intracellular/model_type")) {
    model_type = this->params()->get<std::string>(CP_name+"/intracellular/model_type");
  }
  if (model_type == "deterministic_boolean") model_type = "boolean";
  if (model_type != "none" && model_type != "boolean" &&
      model_type != "stochastic_boolean") {
    model_type = "none";
  }
  regulatory_model_type_ = model_type;
  if (model_type == "none") return;
  // update interval
  if (this->params()->have_parameter<int>(CP_name+"/intracellular/grn_update_interval")) {
    regulatory_update_interval_ =
      std::max(1, this->params()->get<int>(CP_name+"/intracellular/grn_update_interval"));
  }
  regulatory_parameters_.backend = model_type;
  regulatory_parameters_.update_interval = regulatory_update_interval_;
  regulatory_parameters_.boolean["noise"] =
    this->params()->have_parameter<bool>(CP_name+"/intracellular/grn_noise")
    ? this->params()->get<bool>(CP_name+"/intracellular/grn_noise") : false;
  // compact helper for per-phenotype GRN parameters
  auto read_grn = [&] (const std::string& key, double fallback) -> double {
    const std::string pname = CP_name + "/grn/" + key;
    return this->params()->have_parameter<double>(pname)
      ? this->params()->get<double>(pname) : fallback;
  };
  // thresholds / scales
  regulatory_parameters_.numeric["ros_high_threshold"] = read_grn("ros_high_threshold", 0.35);
  regulatory_parameters_.numeric["dna_damage_high_threshold"] = read_grn("dna_damage_high_threshold", 0.35);
  regulatory_parameters_.numeric["hypoxia_threshold"] = read_grn("hypoxia_threshold", 0.15);
  regulatory_parameters_.numeric["nutrient_low_threshold"] = read_grn("nutrient_low_threshold", 0.25);
  regulatory_parameters_.numeric["ecm_stiffness_high_threshold"] = read_grn("ecm_stiffness_high_threshold", 1.2);
  regulatory_parameters_.numeric["crowding_high_threshold"] = read_grn("crowding_high_threshold", 0.6);
  regulatory_parameters_.numeric["tgfb_high_threshold"] = read_grn("tgfb_high_threshold", 0.4);
  regulatory_parameters_.numeric["adhesion_low_threshold"] = read_grn("adhesion_low_threshold", 0.2);
  regulatory_parameters_.numeric["apoptosis_hazard_scale"] = read_grn("apoptosis_hazard_scale", 0.15);
  regulatory_parameters_.numeric["necrosis_hazard_scale"] = read_grn("necrosis_hazard_scale", 0.05);
  regulatory_parameters_.numeric["quiescence_hazard_scale"] = read_grn("quiescence_hazard_scale", 0.08);
  regulatory_parameters_.numeric["apoptosis_hazard_base"] = read_grn("apoptosis_hazard_base", 0.0);
  regulatory_parameters_.numeric["necrosis_hazard_base"] = read_grn("necrosis_hazard_base", 0.0);
  regulatory_parameters_.numeric["quiescence_hazard_base"] = read_grn("quiescence_hazard_base", 0.0);
  // stochastic transition rates
  regulatory_parameters_.numeric["rate_on_p53"] = read_grn("p53_on_rate", 1.2);
  regulatory_parameters_.numeric["rate_off_p53"] = read_grn("p53_off_rate", 0.4);
  regulatory_parameters_.numeric["rate_on_NRF2"] = read_grn("nrf2_on_rate", 1.0);
  regulatory_parameters_.numeric["rate_off_NRF2"] = read_grn("nrf2_off_rate", 0.5);
  regulatory_parameters_.numeric["rate_on_Caspase3"] = read_grn("caspase_on_rate", 1.0);
  regulatory_parameters_.numeric["rate_off_Caspase3"] = read_grn("caspase_off_rate", 0.4);
  // initialize selected backend
  if (model_type == "boolean") {
    regulatory_backend_id_ = 1;
    regulatory_model_boolean_.Initialize(regulatory_parameters_);
  } else if (model_type == "stochastic_boolean") {
    regulatory_backend_id_ = 2;
    regulatory_model_stochastic_.Initialize(regulatory_parameters_);
  }
}
// -----------------------------------------------------------------------------
inline
bdm::regulatory::RegulatoryInput
bdm::BiologicalCell::BuildRegulatoryInput(
    const bdm::BiologicalCell::MicroenvironmentState& env) const
{
  bdm::regulatory::RegulatoryInput in;
  in.oxygen = env.local_O2;
  in.nutrient = env.local_nutrient;
  in.extracellular_rons = env.local_rons;
  in.intracellular_ros = ros_internal_;
  in.dna_damage = dna_damage_;
  in.ecm_density = env.ecm_density;
  in.ecm_stiffness = env.ecm_stiffness;
  in.adhesion_signal = env.ecm_adhesion;
  in.crowding = env.local_crowding;
  in.phase = this->GetPhase();
  in.phenotype_id = this->GetPhenotype();
  // Optional external cues
  if (params_ && this->params()->have_parameter<std::vector<std::string>>("substances")) {
    auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
    const bdm::Double3 xyz = this->GetPosition();
    const auto& substances = this->params()->get<std::vector<std::string>>("substances");
    if (std::find(substances.begin(), substances.end(), "TGFb") != substances.end())
      if (auto* dg = rm->GetDiffusionGrid("TGFb"))
        in.tgfb = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
    if (std::find(substances.begin(), substances.end(), "IL6") != substances.end())
      if (auto* dg = rm->GetDiffusionGrid("IL6"))
        in.inflammatory_signal = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
    if (std::find(substances.begin(), substances.end(), "EGF") != substances.end())
      if (auto* dg = rm->GetDiffusionGrid("EGF"))
        in.growth_factor = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
  }
  return in;
}
// -----------------------------------------------------------------------------
inline
void bdm::BiologicalCell::UpdateRegulatoryModel(
    const bdm::BiologicalCell::MicroenvironmentState& env)
{
  if (!this->GetPhenotype() || !params_) return;
  if (regulatory_model_phenotype_id_ != this->GetPhenotype()) {
    ConfigureRegulatoryModel();
  }
  if (regulatory_backend_id_ == 0) return;
  const double dt = this->params()->get<double>("time_step");
  ++regulatory_update_counter_;
  const bool phase_changed = (regulatory_state_.last_phase != this->GetPhase());
  const bool due_update = (regulatory_update_counter_ % regulatory_update_interval_) == 0;
  if (!phase_changed && !due_update) {
    return;
  }
  const bdm::regulatory::RegulatoryInput in = BuildRegulatoryInput(env);
  if (regulatory_backend_id_ == 1) {
    regulatory_model_boolean_.Update(in, &regulatory_state_, &regulatory_output_, dt);
  } else if (regulatory_backend_id_ == 2) {
    regulatory_model_stochastic_.Update(in, &regulatory_state_, &regulatory_output_, dt);
  } else {
    regulatory_model_null_.Update(in, &regulatory_state_, &regulatory_output_, dt);
  }
  regulatory_state_.last_phase = this->GetPhase();
  // Map GRN outputs to intracellular capacities (smoothed to avoid abrupt jumps)
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const double alpha = this->params()->have_parameter<double>(CP_name+"/intracellular/grn_relaxation")
    ? std::clamp(this->params()->get<double>(CP_name+"/intracellular/grn_relaxation"), 0.0, 1.0)
    : 0.5;
  SetAntioxidantCapacity((1.0 - alpha) * antioxidant_capacity_
                         + alpha * std::max(0.0, regulatory_output_.antioxidant_capacity));
  SetRepairCapacity((1.0 - alpha) * repair_capacity_
                    + alpha * std::clamp(regulatory_output_.repair_capacity, 0.0, 1.0));
}
// =============================================================================
// Mechanism 11 supporting methods — added for biologically stricter control
// =============================================================================
inline
bdm::BiologicalCell::MicroenvironmentState
bdm::BiologicalCell::SampleMicroenvironment() const
{
  // Samples all relevant local fields once per timestep to avoid repeated
  // grid lookups. Returns safe defaults for absent substances so existing
  // input files run unchanged.
  MicroenvironmentState st;
  if (!this->GetPhenotype()) return st; // necrotic cells: return defaults
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  const auto& substances = this->params()->get<std::vector<std::string>>("substances");
  const bdm::Double3 xyz = this->GetPosition();
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // --- Oxygen (O2) ---
  if (std::find(substances.begin(), substances.end(), "O2") != substances.end())
    if (auto* dg = rm->GetDiffusionGrid("O2"))
      st.local_O2 = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
  //
  // --- Nutrient / glucose (first matching substance name) ---
  for (const char* nname : {"nutrient", "Nutrient", "glucose", "Glucose"}) {
    const std::string ns(nname);
    if (std::find(substances.begin(), substances.end(), ns) != substances.end())
      if (auto* dg = rm->GetDiffusionGrid(ns)) {
        st.local_nutrient = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
        break;
      }
  }
  //
  // --- ECM density (collagen/fibronectin/laminin scaffold) ---
  if (std::find(substances.begin(), substances.end(), "ECM") != substances.end())
    if (auto* dg = rm->GetDiffusionGrid("ECM"))
      st.ecm_density = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
  //
  // --- ECM stiffness (YAP/TAZ mechanosensing modifier) ---
  // Stub: use "ECM_stiffness" substance if present; default = 1.0 (normal).
  // Biologically: matrix stiffness measured in kPa activates integrin clusters,
  // FAK/Src, RhoA/ROCK, actomyosin tension, and nuclear YAP/TAZ localisation.
  if (std::find(substances.begin(), substances.end(), "ECM_stiffness") != substances.end())
    if (auto* dg = rm->GetDiffusionGrid("ECM_stiffness"))
      st.ecm_stiffness = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
  //
  // --- ECM adhesion ligand density (integrin-binding sites) ---
  // Stub: use "ECM_adhesion" substance if present.
  // Fallback: derive from ecm_density × integrin_sensitivity (phenotype param).
  // Biologically: fibronectin/laminin/vitronectin RGD motifs bind α5β1/αvβ3;
  // low density → detachment → anoikis (BIM/BAD activation).
  if (std::find(substances.begin(), substances.end(), "ECM_adhesion") != substances.end()) {
    if (auto* dg = rm->GetDiffusionGrid("ECM_adhesion"))
      st.ecm_adhesion = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
  } else {
    // Fallback: adhesion signal = ECM density × integrin_sensitivity
    const double intsens =
      this->params()->have_parameter<double>(CP_name+"/ecm/integrin_sensitivity")
      ? this->params()->get<double>(CP_name+"/ecm/integrin_sensitivity") : 1.0;
    st.ecm_adhesion = std::max(0.0, st.ecm_density * intsens);
  }
  //
  // --- Intracellular ROS/RONS stress level ---
  // Represents combined oxidative burden from H2O2/NO2 uptake and
  // mitochondrial superoxide. Used to gate growth and division.
  st.local_rons = ros_internal_;
  //
  // --- Local cell-volume crowding (occupancy ratio) ---
  // Uses can_divide/influence_ratio parameter; default 2.0.
  {
    const double influence_ratio =
      this->params()->have_parameter<double>(CP_name+"/can_divide/influence_ratio")
      ? this->params()->get<double>(CP_name+"/can_divide/influence_ratio") : 2.0;
    st.local_crowding = ComputeLocalOccupancyRatio(xyz, influence_ratio);
  }
  return st;
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::EvaluateIntraSCheckpoint()
{
  // Intra-S checkpoint: returns true if Sy→G2 transition should be BLOCKED.
  // Mechanism: persistent replication stress or unrepaired DNA DSBs during
  // DNA synthesis activate ATR → CHK1 → CDC25A ubiquitination (degradation)
  // → CDK2 inhibition → S-phase arrest.
  // Only active when the relevant threshold parameters are configured.
  if (!this->GetPhenotype()) return false;
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // DNA damage intra-S threshold (ATR-mediated; e.g. stalled replication forks)
  if (this->params()->have_parameter<double>(CP_name+"/checkpoint/IntraS/damage_threshold")) {
    const double thr = this->params()->get<double>(CP_name+"/checkpoint/IntraS/damage_threshold");
    if (dna_damage_ > thr) return true;
  }
  //
  // ROS/replication stress threshold (oxidative damage to replication machinery)
  if (this->params()->have_parameter<double>(CP_name+"/checkpoint/IntraS/ros_threshold")) {
    const double thr = this->params()->get<double>(CP_name+"/checkpoint/IntraS/ros_threshold");
    if (ros_internal_ > thr) return true;
  }
  return false;
}
// -----------------------------------------------------------------------------
inline
bdm::BiologicalCell::CellCycleCheckpointState
bdm::BiologicalCell::EvaluateCellCycleCheckpoints(
    const bdm::BiologicalCell::MicroenvironmentState& env)
{
  // Central checkpoint controller for Mechanism 11.
  // Uses pre-sampled MicroenvironmentState to avoid re-querying diffusion grids.
  // All gates use have_parameter<> guards so absent parameters default to
  // permissive (no blocking), preserving backward compatibility.
  CellCycleCheckpointState st; // all fields default-initialized to permissive
  if (!this->GetPhenotype()) return st; // necrotic cells: no active checkpoints
  //
  auto* params = this->params();
  const std::string& CP_name =
    params->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // ================================================================
  // G1/S checkpoint
  // Molecular annotation (for parameter naming, not GRN):
  //   ATM/ATR → CHK1/CHK2 → p53 → p21 → CDK2/CyclinE inhibition
  //   CDC25A inhibition → CDK2 dephosphorylated → RB hypophosphorylated
  //   → E2F transcription factors blocked → S-phase genes suppressed
  // ================================================================
  // Molecular DDR cascade (ATM/ATR–CHK–p53–p21–Cdc25–CDK) takes precedence
  if (bdm::IsMolecularG1SCheckpointBlocked(this)) {
    st.can_enter_S = false;
  } else {
    // Legacy aggregate DNA damage threshold (when DDR pathway not enabled)
    if (!bdm::IsDdrPathwayEnabled(*params, CP_name) &&
        params->have_parameter<double>(CP_name+"/checkpoint/G1S/damage_threshold")) {
      if (dna_damage_ > params->get<double>(CP_name+"/checkpoint/G1S/damage_threshold"))
        st.can_enter_S = false;
    }
  }
  // O2 gate: hypoxia suppresses E2F via HIF-1α / CDK inhibitor accumulation
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/O2_threshold") &&
      env.local_O2 < params->get<double>(CP_name+"/checkpoint/G1S/O2_threshold"))
    st.can_enter_S = false;
  // Nutrient gate: AMPK activation under nutrient stress → p21/p27 upregulation
  //   → CDK2/CyclinE inhibition → G1 arrest (mTORC1 nutrient sensing)
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/nutrient_min") &&
      env.local_nutrient < params->get<double>(CP_name+"/checkpoint/G1S/nutrient_min"))
    st.can_enter_S = false;
  // ECM density gate: sparse matrix → poor integrin signalling → G1 arrest
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/ECM_threshold") &&
      env.ecm_density < params->get<double>(CP_name+"/checkpoint/G1S/ECM_threshold"))
    st.can_enter_S = false;
  // ECM adhesion gate: low adhesion → integrin/FAK loss → RB hypophosphorylation
  //   Mechanistically distinct from anoikis: promotes reversible G1 arrest
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/adhesion_min") &&
      env.ecm_adhesion < params->get<double>(CP_name+"/checkpoint/G1S/adhesion_min"))
    st.can_enter_S = false;
  // Crowding gate: contact inhibition of proliferation via E-cadherin/RhoA/ROCK
  //   → LATS1/2 kinase activation → YAP/TAZ cytoplasmic sequestration → CIP/KIP
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/crowding_threshold") &&
      env.local_crowding >= params->get<double>(CP_name+"/checkpoint/G1S/crowding_threshold"))
    st.can_enter_S = false;
  //
  // ================================================================
  // G2/M checkpoint
  // ATM/ATR → CHK1/CHK2 → CDC25C inhibition (14-3-3σ sequestration)
  //   → CDK1/CyclinB dephosphorylated/inactive → mitosis entry blocked
  // ================================================================
  if (bdm::IsMolecularG2MCheckpointBlocked(this)) {
    st.can_enter_M = false;
  } else {
    if (!bdm::IsDdrPathwayEnabled(*params, CP_name) &&
        params->have_parameter<double>(CP_name+"/checkpoint/G2M/damage_threshold")) {
      if (dna_damage_ > params->get<double>(CP_name+"/checkpoint/G2M/damage_threshold"))
        st.can_enter_M = false;
    }
  }
  // O2 gate for G2/M: cells need ATP for spindle assembly and chromosome segregation
  if (st.can_enter_M &&
      params->have_parameter<double>(CP_name+"/checkpoint/G2M/O2_threshold") &&
      env.local_O2 < params->get<double>(CP_name+"/checkpoint/G2M/O2_threshold"))
    st.can_enter_M = false;
  //
  // Composite: can_divide requires both G2/M gate cleared AND G1/S gate cleared
  // (can_enter_S false means cell never reached Di phase via checkpoint)
  if (!st.can_enter_M) st.can_divide = false;
  //
  // ================================================================
  // Growth gate
  // ================================================================
  // Nutrient: PI3K/Akt/mTOR → S6K/4EBP1 → protein synthesis
  if (params->have_parameter<double>(CP_name+"/checkpoint/growth/nutrient_min") &&
      env.local_nutrient < params->get<double>(CP_name+"/checkpoint/growth/nutrient_min"))
    st.can_grow = false;
  // O2: oxidative phosphorylation needed for biosynthetic energy
  if (st.can_grow &&
      params->have_parameter<double>(CP_name+"/checkpoint/growth/O2_min") &&
      env.local_O2 < params->get<double>(CP_name+"/checkpoint/growth/O2_min"))
    st.can_grow = false;
  // Crowding: contact inhibition of growth (E-cadherin → LATS/YAP/TAZ)
  if (st.can_grow &&
      params->have_parameter<double>(CP_name+"/can_grow/crowding_threshold") &&
      env.local_crowding >= params->get<double>(CP_name+"/can_grow/crowding_threshold"))
    st.can_grow = false;
  //
  // ================================================================
  // G0/quiescence recommendation
  // p27/Kip1 upregulation; RB hypophosphorylation; reversible arrest
  // ================================================================
  {
    const double crowd_entry =
      params->have_parameter<double>(CP_name+"/quiescence/crowding_threshold")
      ? params->get<double>(CP_name+"/quiescence/crowding_threshold") : -1.0;
    const double o2_quies =
      params->have_parameter<double>(CP_name+"/quiescence/O2_threshold")
      ? params->get<double>(CP_name+"/quiescence/O2_threshold") : -1.0;
    const double nut_quies =
      params->have_parameter<double>(CP_name+"/quiescence/nutrient_threshold")
      ? params->get<double>(CP_name+"/quiescence/nutrient_threshold") : -1.0;
    if ((crowd_entry >= 0.0 && env.local_crowding >= crowd_entry) ||
        (o2_quies   >= 0.0 && env.local_O2       <  o2_quies)    ||
        (nut_quies  >= 0.0 && env.local_nutrient  <  nut_quies))
      st.should_enter_G0 = true;
  }
  //
  // ================================================================
  // Apoptosis and necrosis recommendations (hazard-rate execution is
  // in Mechanism 11 Run(); here we just set the flags)
  // ================================================================
  // Anoikis: low ECM adhesion → BIM/BAD activation → intrinsic apoptosis
  //   Modulated by anoikis_resistance (e.g. via BCL-2/BCL-xL overexpression)
  if (this->GetCanApoptose()) {
    const double anoikis_thr =
      params->have_parameter<double>(CP_name+"/ecm/anoikis_threshold")
      ? params->get<double>(CP_name+"/ecm/anoikis_threshold") : -1.0;
    if (anoikis_thr >= 0.0 && env.ecm_adhesion < anoikis_thr)
      st.should_enter_Ap = true;
  }
  // Necrosis: severe O2 deprivation → ATP collapse → membrane failure
  {
    const double o2_nec =
      params->have_parameter<double>(CP_name+"/can_necrose/O2_threshold")
      ? params->get<double>(CP_name+"/can_necrose/O2_threshold") : -1.0;
    if (o2_nec >= 0.0 && env.local_O2 <= o2_nec) st.should_enter_Nec = true;
  }
  //
  // ================================================================
  // DNA repair gating
  // Repair requires O2 (NHEJ/HR need ATP) and nutrients (repair synthesis).
  // HIF-1α under hypoxia competes for repair factor binding sites.
  // ================================================================
  {
    const double o2_rep =
      params->have_parameter<double>(CP_name+"/repair/O2_min")
      ? params->get<double>(CP_name+"/repair/O2_min") : 0.0;
    const double nut_rep =
      params->have_parameter<double>(CP_name+"/repair/nutrient_min")
      ? params->get<double>(CP_name+"/repair/nutrient_min") : 0.0;
    if (env.local_O2 < o2_rep || env.local_nutrient < nut_rep)
      st.repair_allowed = false;
  }
  // ================================================================
  // GRN output fusion (optional)
  // ================================================================
  if (regulatory_backend_id_ != 0)
    {
      st.can_enter_S = st.can_enter_S && (regulatory_output_.can_enter_S > 0.5);
      st.can_enter_M = st.can_enter_M && (regulatory_output_.can_enter_M > 0.5);
      st.can_grow = st.can_grow && (regulatory_output_.proliferation_signal > 0.25);
      if (regulatory_output_.quiescence_hazard > 0.05)
        st.should_enter_G0 = true;
      if (this->GetCanApoptose() && regulatory_output_.apoptosis_hazard > 0.2)
        st.should_enter_Ap = true;
      if (regulatory_output_.necrosis_hazard > 0.2)
        st.should_enter_Nec = true;
      st.repair_allowed = st.repair_allowed && (regulatory_output_.repair_capacity > 0.1);
    }
  return st;
}
// -----------------------------------------------------------------------------
inline
double bdm::BiologicalCell::ComputeGrowthModulation(
    const bdm::BiologicalCell::MicroenvironmentState& env) const
{
  // Returns a multiplicative growth-rate modifier f_env ∈ [0,1].
  // If no modulation parameters are defined, returns 1.0 (backward-compatible).
  // Each factor represents a distinct signalling pathway; all are independent
  // and the composite modulation is the product of individual factors.
  if (!this->GetPhenotype()) return 0.0;
  auto* params = this->params();
  const std::string& CP_name =
    params->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  double f = 1.0;
  //
  // O2 modulation: mTOR/HIF-1α — hypoxia reduces mTORC1 activity and
  // protein synthesis. f_O2 = O2 / O2_reference, capped at 1.
  if (params->have_parameter<double>(CP_name+"/can_grow/O2_reference")) {
    const double o2_ref = params->get<double>(CP_name+"/can_grow/O2_reference");
    if (o2_ref > 0.0)
      f *= std::min(1.0, std::max(0.0, env.local_O2 / o2_ref));
  }
  //
  // Nutrient modulation: PI3K/Akt/mTOR nutrient sensing (amino acids, glucose).
  // f_nut = nutrient / nut_reference, capped at 1.
  if (params->have_parameter<double>(CP_name+"/can_grow/nutrient_reference")) {
    const double nut_ref = params->get<double>(CP_name+"/can_grow/nutrient_reference");
    if (nut_ref > 0.0)
      f *= std::min(1.0, std::max(0.0, env.local_nutrient / nut_ref));
  }
  //
  // ECM stiffness modulation: YAP/TAZ mechanosensing.
  // Stiff ECM → nuclear YAP/TAZ → TEAD → cyclin D1 / CTGF → proliferative.
  // f_stiff = stiffness / stiffness_reference, capped at 1.
  if (params->have_parameter<double>(CP_name+"/can_grow/ecm_stiffness_reference")) {
    const double stiff_ref = params->get<double>(CP_name+"/can_grow/ecm_stiffness_reference");
    if (stiff_ref > 0.0)
      f *= std::min(1.0, std::max(0.0, env.ecm_stiffness / stiff_ref));
  }
  //
  // Crowding modulation: contact inhibition of growth.
  // E-cadherin → LATS1/2 → YAP/TAZ cytoplasmic sequestration.
  // f_crowd = max(0, 1 - crowding / crowding_max_grow)
  if (params->have_parameter<double>(CP_name+"/can_grow/crowding_max")) {
    const double crowd_max = params->get<double>(CP_name+"/can_grow/crowding_max");
    if (crowd_max > 0.0)
      f *= std::max(0.0, 1.0 - env.local_crowding / crowd_max);
  }
  //
  // ROS/stress modulation: high ROS activates p38/JNK → growth arrest.
  // f_rons = max(0, 1 - rons / rons_max_grow)
  if (params->have_parameter<double>(CP_name+"/can_grow/rons_max")) {
    const double rons_max = params->get<double>(CP_name+"/can_grow/rons_max");
    if (rons_max > 0.0)
      f *= std::max(0.0, 1.0 - env.local_rons / rons_max);
  }
  return std::max(0.0, std::min(1.0, f));
}
// =============================================================================
// MECHANISM 12 — CAP/PAM TREATMENT RESPONSE METHODS
// =============================================================================
// UpdateCAPIntracellular(): RONS-specific intracellular dynamics update.
// Must be called AFTER RunIntracellular() so ros_internal_, dna_damage_,
// and DDR pathway variables (atm_active_, chk1_active_, p53_active_, …)
// are already current for this timestep.
//
// Biological model overview:
//   CAP/PAM → extracellular H2O2 + NO2- (RONS cocktail)
//   → intracellular H2O2 (via AQP3/AQP8) → Fenton → OH· → 8-oxoG, SSB/DSB
//   → intracellular NO2-/ONOO- → nitrosative base lesions
//   → ATM (DSB sensor) + ATR (ssDNA/fork-stall sensor) → CHK1/CHK2 phosphorylation
//   → γH2AX foci form at DSB sites (ATM/ATR → H2AX-Ser139)
//   → CHK1 phosphorylation reported in EGI-1 and HuCCT1 after CAP/PAM
//   → p53 phosphorylation/accumulation → p21 (CDKN1A) → CDK inhibition
//   → G1/S and G2/M cell-cycle arrest (primary response before apoptosis)
//   → after sufficient damage and failed repair: BAX/BCL-2 imbalance →
//     cytochrome C → Apaf-1 → caspase-9 → caspase-3 (cleaved)
//   → PARP-1 cleavage (caspase-3 substrate) seals apoptosis fate
//   → Annexin V+ / 7-AAD+: early/late apoptosis detected by flow cytometry
//   → EGI-1 apoptosis ~72h; HuCCT1 ~48h (cell-line encoded as parameters, not names)
//   → In vivo: 8-oxoG IHC + cleaved caspase-3 IHC positive in CAP-treated xenografts
//   → Primary hepatocytes: weaker CHK1/p53 activation, no PARP cleavage observed
//     (represented as phenotype-specific parameter set with low sensitivity, not hardcoded)
// =============================================================================
inline
void bdm::BiologicalCell::UpdateCAPIntracellular()
{
  // Only viable cells with non-zero phenotype maintain CAP dynamics
  if (!this->GetPhenotype()) return;
  // Access resource manager and parameters
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  const double dt = this->params()->get<double>("time_step");
  const std::string& CP_name =
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const bdm::Double3 xyz = this->GetPosition();
  const auto& substances = this->params()->get<std::vector<std::string>>("substances");
  // ================================================================
  // STEP A — Sample extracellular RONS fields
  // ================================================================
  // Total extracellular RONS (H2O2 + NO2_); species-resolved where available.
  // When RONS fields are absent (e.g. untreated controls), ext_rons_total = 0.
  double ext_h2o2 = 0.0, ext_no2 = 0.0;
  if (std::find(substances.begin(), substances.end(), "H2O2") != substances.end())
    if (auto* dg = rm->GetDiffusionGrid("H2O2"))
      ext_h2o2 = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
  if (std::find(substances.begin(), substances.end(), "NO2_") != substances.end())
    if (auto* dg = rm->GetDiffusionGrid("NO2_"))
      ext_no2 = std::max(0.0, GetInterpolatedValue(dg, xyz, this->params()));
  const double ext_rons_total = ext_h2o2 + ext_no2;
  // ================================================================
  // STEP B — Update cap_dose_integral_ and time_since_cap_
  //   cap_dose_integral_ = ∫ ext_rons_total dt (area under curve)
  //   time_since_cap_ counts steps since first non-zero RONS exposure.
  // ================================================================
  cap_dose_integral_ += dt * ext_rons_total;
  if (ext_rons_total > 0.0)
    ++time_since_cap_;
  // ================================================================
  // STEP C — Update rns_internal_ (reactive nitrogen species)
  //   Tracks NO2--derived intracellular RNS separately from ROS.
  //   NO2- enters via ion channels; reacts with O2-· to form ONOO-.
  //   Separate ODE from ros_internal_ (RunIntracellular already handles H2O2 → ROS).
  //
  //   ODE: drns/dt = k_rns_uptake * ext_no2 * membrane_permeability
  //                - k_rns_clear * antioxidant_capacity * rns
  //   Use exponential integrator for stability.
  // ================================================================
  {
    const double k_rns_uptake =
      this->params()->have_parameter<double>(CP_name+"/cap/rns/k_uptake")
      ? this->params()->get<double>(CP_name+"/cap/rns/k_uptake") : 0.0;
    const double k_rns_clear =
      this->params()->have_parameter<double>(CP_name+"/cap/rns/k_clearance")
      ? this->params()->get<double>(CP_name+"/cap/rns/k_clearance") : 0.0;
    // Read membrane_permeability from parameter if configured;
    // otherwise use the cell's current membrane_permeability_ state.
    const double mp =
      this->params()->have_parameter<double>(CP_name+"/cap/membrane_permeability")
      ? this->params()->get<double>(CP_name+"/cap/membrane_permeability")
      : membrane_permeability_;
    const double prod_rns = k_rns_uptake * ext_no2 * mp;
    const double decay_rns = k_rns_clear * antioxidant_capacity_;
    if (decay_rns > 1.0e-12) {
      const double ss_rns = prod_rns / decay_rns;
      rns_internal_ = ss_rns + (rns_internal_ - ss_rns) * std::exp(-decay_rns * dt);
    } else {
      rns_internal_ += dt * prod_rns;
    }
    if (rns_internal_ < 0.0) rns_internal_ = 0.0;
  }
  // ================================================================
  // STEP D — Update membrane_permeability_
  //   High ROS causes lipid peroxidation → reduced membrane integrity.
  //   Conversely, moderate peroxidation may initially increase permeability.
  //   Simple exponential relaxation toward baseline, perturbed by ros_internal_.
  //   Only active when cap/membrane_permeability/k_damage is defined.
  // ================================================================
  {
    const double k_mp_damage =
      this->params()->have_parameter<double>(CP_name+"/cap/membrane_permeability/k_damage")
      ? this->params()->get<double>(CP_name+"/cap/membrane_permeability/k_damage") : 0.0;
    const double k_mp_repair =
      this->params()->have_parameter<double>(CP_name+"/cap/membrane_permeability/k_repair")
      ? this->params()->get<double>(CP_name+"/cap/membrane_permeability/k_repair") : 0.0;
    if (k_mp_damage > 0.0 || k_mp_repair > 0.0) {
      // membrane_permeability_ decreases with sustained high ROS
      const double ros_ref =
        this->params()->have_parameter<double>(CP_name+"/cap/membrane_permeability/ros_reference")
        ? this->params()->get<double>(CP_name+"/cap/membrane_permeability/ros_reference") : 1.0;
      const double stress = std::min(1.0, std::max(0.0, ros_internal_ / std::max(1.0e-9, ros_ref)));
      membrane_permeability_ += dt * (-k_mp_damage * stress * membrane_permeability_
                                      + k_mp_repair * (1.0 - membrane_permeability_));
      membrane_permeability_ = std::clamp(membrane_permeability_, 0.0, 1.0);
    }
  }
  // ================================================================
  // STEP E — Update repair_capacity_
  //   Repair enzymes (OGG1/APE1 for BER, Ku70/Ku80 for NHEJ, RAD51 for HR)
  //   can be inactivated by sustained oxidative stress.
  //   Simple logistic suppression by ros_internal_ / rns_internal_.
  //   Only active when cap/repair_capacity parameters are configured.
  // ================================================================
  {
    const double k_rc_suppress =
      this->params()->have_parameter<double>(CP_name+"/cap/repair_capacity/k_suppress")
      ? this->params()->get<double>(CP_name+"/cap/repair_capacity/k_suppress") : 0.0;
    const double k_rc_recover =
      this->params()->have_parameter<double>(CP_name+"/cap/repair_capacity/k_recovery")
      ? this->params()->get<double>(CP_name+"/cap/repair_capacity/k_recovery") : 0.0;
    const double rc_baseline =
      this->params()->have_parameter<double>(CP_name+"/cap/repair_capacity/baseline")
      ? this->params()->get<double>(CP_name+"/cap/repair_capacity/baseline") : 1.0;
    if (k_rc_suppress > 0.0 || k_rc_recover > 0.0) {
      const double rons_total_int = ros_internal_ + rns_internal_;
      const double ros_ref =
        this->params()->have_parameter<double>(CP_name+"/cap/repair_capacity/ros_reference")
        ? this->params()->get<double>(CP_name+"/cap/repair_capacity/ros_reference") : 1.0;
      const double stress = std::min(1.0, std::max(0.0, rons_total_int / std::max(1.0e-9, ros_ref)));
      repair_capacity_ += dt * (-k_rc_suppress * stress * repair_capacity_
                                + k_rc_recover * (rc_baseline - repair_capacity_));
      repair_capacity_ = std::clamp(repair_capacity_, 0.0, 1.0);
    }
  }
  // ================================================================
  // STEP F — Compute CAP-specific marker proxies
  //   These are normalized observables corresponding to experimental readouts.
  //   They are computed from ros_internal_, dna_damage_, and DDR pathway state.
  // ================================================================
  // F1 — 8-oxoguanine proxy (8-oxoG):
  //   Oxidative DNA base lesion; formed when OH· (from H2O2 Fenton reaction) attacks
  //   guanine N7/C8. Accumulates proportional to cumulative oxidative stress.
  //   Repaired by BER (OGG1 → APE1 → pol-β → ligase).
  //   Simple linear combination: 8-oxoG ≈ α_ros * ros + α_dam * dna_damage.
  {
    const double alpha_ros =
      this->params()->have_parameter<double>(CP_name+"/cap/markers/8oxoG/alpha_ros")
      ? this->params()->get<double>(CP_name+"/cap/markers/8oxoG/alpha_ros") : 0.5;
    const double alpha_dam =
      this->params()->have_parameter<double>(CP_name+"/cap/markers/8oxoG/alpha_damage")
      ? this->params()->get<double>(CP_name+"/cap/markers/8oxoG/alpha_damage") : 0.5;
    oxidative_damage_8oxoG_proxy_ =
      std::clamp(alpha_ros * ros_internal_ + alpha_dam * dna_damage_, 0.0, 1.0);
  }
  // F2 — γH2AX proxy (DSB marker):
  //   H2AX phosphorylation at Ser139 by ATM/ATR at DSB sites (γH2AX foci).
  //   Rapid response: γH2AX forms within minutes of DSB induction.
  //   Reported positive in CAP/PAM-treated EGI-1 and HuCCT1 (immunofluorescence).
  //   Sigmoid function of dna_damage_ captures threshold-like foci formation.
  //   γH2AX_proxy = sigmoid(dna_damage / half_sat)
  {
    const double half_sat =
      this->params()->have_parameter<double>(CP_name+"/cap/markers/gammaH2AX/half_saturation")
      ? this->params()->get<double>(CP_name+"/cap/markers/gammaH2AX/half_saturation") : 0.3;
    const double steepness =
      this->params()->have_parameter<double>(CP_name+"/cap/markers/gammaH2AX/steepness")
      ? this->params()->get<double>(CP_name+"/cap/markers/gammaH2AX/steepness") : 10.0;
    if (half_sat > 1.0e-9)
      dsb_damage_gammaH2AX_proxy_ =
        1.0 / (1.0 + std::exp(-steepness * (dna_damage_ - half_sat)));
    else
      dsb_damage_gammaH2AX_proxy_ = 1.0;
  }
  // F3 — pCHK1 proxy: already stored as chk1_active_ (computed in UpdateDdrPathway)
  // F4 — p53 proxy: already stored as p53_active_ (computed in UpdateDdrPathway)
  //   CHK1 phosphorylation and p53 phosphorylation/accumulation reported in
  //   CAP/PAM-treated CCA cells (CHK1-Ser345; p53-Ser15 are canonical markers).
  //   Primary hepatocytes show weaker CHK1/p53 activation under comparable conditions.
  // F5 — PARP cleavage proxy:
  //   PARP-1 is cleaved by activated caspase-3 → 89 kDa fragment.
  //   Cleavage only occurs AFTER caspase-3 is activated (apoptosis execution phase).
  //   Therefore, parp_cleavage_proxy_ only increases after caspase3_activation_proxy_.
  //   Reported as western blot marker in CAP-treated CCA cells.
  {
    const double k_parp_cleavage =
      this->params()->have_parameter<double>(CP_name+"/cap/markers/PARP/k_cleavage")
      ? this->params()->get<double>(CP_name+"/cap/markers/PARP/k_cleavage") : 2.0;
    // PARP cleavage driven by caspase-3 activation (execution phase only)
    const double drive = k_parp_cleavage * caspase3_activation_proxy_ * (1.0 - parp_cleavage_proxy_);
    parp_cleavage_proxy_ = std::clamp(parp_cleavage_proxy_ + dt * drive, 0.0, 1.0);
  }
  // F6 — Cleaved caspase-3 proxy:
  //   Executioner caspase; activated by intrinsic (mitochondrial) pathway:
  //   BAX/BCL-2 imbalance → cytochrome C release → Apaf-1/caspase-9 apoptosome
  //   → caspase-9 cleaves and activates caspase-3.
  //   Increases ONLY when apoptosis_commitment_state_ exceeds commitment threshold.
  //   Reported by IHC in CAP-treated CCA xenografts.
  {
    const double commit_thr =
      this->params()->have_parameter<double>(CP_name+"/cap/apoptosis/commitment_threshold")
      ? this->params()->get<double>(CP_name+"/cap/apoptosis/commitment_threshold") : 0.7;
    const double k_casp3 =
      this->params()->have_parameter<double>(CP_name+"/cap/markers/caspase3/k_activation")
      ? this->params()->get<double>(CP_name+"/cap/markers/caspase3/k_activation") : 1.5;
    if (apoptosis_commitment_state_ >= commit_thr) {
      const double drive = k_casp3 * apoptosis_commitment_state_ * (1.0 - caspase3_activation_proxy_);
      caspase3_activation_proxy_ = std::clamp(caspase3_activation_proxy_ + dt * drive, 0.0, 1.0);
    }
  }
  // ================================================================
  // STEP G — Update apoptosis_commitment_state_
  //   Represents progressive accumulation of pro-apoptotic signals:
  //   p53/BAX accumulation, BCL-2 degradation, mitochondrial membrane
  //   potential collapse, cytochrome C release.
  //
  //   commitment increases when:
  //     (a) checkpoint_activation (max of pCHK1, p53) > apoptosis_checkpoint_threshold
  //         AND dna_damage > apoptosis_damage_threshold
  //     (b) arrest_time > max_repair_time AND dna_damage still elevated
  //         (irreparable damage → forced commitment)
  //
  //   commitment decreases slowly when damage/stress falls below threshold
  //   (repair window: cell can recover if commitment is still low)
  //
  //   This models:
  //   - The temporal delay between DNA damage detection and apoptosis onset
  //   - The repair window: cells with low damage/commitment can recover
  //   - The commitment threshold: once crossed, apoptosis is irreversible
  // ================================================================
  {
    const double checkpoint_act = std::max(chk1_active_, p53_active_);
    const double damage_for_commit =
      this->params()->have_parameter<double>(CP_name+"/cap/apoptosis/damage_threshold")
      ? this->params()->get<double>(CP_name+"/cap/apoptosis/damage_threshold")
      : (this->params()->have_parameter<double>(CP_name+"/intracellular/damage/threshold")
         ? this->params()->get<double>(CP_name+"/intracellular/damage/threshold") * 0.6
         : 0.6);
    const double chk_for_commit =
      this->params()->have_parameter<double>(CP_name+"/cap/apoptosis/checkpoint_threshold")
      ? this->params()->get<double>(CP_name+"/cap/apoptosis/checkpoint_threshold") : 0.5;
    const double k_commit =
      this->params()->have_parameter<double>(CP_name+"/cap/apoptosis/k_commitment")
      ? this->params()->get<double>(CP_name+"/cap/apoptosis/k_commitment") : 0.1;
    const double k_recover =
      this->params()->have_parameter<double>(CP_name+"/cap/apoptosis/k_recovery")
      ? this->params()->get<double>(CP_name+"/cap/apoptosis/k_recovery") : 0.05;
    const int max_repair_time =
      this->params()->have_parameter<int>(CP_name+"/cap/repair/max_repair_time")
      ? this->params()->get<int>(CP_name+"/cap/repair/max_repair_time")
      : (this->params()->have_parameter<int>(CP_name+"/phase_dwell/max_arrest_time")
         ? this->params()->get<int>(CP_name+"/phase_dwell/max_arrest_time") : 999999);
    // Is the cell in a pro-apoptotic state?
    const bool damage_driven =
      (dna_damage_ > damage_for_commit) && (checkpoint_act > chk_for_commit);
    const bool arrest_driven =
      (arrest_time_ > max_repair_time) && (dna_damage_ > damage_for_commit * 0.5);
    if (damage_driven || arrest_driven) {
      // Commitment accumulation: sigmoid-shaped drive proportional to damage excess
      const double excess_damage = std::max(0.0, dna_damage_ - damage_for_commit);
      const double drive = k_commit * (excess_damage + 0.5 * checkpoint_act)
                           * (1.0 - apoptosis_commitment_state_);
      apoptosis_commitment_state_ = std::clamp(apoptosis_commitment_state_ + dt * drive, 0.0, 1.0);
    } else if (apoptosis_commitment_state_ > 0.0 && !damage_driven && !arrest_driven) {
      // Recovery: commitment decays when conditions improve (repair window still open)
      const double commit_thr =
        this->params()->have_parameter<double>(CP_name+"/cap/apoptosis/commitment_threshold")
        ? this->params()->get<double>(CP_name+"/cap/apoptosis/commitment_threshold") : 0.7;
      if (apoptosis_commitment_state_ < commit_thr) {
        // Below commitment threshold: full recovery possible
        apoptosis_commitment_state_ = std::clamp(
          apoptosis_commitment_state_ - dt * k_recover, 0.0, 1.0);
      }
      // Above commitment threshold: no recovery (apoptosis irreversible)
    }
  }
}
// =============================================================================
// EvaluateCAPCheckpointState(): Mechanism 12 checkpoint controller.
// Returns CAPCheckpointState encoding all gate decisions for CAP/PAM-treated cells.
//
// Uses the pre-sampled MicroenvironmentState (O2, nutrient, ECM, crowding)
// combined with the intracellular RONS/damage state (ros_internal_, dna_damage_,
// DDR pathway variables, apoptosis_commitment_state_) to make checkpoint decisions.
//
// Biology:
//   G1/S block: p53/p21-CDK2/CyclinE axis (block DNA replication with damage)
//   Intra-S slowing: ATR–CHK1–CDC25A ubiquitination (stalled replication)
//   G2/M block: CHK1/CHK2–CDC25C–CDK1/CyclinB (primary CAP arrest gate)
//   Division block: never with unresolved damage > G2M threshold
//   Repair window: O2/nutrient must be sufficient for NER/NHEJ/HR
//   Apoptosis: commitment_state > threshold OR arrest > max_repair_time + damage
//   Necrosis: extreme ROS/energy failure (high ros + hypoxia + membrane damage)
// =============================================================================
inline
bdm::BiologicalCell::CAPCheckpointState
bdm::BiologicalCell::EvaluateCAPCheckpointState(
    const bdm::BiologicalCell::MicroenvironmentState& env)
{
  CAPCheckpointState st; // all fields default-initialized to permissive/safe
  if (!this->GetPhenotype()) return st;
  auto* params = this->params();
  const std::string& CP_name =
    params->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  const double dt = params->get<double>("time_step");
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  // ================================================================
  // G1/S checkpoint
  //   Primary gate: p53/p21-mediated CDK2/CyclinE inhibition.
  //   The molecular DDR pathway (UpdateDdrPathway) computes p21_level_ and
  //   cdk_activity_ from dna_damage_. IsMolecularG1SCheckpointBlocked() uses these.
  //   Additional fallback: legacy aggregate damage threshold.
  //
  //   CAP context: CAP-induced p53 phosphorylation and CHK1 activation block
  //   G1/S transition. Arrest appears before apoptosis (DDR response first).
  // ================================================================
  if (bdm::IsMolecularG1SCheckpointBlocked(this))
    st.can_enter_S = false;
  // Fallback: legacy aggregate damage threshold (when DDR not fully configured)
  if (st.can_enter_S &&
      !bdm::IsDdrPathwayEnabled(*params, CP_name) &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/damage_threshold") &&
      dna_damage_ > params->get<double>(CP_name+"/checkpoint/G1S/damage_threshold"))
    st.can_enter_S = false;
  // CAP-specific G1/S damage threshold (separate from legacy, if defined)
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/cap/checkpoint/G1S/damage_threshold") &&
      dna_damage_ > params->get<double>(CP_name+"/cap/checkpoint/G1S/damage_threshold"))
    st.can_enter_S = false;
  // Hypoxia gate: HIF-1α suppresses CDK2/CyclinD → G1 arrest under severe hypoxia
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/O2_threshold") &&
      env.local_O2 < params->get<double>(CP_name+"/checkpoint/G1S/O2_threshold"))
    st.can_enter_S = false;
  // Nutrient gate: AMPK → p27/p21 under nutrient stress
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/nutrient_min") &&
      env.local_nutrient < params->get<double>(CP_name+"/checkpoint/G1S/nutrient_min"))
    st.can_enter_S = false;
  // Crowding gate: contact inhibition
  if (st.can_enter_S &&
      params->have_parameter<double>(CP_name+"/checkpoint/G1S/crowding_threshold") &&
      env.local_crowding >= params->get<double>(CP_name+"/checkpoint/G1S/crowding_threshold"))
    st.can_enter_S = false;
  // ================================================================
  // G2/M checkpoint — THE CRITICAL CAP ARREST GATE
  //   CHK1 (ATR-driven) and CHK2 (ATM-driven) phosphorylate CDC25C → 14-3-3σ
  //   sequestration → CDC25C inactive → CDK1/CyclinB cannot be dephosphorylated
  //   → CDK1/CyclinB complex stays inactive → mitosis entry BLOCKED.
  //   Division with unresolved DSBs = mitotic catastrophe → lethal.
  //   γH2AX-positive cells must NOT divide: enforced here absolutely.
  //
  //   CAP context: most prominently reported arrest checkpoint in CCA cells.
  //   γH2AX-positive cells arrested; arrest duration correlates with dose.
  // ================================================================
  if (bdm::IsMolecularG2MCheckpointBlocked(this))
    st.can_enter_M = false;
  // Fallback: legacy damage threshold
  if (st.can_enter_M &&
      !bdm::IsDdrPathwayEnabled(*params, CP_name) &&
      params->have_parameter<double>(CP_name+"/checkpoint/G2M/damage_threshold") &&
      dna_damage_ > params->get<double>(CP_name+"/checkpoint/G2M/damage_threshold"))
    st.can_enter_M = false;
  // CAP-specific G2/M damage threshold
  if (st.can_enter_M &&
      params->have_parameter<double>(CP_name+"/cap/checkpoint/G2M/damage_threshold") &&
      dna_damage_ > params->get<double>(CP_name+"/cap/checkpoint/G2M/damage_threshold"))
    st.can_enter_M = false;
  // γH2AX absolute block: cells with significant DSB marker must never enter mitosis
  {
    const double gh2ax_div_thr =
      params->have_parameter<double>(CP_name+"/cap/checkpoint/G2M/gammaH2AX_block_threshold")
      ? params->get<double>(CP_name+"/cap/checkpoint/G2M/gammaH2AX_block_threshold") : 0.5;
    if (dsb_damage_gammaH2AX_proxy_ > gh2ax_div_thr)
      st.can_enter_M = false;
  }
  // Hypoxia gate
  if (st.can_enter_M &&
      params->have_parameter<double>(CP_name+"/checkpoint/G2M/O2_threshold") &&
      env.local_O2 < params->get<double>(CP_name+"/checkpoint/G2M/O2_threshold"))
    st.can_enter_M = false;
  // If G2/M is blocked, division is blocked
  if (!st.can_enter_M) { /* can_divide is already constrained by can_enter_M in Mech 12 */ }
  // ================================================================
  // Intra-S checkpoint
  //   ATR–CHK1–CDC25A ubiquitination during S-phase.
  //   Moderate damage during replication → stalled forks → ATR activation.
  // ================================================================
  if (phase_ == Phase::Sy) {
    const double intraS_thr =
      params->have_parameter<double>(CP_name+"/cap/checkpoint/IntraS/damage_threshold")
      ? params->get<double>(CP_name+"/cap/checkpoint/IntraS/damage_threshold")
      : (params->have_parameter<double>(CP_name+"/checkpoint/IntraS/damage_threshold")
         ? params->get<double>(CP_name+"/checkpoint/IntraS/damage_threshold") : 1.0e99);
    const double intraS_ros_thr =
      params->have_parameter<double>(CP_name+"/cap/checkpoint/IntraS/ros_threshold")
      ? params->get<double>(CP_name+"/cap/checkpoint/IntraS/ros_threshold")
      : (params->have_parameter<double>(CP_name+"/checkpoint/IntraS/ros_threshold")
         ? params->get<double>(CP_name+"/checkpoint/IntraS/ros_threshold") : 1.0e99);
    if (dna_damage_ > intraS_thr || ros_internal_ > intraS_ros_thr)
      st.intra_s_blocked = true;
  }
  // ================================================================
  // General arrest flag
  // ================================================================
  st.must_arrest = (!st.can_enter_S || !st.can_enter_M || st.intra_s_blocked);
  // ================================================================
  // Repair allowance
  //   DNA repair (NHEJ/HR/BER) requires O2 (for NADPH, ATP) and metabolic energy.
  //   Under severe hypoxia, HIF-1α competes for chromatin remodelling and
  //   repair-factor recruitment (RAD51/OGG1/XRCC1 are O2-sensitive).
  // ================================================================
  {
    const double o2_rep =
      params->have_parameter<double>(CP_name+"/cap/repair/O2_min")
      ? params->get<double>(CP_name+"/cap/repair/O2_min")
      : (params->have_parameter<double>(CP_name+"/repair/O2_min")
         ? params->get<double>(CP_name+"/repair/O2_min") : 0.0);
    const double nut_rep =
      params->have_parameter<double>(CP_name+"/cap/repair/nutrient_min")
      ? params->get<double>(CP_name+"/cap/repair/nutrient_min")
      : (params->have_parameter<double>(CP_name+"/repair/nutrient_min")
         ? params->get<double>(CP_name+"/repair/nutrient_min") : 0.0);
    if (env.local_O2 < o2_rep || env.local_nutrient < nut_rep)
      st.repair_allowed = false;
  }
  // ================================================================
  // Apoptosis commitment flag
  //   must_enter_Ap = true when:
  //   (a) apoptosis_commitment_state_ >= commitment_threshold (slow accumulation)
  //   (b) arrest_time_ > max_repair_time AND damage remains elevated
  //   (c) dna_damage_ > hard_apoptosis_threshold (extreme damage, fast path)
  //
  //   The stochastic apoptosis hazard is applied in the Mechanism 12 Run()
  //   (rate-to-probability: P = 1 - exp(-k * dt)).
  // ================================================================
  {
    const double commit_thr =
      params->have_parameter<double>(CP_name+"/cap/apoptosis/commitment_threshold")
      ? params->get<double>(CP_name+"/cap/apoptosis/commitment_threshold") : 0.7;
    if (apoptosis_commitment_state_ >= commit_thr)
      st.must_enter_Ap = true;
    // Prolonged arrest path: max_repair_time exceeded AND damage still high
    const int max_repair_time =
      params->have_parameter<int>(CP_name+"/cap/repair/max_repair_time")
      ? params->get<int>(CP_name+"/cap/repair/max_repair_time")
      : (params->have_parameter<int>(CP_name+"/phase_dwell/max_arrest_time")
         ? params->get<int>(CP_name+"/phase_dwell/max_arrest_time") : 999999);
    const double arrest_ap_damage_thr =
      params->have_parameter<double>(CP_name+"/cap/apoptosis/damage_threshold")
      ? params->get<double>(CP_name+"/cap/apoptosis/damage_threshold")
      : (params->have_parameter<double>(CP_name+"/intracellular/damage/threshold")
         ? params->get<double>(CP_name+"/intracellular/damage/threshold") * 0.4 : 0.4);
    if (arrest_time_ > max_repair_time && dna_damage_ > arrest_ap_damage_thr)
      st.must_enter_Ap = true;
    // Hard damage threshold (extreme unrepaired damage → immediate apoptosis path)
    const double hard_ap_thr =
      params->have_parameter<double>(CP_name+"/cap/apoptosis/hard_damage_threshold")
      ? params->get<double>(CP_name+"/cap/apoptosis/hard_damage_threshold")
      : (params->have_parameter<double>(CP_name+"/intracellular/damage/threshold")
         ? params->get<double>(CP_name+"/intracellular/damage/threshold") : 1.0e99);
    if (dna_damage_ > hard_ap_thr)
      st.must_enter_Ap = true;
    // cap_dose_exceeded: track chronic exposure signal (not necessarily apoptotic alone)
    const double dose_thr =
      params->have_parameter<double>(CP_name+"/cap/apoptosis/dose_integral_threshold")
      ? params->get<double>(CP_name+"/cap/apoptosis/dose_integral_threshold") : 1.0e99;
    if (cap_dose_integral_ > dose_thr)
      st.cap_dose_exceeded = true;
  }
  // ================================================================
  // Necrosis flag
  //   Extreme RONS + hypoxia + membrane damage → energy collapse.
  //   Biologically: peroxynitrite (ONOO-) destroys mitochondrial membranes;
  //   ATP depletion triggers oncosis → HMGB1/DAMPs release.
  //   Separate from apoptosis (necrosis is passive/uncontrolled).
  // ================================================================
  if (this->GetCanApoptose()) {
    // ROS necrosis: extreme intracellular ROS or RNS
    const double ros_nec_thr =
      params->have_parameter<double>(CP_name+"/can_necrose/ros_threshold")
      ? params->get<double>(CP_name+"/can_necrose/ros_threshold") : 1.0e99;
    if (ros_internal_ > ros_nec_thr || rns_internal_ > ros_nec_thr * 0.8)
      st.must_enter_Nec = true;
    // Severe hypoxia necrosis
    const double o2_nec_thr =
      params->have_parameter<double>(CP_name+"/can_necrose/O2_threshold")
      ? params->get<double>(CP_name+"/can_necrose/O2_threshold") : -1.0;
    if (o2_nec_thr >= 0.0 && env.local_O2 <= o2_nec_thr)
      st.must_enter_Nec = true;
  }
  // ================================================================
  // GRN output fusion (optional)
  // ================================================================
  if (regulatory_backend_id_ != 0)
    {
      st.can_enter_S = st.can_enter_S && (regulatory_output_.can_enter_S > 0.5);
      st.can_enter_M = st.can_enter_M && (regulatory_output_.can_enter_M > 0.5);
      st.must_arrest = st.must_arrest || (regulatory_output_.can_enter_S <= 0.5)
                      || (regulatory_output_.can_enter_M <= 0.5);
      st.repair_allowed = st.repair_allowed && (regulatory_output_.repair_capacity > 0.1);
      if (this->GetCanApoptose() && regulatory_output_.apoptosis_hazard > 0.2)
        st.must_enter_Ap = true;
      if (regulatory_output_.necrosis_hazard > 0.2)
        st.must_enter_Nec = true;
    }
  return st;
}
// =============================================================================
#endif // _BIOLOGICAL_CELL_INLINE_H_
// =============================================================================
