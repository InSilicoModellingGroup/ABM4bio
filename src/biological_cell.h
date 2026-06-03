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
#ifndef _BIOLOGICAL_CELL_H_
#define _BIOLOGICAL_CELL_H_
// =============================================================================
#include "./global.h"
#include "./biology.h"
#include "./biochemical.h"
#include "./cell_protrusion.h"
#include "./obstacles.h"
#include "./io_flux.h"
#include "./regulatory/regulatory_model-inline.h"
namespace bdm { void UpdateDdrPathway(BiologicalCell* cell); }
// =============================================================================
namespace bdm {
// =============================================================================
class BiologicalCell : public bdm::neuroscience::NeuronSoma {
BDM_AGENT_HEADER(BiologicalCell, bdm::neuroscience::NeuronSoma, 1);
//
public:
  // local enumerator that monitors the phase of a cell's circle
  enum Phase {
    Ap =-1,
    I0 =0, G1 =1, Sy =2, G2 =3, Di =4, Tr =5
  };
  // -----------------------------------------------------------------------
  // Mechanism 11: Snapshot of local microenvironment sampled once per step.
  // Sampled by SampleMicroenvironment() to avoid redundant grid lookups.
  // All values use safe defaults so old input files still run without change.
  struct MicroenvironmentState {
    double local_O2       = 1.0; ///< O2 concentration (normalized, [0,1]+)
    double local_nutrient = 1.0; ///< nutrient/glucose level (normalized)
    double local_crowding = 0.0; ///< local cell-volume occupancy ratio [0,inf)
    double ecm_density    = 1.0; ///< ECM (collagen/fibronectin) field value
    // ECM_stiffness field ("ECM_stiffness" substance if present):
    // represents effective mechanosensing input to YAP/TAZ / ROCK signalling.
    double ecm_stiffness  = 1.0; ///< effective ECM stiffness modifier
    // ECM_adhesion_ligand_density: derived from ecm_density * integrin_sensitivity
    // or from an "ECM_adhesion" substance field if present.
    // Controls focal-adhesion / FAK/Src / anoikis gating.
    double ecm_adhesion   = 1.0; ///< effective adhesion ligand density
    double local_rons     = 0.0; ///< intracellular ROS/RONS (from ros_internal_)
  };
  // -----------------------------------------------------------------------
  // Mechanism 11: Central cell-cycle checkpoint state from one consolidated
  // evaluation. Used to gate all proliferation/survival decisions without
  // re-querying the microenvironment multiple times per step.
  // Molecular annotations are for documentation / parameter naming only.
  struct CellCycleCheckpointState {
    // G1/S gate: ATM/ATR → CHK1/CHK2 → CDC25A inhibition → CDK2/CyclinE arrest
    //            → RB/E2F transcription block
    bool can_enter_S      = true;
    // G2/M gate: ATM/ATR → CHK1/CHK2 → CDC25C inhibition → CDK1/CyclinB arrest
    bool can_enter_M      = true;
    // Growth gate: nutrient/O2/crowding sufficient for biomass increase
    //   PI3K/Akt/mTOR pathway; YAP/TAZ mechanosensing; contact inhibition
    bool can_grow         = true;
    // Division gate: composite (can_enter_M + volume threshold + crowding)
    bool can_divide       = true;
    // G0/quiescence recommended by microenvironment (crowding, hypoxia,
    //   nutrient depletion → p27/Kip1 upregulation; RB hypophosphorylation)
    bool should_enter_G0  = false;
    // Apoptosis commitment (anoikis: BIM via integrin/FAK loss;
    //   or stress-induced intrinsic pathway via p53/BAX/BCL-2)
    bool should_enter_Ap  = false;
    // Necrosis commitment (severe energy collapse; membrane integrity failure)
    bool should_enter_Nec = false;
    // DNA repair permissible (requires O2 for NHEJ; nutrients for synthesis)
    bool repair_allowed   = true;
  };
  // -----------------------------------------------------------------------
  // Mechanism 12: CAP/PAM-specific checkpoint state from EvaluateCAPCheckpointState().
  // Encodes all gate decisions driven by RONS-induced DNA damage, DDR signalling,
  // accumulated arrest time, and apoptosis commitment.
  // Distinct from CellCycleCheckpointState (Mechanism 11).
  struct CAPCheckpointState {
    // G1/S gate: p53/p21-mediated CDK2/CyclinE inhibition after RONS-induced damage.
    // Molecular: ATM/ATR → CHK1/CHK2 → CDC25A degradation → CDK2 inactive
    //            → RB hypophosphorylated → E2F transcription factors blocked
    //            p53 → p21 (CDKN1A) → CDK2/CyclinE inhibition
    bool can_enter_S      = true;
    // G2/M gate: CHK1-CDC25C-CDK1/CyclinB inhibition.
    // Most important for CAP response: γH2AX-positive cells remain arrested
    // until damage resolved. Division with unresolved DSBs = mitotic catastrophe.
    bool can_enter_M      = true;
    // Intra-S slowing: ATR–CHK1–CDC25A ubiquitination during replication stress.
    // Returns true when Sy→G2 should be blocked (damage arises during S phase).
    bool intra_s_blocked  = false;
    // General arrest flag: true when any checkpoint is actively blocking.
    bool must_arrest      = false;
    // DNA repair permitted: NHEJ/HR/NER require O2 and metabolic energy.
    // HIF-1α under hypoxia competes for repair factor binding (RAD51/Ku70/OGG1).
    bool repair_allowed   = true;
    // Apoptosis commitment: accumulated damage/arrest threshold crossed.
    // Models delayed apoptosis observed in EGI-1 (~72h) and HuCCT1 (~48h) CCA lines.
    // Triggered when: apoptosis_commitment_state > threshold
    //              OR arrest_time > max_repair_time AND dna_damage still elevated.
    // Downstream: p53 → BAX/BCL-2 imbalance → cytochrome C → Apaf-1 → caspase-9/3.
    bool must_enter_Ap    = false;
    // Necrosis: extreme energy/membrane collapse under very high RONS or hypoxia.
    // Biologically: high peroxynitrite → mitochondrial membrane rupture → HMGB1/DAMPs.
    bool must_enter_Nec   = false;
    // CAP dose integral exceeded safe threshold: signals chronic treatment exposure
    bool cap_dose_exceeded = false;
  };
//
public:
  BiologicalCell() {}
  explicit BiologicalCell(int p, const bdm::Double3& xyz) : bdm::neuroscience::NeuronSoma(xyz) {
    phenotype_ = p;
    phase_ = BiologicalCell::Phase::I0;
    age_ = 1;
    polarize_ = eye();
    can_apoptose_ = can_grow_ = can_divide_ = can_migrate_ = can_transform_ = can_polarize_ = can_protrude_ = false;
    trail_ = 0.0;
    n_divisions_ = n_trasformations_ = n_protrusions_ = 0;
    params_ = 0; // nullify pointer...
  }
  //
  void Initialize(const bdm::NewAgentEvent& event) override {
    NeuronSoma::Initialize(event);
    // if cell divides then attributes have to be initialized
    if (auto* mother = dynamic_cast<BiologicalCell*>(event.existing_agent))
      {
        if (event.GetUid() == bdm::CellDivisionEvent::kUid)
          {
            phenotype_ = mother->GetPhenotype();
            phase_ = BiologicalCell::Phase::I0;
            SetAge(); // ...age cannot be inherited
            polarize_ = mother->GetPolarization();
            can_apoptose_  = mother->GetCanApoptose();
            can_grow_      = mother->GetCanGrow();
            can_divide_    = mother->GetCanDivide();
            can_migrate_   = mother->GetCanMigrate();
            can_transform_ = mother->GetCanTransform();
            can_polarize_  = mother->GetCanPolarize();
            can_protrude_  = mother->GetCanProtrude();
            ResetTrail(); // ...trail cannot be inherited
            n_divisions_ = 0; mother->IncrementNumberOfDivisions();
            n_trasformations_ = 0; // ...index is initialized
            n_protrusions_ = 0; // ...index is initialized
            params_ = mother->params_; // copy parameters pointer...
            // daughter cells start fresh on phase timer and arrest state
            phase_age_ = 0;
            arrest_time_ = 0;
            is_quiescent_ = false;
            regulatory_backend_id_ = mother->regulatory_backend_id_;
            regulatory_model_phenotype_id_ = mother->regulatory_model_phenotype_id_;
            regulatory_update_interval_ = mother->regulatory_update_interval_;
            regulatory_update_counter_ = 0;
            regulatory_model_type_ = mother->regulatory_model_type_;
            regulatory_parameters_ = mother->regulatory_parameters_;
            regulatory_state_ = mother->regulatory_state_;
            regulatory_output_ = mother->regulatory_output_;
            // Inherit a configurable fraction of parent's intracellular ROS and DNA damage.
            // Biologically, daughter cells can receive some oxidative burden and unrepaired
            // lesions from the parent. Default 0.0 = daughters start clean (backward-compatible).
            {
              const std::string& CP_p =
                params_->get<std::string>("phenotype_ID/"+std::to_string(phenotype_));
              const double frac =
                params_->have_parameter<double>(CP_p+"/intracellular/damage/daughter_partition_factor")
                ? params_->get<double>(CP_p+"/intracellular/damage/daughter_partition_factor")
                : 0.0;
              ros_internal_         = mother->ros_internal_ * frac;
              antioxidant_capacity_ = 1.0;
              dna_damage_           = mother->dna_damage_   * frac;
              atm_active_           = mother->atm_active_   * frac;
              atr_active_           = mother->atr_active_   * frac;
              chk1_active_          = mother->chk1_active_  * frac;
              chk2_active_          = mother->chk2_active_  * frac;
              p53_active_           = mother->p53_active_   * frac;
              p21_level_            = mother->p21_level_    * frac;
              cdc25_active_         = mother->cdc25_active_ * frac
                                      + (1.0 - frac) * 1.0;
              cdk_activity_         = mother->cdk_activity_ * frac
                                      + (1.0 - frac) * 1.0;
              // Mechanism 12 — CAP/PAM state inheritance
              // RNS: daughter inherits same fraction as ROS
              rns_internal_             = mother->rns_internal_             * frac;
              // Dose integral: daughter starts fresh (the treatment was applied to the parent)
              cap_dose_integral_        = 0.0;
              // Time since CAP: inherited to preserve treatment timing context
              time_since_cap_           = mother->time_since_cap_;
              // Marker proxies start clean (proxies are computed, not inherited)
              oxidative_damage_8oxoG_proxy_  = mother->oxidative_damage_8oxoG_proxy_  * frac;
              dsb_damage_gammaH2AX_proxy_    = mother->dsb_damage_gammaH2AX_proxy_    * frac;
              // Execution markers reset to zero — daughters start uncommitted
              parp_cleavage_proxy_      = 0.0;
              caspase3_activation_proxy_= 0.0;
              apoptosis_commitment_state_= mother->apoptosis_commitment_state_ * frac;
              // Physical/capacity parameters: inherited fully from parent phenotype
              membrane_permeability_    = mother->membrane_permeability_;
              repair_capacity_          = mother->repair_capacity_;
              // Arrest phase resets: daughters are not checkpoint-arrested at birth
              cap_arrest_phase_         = 0;
              cap_recovered_from_arrest_count_ = 0;
            }
            // Rebuild backend selection for the daughter phenotype and parameters.
            ConfigureRegulatoryModel();
            CheckAndFixDiameter(); mother->CheckAndFixDiameter();
          }
        else
          ABORT_("an exception is caught");
      }
  }
  //
  void SetPhenotype(int p) { phenotype_ = p; }
  int GetPhenotype() const { return phenotype_; }
  //
  void SetPhase(int p) { phase_ = static_cast<BiologicalCell::Phase>(p); }
  int GetPhase() const { return phase_; }
  //
  void SetAge(unsigned int a =1) { age_ = a; }
  int  GetAge() const { return age_; }
  void IncrementAge() { age_++; }
  //
  void IncrementPhaseAge() { ++phase_age_; }
  int  GetPhaseAge() const { return phase_age_; }
  void ResetPhaseAge() { phase_age_ = 0; }
  //
  void IncrementArrestTime() { ++arrest_time_; }
  int  GetArrestTime() const { return arrest_time_; }
  void ResetArrestTime() { arrest_time_ = 0; }
  //
  void SetQuiescent(bool q) { is_quiescent_ = q; }
  bool IsQuiescent() const { return is_quiescent_; }
  //
  double GetROSInternal() const { return ros_internal_; }
  double GetDNADamage() const { return dna_damage_; }
  double GetAntioxidantCapacity() const { return antioxidant_capacity_; }
  double GetAtmActive() const { return atm_active_; }
  double GetAtrActive() const { return atr_active_; }
  double GetChk1Active() const { return chk1_active_; }
  double GetChk2Active() const { return chk2_active_; }
  double GetP53Active() const { return p53_active_; }
  double GetP21Level() const { return p21_level_; }
  double GetCdc25Active() const { return cdc25_active_; }
  double GetCdkActivity() const { return cdk_activity_; }
  const regulatory::RegulatoryOutput& GetRegulatoryOutput() const {
    return regulatory_output_;
  }
  double GetRegulatoryNodeActivity(const std::string& node_name) const;
  bool IsRegulatoryModelActive() const { return regulatory_backend_id_ != 0; }
  std::string GetRegulatoryModelType() const { return regulatory_model_type_; }
  // Mechanism 12 — CAP/PAM-specific state getters
  double GetRNSInternal() const { return rns_internal_; }
  double GetCapDoseIntegral() const { return cap_dose_integral_; }
  int    GetTimeSinceCap() const { return time_since_cap_; }
  double GetOxidativeDamage8OxoGProxy() const { return oxidative_damage_8oxoG_proxy_; }
  double GetDSBDamageGammaH2AXProxy() const { return dsb_damage_gammaH2AX_proxy_; }
  double GetPARPCleavageProxy() const { return parp_cleavage_proxy_; }
  double GetCaspase3ActivationProxy() const { return caspase3_activation_proxy_; }
  double GetApoptosisCommitmentState() const { return apoptosis_commitment_state_; }
  double GetMembranePermeability() const { return membrane_permeability_; }
  double GetRepairCapacity() const { return repair_capacity_; }
  void   SetAntioxidantCapacity(double value) {
    antioxidant_capacity_ = std::max(0.0, value);
  }
  void   SetRepairCapacity(double value) {
    repair_capacity_ = std::clamp(value, 0.0, 1.0);
  }
  void   SetCapArrestPhase(int phase) { cap_arrest_phase_ = phase; }
  int    GetCapArrestPhase() const { return cap_arrest_phase_; }
  void   IncrementCapRecoveredCount() { ++cap_recovered_from_arrest_count_; }
  int    GetCapRecoveredCount() const { return cap_recovered_from_arrest_count_; }
  void   ResetCapRecoveredCount() { cap_recovered_from_arrest_count_ = 0; }
  void   IncrementTimeSinceCap() { ++time_since_cap_; }
  //
  void SetPolarization(const bdm::Double3x3& p) { polarize_ = p; }
  const bdm::Double3x3& GetPolarization() const { return polarize_; }
  const double& GetPolarization(size_t i, size_t j) const { return polarize_[i][j]; }
  //
  void SetCanApoptose(bool apoptoses) { can_apoptose_ = apoptoses; }
  bool GetCanApoptose() const { return can_apoptose_; }
  //
  void SetCanGrow(bool grows) { can_grow_ = grows; }
  bool GetCanGrow() const { return can_grow_; }
  //
  void SetCanDivide(bool divides) { can_divide_ = divides; }
  bool GetCanDivide() const { return can_divide_; }
  //
  void SetCanMigrate(bool migrates) { can_migrate_ = migrates; }
  bool GetCanMigrate() const { return can_migrate_; }
  //
  void SetCanTransform(bool transforms) { can_transform_ = transforms; }
  bool GetCanTransform() const { return can_transform_; }
  //
  void SetCanPolarize(bool polarizes) { can_polarize_ = polarizes; }
  bool GetCanPolarize() const { return can_polarize_; }
  //
  void SetCanProtrude(bool protrudes) { can_protrude_ = protrudes; }
  bool GetCanProtrude() const { return can_protrude_; }
  //
  void ResetTrail() { trail_ = 0.0; }
  void UpdateTrail(double d) { trail_ += d; }
  double GetTrail() const { return trail_; }
  //
  const bdm::Double3& GetActiveDisplacement () const { return active_displacement_; }
  const bdm::Double3& GetPassiveDisplacement() const { return passive_displacement_; }
  const bdm::Double3 GetDisplacement() const { return active_displacement_+passive_displacement_; }
  const double GetDisplacement(size_t i) const { return active_displacement_[i]+passive_displacement_[i]; }
  //
  void IncrementNumberOfDivisions() { ++n_divisions_; }
  int GetNumberOfDivisions() const { return n_divisions_; }
  //
  void IncrementNumberOfTrasformations() { ++n_trasformations_; }
  int GetNumberOfTrasformations() const { return n_trasformations_; }
  //
  void IncrementNumberOfProtrusions() { ++n_protrusions_; }
  int GetNumberOfProtrusions() const { return n_protrusions_; }
  //
  void SetParametersPointer(Parameters* p) { params_ = p; }
  Parameters* params() const { return params_; }
  //
  void RunBiochemics();
  void RunIntracellular();
  bool RunECMInteraction();      // returns true if anoikis was triggered
  bool EvaluateG1SCheckpoint();  // returns true if G1->S transition is blocked
  bool EvaluateG2MCheckpoint();  // returns true if G2->M transition is blocked
  bool CheckNecrosis();          // returns true if necrotic transformation occurred
  bool CheckApoptosisByDamage();
  bool CheckPositionValidity();
  bool CheckApoptosisAging();
  bool CheckApoptosis();
  bool CheckAfterApoptosis();
  bool CheckQuiescenceAfterDivision();
  bool CheckMigration();
  bool CheckTransformation();
  bool CheckPolarization();
  bool CheckProtrusion();
  bool CheckGrowth();
  bool CheckTransformationAndDivision();
  bool CheckAsymmetricDivision();
  bool CheckDivision();
  // Returns phi = occupied_area_or_volume / influence_area_or_volume in [0, inf).
  // influence_ratio scales the cell diameter to define the local search radius.
  // Returns 0.0 if influence_ratio <= 0.
  double ComputeLocalOccupancyRatio(const bdm::Double3& position, double influence_ratio) const;
  void Set2DeleteProtrusions();
  // Mechanism 11 — new supporting methods (backward-compatible, all optional)
  // -----------------------------------------------------------------------
  /// Sample local microenvironment fields once per timestep. Safe defaults
  /// are returned for any field whose substance/grid is absent, so all
  /// existing input files continue to run unchanged.
  MicroenvironmentState SampleMicroenvironment() const;
  /// Central cell-cycle checkpoint controller (Mechanism 11).
  /// Uses a pre-sampled MicroenvironmentState to avoid redundant grid queries.
  /// Returns a CellCycleCheckpointState encoding all gate decisions.
  CellCycleCheckpointState EvaluateCellCycleCheckpoints(
      const MicroenvironmentState& env);
  /// Intra-S checkpoint: returns true if Sy→G2 should be blocked.
  /// Phenomenological implementation of ATR–CHK1–CDC25A inhibition during
  /// S-phase when replication stress or DNA DSBs are elevated.
  bool EvaluateIntraSCheckpoint();
  /// Multiplicative growth-rate modifier ∈ [0,1] encoding microenvironmental
  /// support for biomass accumulation. Returns 1.0 if no modulation
  /// parameters are defined (backward-compatible).
  /// Factors: O2 (mTOR/HIF-1α), nutrient (PI3K/Akt/mTOR), ECM stiffness
  /// (YAP/TAZ mechanosensing), crowding (contact inhibition), ROS/stress.
  double ComputeGrowthModulation(const MicroenvironmentState& env) const;
  // -----------------------------------------------------------------------
  // Mechanism 12 — CAP/PAM-specific methods
  // -----------------------------------------------------------------------
  /// Update CAPP/RONS intracellular dynamics for Mechanism 12:
  ///   - rns_internal_ (reactive nitrogen species, NO2-/peroxynitrite tracking)
  ///   - cap_dose_integral_ (area-under-curve of extracellular RONS)
  ///   - time_since_cap_ increment
  ///   - Marker proxies: 8-oxoG, γH2AX (DSB), pCHK1, p53, PARP, caspase-3
  ///   - apoptosis_commitment_state_ accumulation
  ///   - membrane_permeability_ and repair_capacity_ updates
  ///
  /// Must be called AFTER RunIntracellular() so that ros_internal_,
  /// dna_damage_, and DDR pathway variables are already up-to-date.
  ///
  /// Biological model:
  ///   ROS accumulates from H2O2/ONOO- uptake (aquaporins AQP3/AQP8).
  ///   DNA damage includes 8-oxoG (base oxidation) and DSBs (strand breaks).
  ///   ATM (DSB sensor) and ATR (ssDNA/fork stall sensor) activate CHK1/CHK2.
  ///   γH2AX marks DSB sites; pCHK1 and p53-Ser15 phosphorylation are DDR markers.
  ///   PARP cleavage and caspase-3 activation are apoptosis execution markers.
  ///   Cell-line sensitivity is encoded in phenotype parameter sets, not hardcoded.
  void UpdateCAPIntracellular();
  /// Central checkpoint controller for CAP/PAM-treated cells (Mechanism 12).
  /// Returns a CAPCheckpointState encoding:
  ///   - G1/S gate (p53/p21-CDK2 axis)
  ///   - G2/M gate (CHK1-CDC25C-CDK1 axis — primary arrest gate for CAP)
  ///   - Intra-S slowing (ATR-CHK1 during replication)
  ///   - Repair allowance (O2/nutrient requirements for NHEJ/HR)
  ///   - Apoptosis commitment (delayed, dose/time-dependent)
  ///   - Necrosis flag (extreme RONS/energy collapse)
  ///
  /// Distinct from EvaluateCellCycleCheckpoints() (Mechanism 11).
  CAPCheckpointState EvaluateCAPCheckpointState(
      const MicroenvironmentState& env);
  // -----------------------------------------------------------------------
  // Modular intracellular regulation interface (Boolean / stochastic Boolean)
  // -----------------------------------------------------------------------
  void ConfigureRegulatoryModel();
  void UpdateRegulatoryModel(const MicroenvironmentState& env);
  regulatory::RegulatoryInput BuildRegulatoryInput(
      const MicroenvironmentState& env) const;
  //
//
private:
  //
  void CheckAndFixDiameter();
  bool CheckProtrusionAxis(bdm::Double3 axis);
  //
//
private:
  // index to designate the cell phenotype
  int phenotype_ = 0; // WARNING: phenotype ID must be >=0
  // cell circle phase
  BiologicalCell::Phase phase_ = BiologicalCell::Phase::I0;
  // cell age (non-fractional time)
  int age_ = 1;
  // cell polarization axes
  bdm::Double3x3 polarize_ = eye();
  // flags to designate (individual) cell behaviour
  bool can_apoptose_, can_grow_, can_divide_, can_migrate_, can_transform_, can_polarize_, can_protrude_;
  // total cell trail (displacement) between user-defined time points
  double trail_ = 0.0;
  bdm::Double3 active_displacement_ = {0.0, 0.0, 0.0};
  bdm::Double3 passive_displacement_ = {0.0, 0.0, 0.0};
  // index to keep track of the (individual) cell divisions & trasformations
  // and total number of filopodium or/and neurite (outgrowth) protrusions
  int n_divisions_ = 0, n_trasformations_ = 0, n_protrusions_ = 0;
  // pointer to all simulation parameters
  mutable
  Parameters* params_ = 0;
  // intracellular states
  double ros_internal_ = 0.0;
  double antioxidant_capacity_ = 1.0;
  double dna_damage_ = 0.0;
  // DDR pathway (ATM/ATR–CHK–p53–p21–Cdc25–CDK), normalized activity in [0, 1]
  double atm_active_ = 0.0;
  double atr_active_ = 0.0;
  double chk1_active_ = 0.0;
  double chk2_active_ = 0.0;
  double p53_active_ = 0.0;
  double p21_level_ = 0.0;
  double cdc25_active_ = 1.0;
  double cdk_activity_ = 1.0;
  friend void UpdateDdrPathway(BiologicalCell* cell);
  // cell-cycle phase timer: number of time steps spent in current phase
  int phase_age_ = 0;
  // checkpoint arrest timer: number of time steps spent arrested at a checkpoint
  int arrest_time_ = 0;
  // G0 quiescence flag: true when cell has entered nutrient/crowding-driven quiescence
  bool is_quiescent_ = false;
  // -----------------------------------------------------------------------
  // Mechanism 12 — CAP/PAM-specific intracellular state variables
  // -----------------------------------------------------------------------
  // rns_internal_: intracellular reactive nitrogen species (NO2-/peroxynitrite)
  // Biologically: NO2- from CAP/PAM enters via aquaporins and can react with
  // superoxide O2·- to form peroxynitrite (ONOO-), which nitrosylates proteins
  // and DNA bases (8-nitroguanine, DNA strand nicking).
  double rns_internal_ = 0.0;
  // cap_dose_integral_: time-integrated extracellular RONS exposure (dose AUC).
  // Represents cumulative reactive species burden analogous to LQ dose integral.
  // Used to model dose-dependent survival reduction observed in EGI-1 / HuCCT1.
  double cap_dose_integral_ = 0.0;
  // time_since_cap_: time steps since first CAP/PAM exposure (non-zero RONS detected).
  // Used to model time-dependent delayed apoptosis (EGI-1 ~72h, HuCCT1 ~48h
  // under tested PAM conditions). Cell-line timing encoded via sensitivity params.
  int time_since_cap_ = 0;
  // oxidative_damage_8oxoG_proxy_: normalised 8-oxoguanine level ∈ [0,1].
  // 8-oxoG arises from OH· (Fenton/Haber-Weiss) oxidation of guanine in DNA.
  // Reported by IHC (anti-8-OHdG antibody) in in vivo CAP-treated CCA xenografts.
  // Drives NER (OGG1/APE1/PCNA base excision repair pathway).
  double oxidative_damage_8oxoG_proxy_ = 0.0;
  // dsb_damage_gammaH2AX_proxy_: normalised γH2AX foci count ∈ [0,1].
  // γH2AX (H2AX-Ser139 phosphorylation) is the canonical DSB marker.
  // Rapid ATM/ATR-mediated phosphorylation at DSB sites (≤1h after damage).
  // Reported positive in CAP/PAM-treated EGI-1 and HuCCT1 cells (immunofluorescence).
  double dsb_damage_gammaH2AX_proxy_ = 0.0;
  // parp_cleavage_proxy_: normalised cleaved PARP-1 (89 kDa fragment) level ∈ [0,1].
  // PARP-1 is cleaved by executioner caspase-3 (-7) during apoptosis.
  // Cleavage inactivates PARP-1's DNA repair function, sealing apoptosis fate.
  // Reported as western blot marker in CAP-treated CCA apoptosis assays.
  // Only increases after apoptosis commitment (not during DDR/arrest phase).
  double parp_cleavage_proxy_ = 0.0;
  // caspase3_activation_proxy_: normalised cleaved caspase-3 level ∈ [0,1].
  // Caspase-3 is the primary executioner caspase; activated by caspase-9
  // downstream of cytochrome-C/Apaf-1 apoptosome (intrinsic/mitochondrial path).
  // Reported by IHC in in vivo CAP-treated CCA xenografts (cleaved caspase-3 Ab).
  // Used as definitive apoptosis readout alongside Annexin V / 7-AAD staining.
  double caspase3_activation_proxy_ = 0.0;
  // apoptosis_commitment_state_: continuous commitment accumulation ∈ [0,1].
  // Integrates pro-apoptotic signals over time (damage history, arrest duration).
  // Once commitment crosses apoptosis_commitment_threshold, cell enters Ap phase.
  // Models delayed apoptosis: arrest appears first, apoptosis occurs later.
  // Represents p53/BAX–BCL-2 imbalance slowly tipping toward irreversible fate.
  double apoptosis_commitment_state_ = 0.0;
  // membrane_permeability_: RONS uptake modulator ∈ [0,1].
  // Reflects aquaporin (AQP3/AQP8) density and lipid bilayer integrity.
  // Reduced by membrane damage or lipid peroxidation from high ROS.
  // Models inter-cell variability in RONS uptake (contributes to dose-response variance).
  double membrane_permeability_ = 1.0;
  // repair_capacity_: cell's active DNA repair ability ∈ [0,1].
  // Aggregates NHEJ (Ku70/Ku80/DNA-PKcs for DSBs), BER (OGG1/APE1 for 8-oxoG),
  // and HR (RAD51/BRCA1/BRCA2 for complex DSBs in S/G2 phase).
  // Reduced by sustained oxidative damage to repair enzymes.
  // Biologically, primary hepatocytes have higher repair capacity than CCA cells.
  // Phenotype-configurable: repair_capacity is a parameter, not a cell-line name.
  double repair_capacity_ = 1.0;
  // cap_arrest_phase_: tracks which checkpoint is responsible for current arrest.
  //   0 = not checkpoint-arrested
  //   1 = arrested at G1/S (p53/p21 gate; CDK2/CyclinE inhibition)
  //   2 = arrested at intra-S (ATR–CHK1–CDC25A; replication stress)
  //   3 = arrested at G2/M (CHK1/CHK2–CDC25C–CDK1/CyclinB inhibition)
  // Used for statistics export (number of cells blocked at each checkpoint).
  int cap_arrest_phase_ = 0;
  // cap_recovered_from_arrest_count_: number of successful recoveries where a
  // checkpoint-arrested cell repaired sufficiently and re-entered progression.
  int cap_recovered_from_arrest_count_ = 0;
  // -----------------------------------------------------------------------
  // Modular intracellular regulatory layer (GRN)
  // -----------------------------------------------------------------------
  int regulatory_backend_id_ = 0;  // 0 none, 1 boolean, 2 stochastic_boolean
  int regulatory_model_phenotype_id_ = -1;
  int regulatory_update_interval_ = 1;
  int regulatory_update_counter_ = 0;
  std::string regulatory_model_type_ = "none";
  regulatory::RegulatoryParameters regulatory_parameters_;
  regulatory::RegulatoryState regulatory_state_;
  regulatory::RegulatoryOutput regulatory_output_;
  regulatory::NullRegulatoryModel regulatory_model_null_;
  regulatory::BooleanRegulatoryModel regulatory_model_boolean_;
  regulatory::StochasticBooleanRegulatoryModel regulatory_model_stochastic_;
  // list of cell protrusions (filopodia or neurites)
  std::vector<bdm::Double3> protrusions_;
};
// =============================================================================
} // ...end of namespace
// =============================================================================
#endif // _BIOLOGICAL_CELL_H_
// =============================================================================
