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
            }
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
  // list of cell protrusions (filopodia or neurites)
  std::vector<bdm::Double3> protrusions_;
};
// =============================================================================
} // ...end of namespace
// =============================================================================
#endif // _BIOLOGICAL_CELL_H_
// =============================================================================
