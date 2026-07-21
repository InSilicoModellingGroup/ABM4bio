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
  
  // =============================================================================
  // Cell-matrix mechanics: persistent attachment record
  // =============================================================================

  struct AttachmentRecord {
    // Persistent one-based scaffold/FEM node identifier.
    int node_id = -1;

    // Current attachment coordinate on the deformed scaffold.
    bdm::Double3 position = {0.0, 0.0, 0.0};

    // Most recently calculated local scaffold stiffness, k_ecm.
    double k_ecm = 0.0;

    // True only after k_ecm has been calculated by the FEM.
    bool has_valid_k_ecm = false;

    // True when the attachment was formed by the ABM after the previous FEM step.
    bool newly_formed = false;
  };

  // =============================================================================
  // Cell-matrix mechanics: attachment-state queries
  // =============================================================================

  std::size_t GetNumberOfAttachments() const {
    /*
    * Function goal
    * -------------
    * Return the number of persistent attachment records stored by the cell.
    */

    return attachment_records_.size();
  }


  bool HasAttachments() const {
    /*
    * Function goal
    * -------------
    * Return true when the cell has at least one persistent attachment record.
    *
    * A valid k_ecm value is not required because newly formed attachments have
    * not yet been processed by the FEM.
    */

    return !attachment_records_.empty();
  }

  
  
  // -----------------------------------------------------------------------------
  // Clear all attachment data
  // -----------------------------------------------------------------------------

  void ClearAttachments() {
  /*
   * Function goal
   * -------------
   * Remove all persistent attachments from this cell.
   */

  attachment_records_.clear();
}

  // =============================================================================
  // Cell-matrix mechanics: persistent attachment-record interface
  // =============================================================================

  // -----------------------------------------------------------------------------
  // Store attachment records
  // -----------------------------------------------------------------------------

  void SetAttachmentRecords(
      const std::vector<AttachmentRecord>& records) {
    /*
    * Function goal
    * -------------
    * Store the complete persistent attachment state for this cell.
    *
    * Each record keeps the scaffold node ID, attachment coordinate and local
    * stiffness together, preventing the attachment data from becoming
    * misaligned.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the records
    // -------------------------------------------------------------------------

    for (std::size_t i = 0; i < records.size(); ++i) {
      ASSERT_(
        records[i].node_id > 0,
        "AttachmentRecord contains a non-positive scaffold node ID"
      );

      // A cell must not attach to the same scaffold node more than once.
      for (std::size_t j = i + 1; j < records.size(); ++j) {
        ASSERT_(
          records[i].node_id != records[j].node_id,
          "AttachmentRecord contains duplicate scaffold node ID "
          + std::to_string(records[i].node_id)
        );
      }
    }

    // -------------------------------------------------------------------------
    // Step 2: Replace the current attachment state
    // -------------------------------------------------------------------------

    attachment_records_ = records;
  }

  // -----------------------------------------------------------------------------
  // Retrieve attachment records
  // -----------------------------------------------------------------------------

  const std::vector<AttachmentRecord>& GetAttachmentRecords() const {
    return attachment_records_;
  }

  const AttachmentRecord& GetAttachmentRecord(
      const std::size_t index) const {
    /*
    * Function goal
    * -------------
    * Retrieve one attachment record using its position in the cell's
    * attachment-record vector.
    */

    ASSERT_(
      index < attachment_records_.size(),
      "AttachmentRecord index is out of range"
    );

    return attachment_records_[index];
  }

  std::size_t GetNumberOfAttachmentRecords() const {
    return attachment_records_.size();
  }


  // -----------------------------------------------------------------------------
  // Attachment-record state queries
  // -----------------------------------------------------------------------------

  bool HasAttachmentRecords() const {
    return !attachment_records_.empty();
  }

  bool HaveAllAttachmentRecordsValidKecm() const {
    /*
    * Function goal
    * -------------
    * Return true when every stored attachment record has a FEM-calculated
    * k_ecm value.
    */

    if (attachment_records_.empty()) {
      return false;
    }

    for (const auto& record : attachment_records_) {
      if (!record.has_valid_k_ecm) {
        return false;
      }
    }

    return true;
  }

  bool HasValidAttachmentRecordsForMechanics() const {
    /*
    * Function goal
    * -------------
    * Return true when the cell has attachment records and every record contains
    * the valid FEM-calculated mechanics data required for contraction.
    */

    // A cell without attachment records cannot contract.
    if (attachment_records_.empty()) {
      return false;
    }

    // Every attachment must have a valid persistent ID and k_ecm value.
    for (const auto& record : attachment_records_) {
      if (record.node_id <= 0) {
        return false;
      }

      if (!record.has_valid_k_ecm) {
        return false;
      }
    }

    return true;
  }



  // -----------------------------------------------------------------------------
  // Clear attachment records
  // -----------------------------------------------------------------------------

  void ClearAttachmentRecords() {
    attachment_records_.clear();
  }

  // Matrix stiffness perceived by a cell
  void SetKce(double k_ce) {
    k_ce_ = k_ce;
  }
  double GetKce() const {
    return k_ce_;
  }
  // Contractile force of the cell is adaptable
  void SetContractileForce(double contractile_force) {
    contractile_force_ = contractile_force;
  }
  double GetContractileForce() const {
    return contractile_force_;
  }
  // Cells moved due to mechanics
  void SetMovedDueToMechanics(bool moved_due_to_mechanics) {
    moved_due_to_mechanics_ = moved_due_to_mechanics;
  }
  bool GetMovedDueToMechanics() const {
    return moved_due_to_mechanics_;
  }
  void ClearMovedDueToMechanics() {
    moved_due_to_mechanics_ = false;
  }
  // Cell state for coupled fem solver (cell-matrix mechanics)
  void SetCellState(const std::string& state) {
    ASSERT_(
        state == "attach" || state == "contract",
        "Invalid cell_state. Expected 'attach' or 'contract'."
    );
    cell_state_ = state;
  }
  const std::string& GetCellState() const {
    return cell_state_;
  }
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
  bdm::Double3 GetActiveDisplacement () const { return active_displacement_; }
  bdm::Double3 GetPassiveDisplacement() const { return passive_displacement_; }
  bdm::Double3 GetDisplacement() const { return active_displacement_+passive_displacement_; }
  double GetDisplacement(size_t i) const { return active_displacement_[i]+passive_displacement_[i]; }
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
  // list of cell protrusions (filopodia or neurites)
  std::vector<bdm::Double3> protrusions_;
  // FEM mechanics state
  std::string cell_state_ = "attach";
  double k_ce_ = 0.0;
  double contractile_force_ = 0.0;
  bool moved_due_to_mechanics_ = false;
  
  // =============================================================================
  // Cell-matrix mechanics: attachment data
  // =============================================================================

  // Persistent attachment records.
  //
  // Each record stores the one-based scaffold node ID, current attachment
  // coordinate and associated local mechanics state.
  std::vector<AttachmentRecord> attachment_records_;


};
// =============================================================================
} // ...end of namespace
// =============================================================================
#endif // _BIOLOGICAL_CELL_H_
// =============================================================================
