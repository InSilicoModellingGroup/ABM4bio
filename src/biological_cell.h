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

  // =============================================================================
  // Cell-matrix mechanics: attachment lifecycle
  // =============================================================================

  enum class CellMatrixLifecycleStatus {
    kNeedsInitialAttachment,
    kEstablished,
    kInitialAttachmentRetry,
    kMarkedForDeath
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
  // Cell-matrix mechanics: internal lifecycle state
  // =============================================================================

  void SetCellMatrixLifecycleStatus(
      const CellMatrixLifecycleStatus status) {
    /*
     * Function goal
     * -------------
     * Store the internal ABM attachment-lifecycle status of this cell.
     *
     * This status is separate from cell_state_, which is used only for FEM
     * communication through the "attach" and "contract" strings.
     */

    cell_matrix_lifecycle_status_ = status;
  }

  CellMatrixLifecycleStatus GetCellMatrixLifecycleStatus() const {
    return cell_matrix_lifecycle_status_;
  }

  bool RequiresMechanicsRecalculation() const {
    /*
     * Function goal
     * -------------
     * Return whether the current attachment set requires a new FEM mechanics
     * calculation.
     */

    return requires_mechanics_recalculation_;
  }

  void MarkMechanicsForRecalculation() {
    /*
     * Function goal
     * -------------
     * Record that the attachment set has changed and its cell-level mechanics
     * must be recalculated by the FEM.
     *
     * Existing k_ce, contractile force and retained attachment k_ecm values are
     * deliberately preserved until the recalculation is completed.
     */

    requires_mechanics_recalculation_ = true;
  }

  void ClearMechanicsRecalculationRequirement() {
    /*
     * Function goal
     * -------------
     * Record that the current attachment set has completed its FEM mechanics
     * calculation.
     */

    requires_mechanics_recalculation_ = false;
  }

  // -----------------------------------------------------------------------------
  // Validate the internal cell-matrix state
  // -----------------------------------------------------------------------------

  void ValidateCellMatrixState() const {
    /*
     * Function goal
     * -------------
     * Confirm that the attachment lifecycle, attachment records and mechanics
     * recalculation state are internally consistent.
     */

    const CellMatrixLifecycleStatus lifecycle_status =
        cell_matrix_lifecycle_status_;

    const std::size_t attachment_count =
        attachment_records_.size();

    // -------------------------------------------------------------------------
    // Step 1: Validate cells awaiting initial FEM attachment
    // -------------------------------------------------------------------------

    if (lifecycle_status ==
            CellMatrixLifecycleStatus::kNeedsInitialAttachment ||
        lifecycle_status ==
            CellMatrixLifecycleStatus::kInitialAttachmentRetry) {

      ASSERT_(
          attachment_count == 0,
          "A cell awaiting initial FEM attachment cannot already contain "
          "persistent attachment records"
      );

      ASSERT_(
          !requires_mechanics_recalculation_,
          "A cell awaiting initial FEM attachment cannot require mechanics "
          "recalculation"
      );
    }

    // -------------------------------------------------------------------------
    // Step 2: Validate established cells
    // -------------------------------------------------------------------------

    if (lifecycle_status ==
        CellMatrixLifecycleStatus::kEstablished) {

      ASSERT_(
          attachment_count > 0,
          "An established cell must retain at least one persistent attachment"
      );
    }

    // -------------------------------------------------------------------------
    // Step 3: Validate individual attachment records
    // -------------------------------------------------------------------------

    for (const auto& record : attachment_records_) {
      ASSERT_(
          record.node_id > 0,
          "Cell-matrix state contains a non-positive attachment node ID"
      );

      ASSERT_(
          !(record.newly_formed && record.has_valid_k_ecm),
          "A newly formed ABM attachment cannot already contain valid FEM k_ecm"
      );

      if (record.newly_formed) {
        ASSERT_(
            requires_mechanics_recalculation_,
            "A cell with a newly formed attachment must require mechanics "
            "recalculation"
        );
      }
    }
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
  // Update established attachments from the ABM
  // -----------------------------------------------------------------------------

  void UpdateAttachmentRecordsFromAbm(
      const std::vector<AttachmentRecord>& records) {
    /*
     * Function goal
     * -------------
     * Replace an established cell's attachment set after an ABM-controlled
     * attachment change.
     *
     * Existing k_ce, contractile force and retained attachment k_ecm values are
     * preserved. The new attachment set is marked as requiring a FEM mechanics
     * recalculation.
     */

    // -------------------------------------------------------------------------
    // Step 1: Confirm that the ABM owns this attachment set
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Only an established cell may update attachments through the ABM"
    );

    // A previously attached cell must always retain at least one attachment.
    ASSERT_(
        !records.empty(),
        "An established cell cannot lose all persistent attachments"
    );

    // -------------------------------------------------------------------------
    // Step 2: Store the updated attachment records
    // -------------------------------------------------------------------------

    SetAttachmentRecords(records);

    // -------------------------------------------------------------------------
    // Step 3: Mark the cell-level mechanics as outdated
    // -------------------------------------------------------------------------

    MarkMechanicsForRecalculation();
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

  void ClearAttachmentRecords() {
    /*
    * Function goal
    * -------------
    * Remove all persistent attachment records from this cell.
    */

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

  // -----------------------------------------------------------------------------
  // Follow a retained single scaffold attachment
  // -----------------------------------------------------------------------------

  void FollowSingleAttachmentScaffold(
      const ObstacleScaffold& scaffold,
      const bdm::Double3& previous_attachment_position);
  
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
  bdm::Double3 CalculateSingleAttachmentStrutDirection(
      const ObstacleScaffold& scaffold,
      int attachment_node_id) const;
  bdm::Double3 CalculateSingleAttachmentRadialDirection(
    const bdm::Double3& previous_attachment_position,
    const bdm::Double3& strut_direction) const;
  bdm::Double3 CalculateSingleAttachmentTargetPosition(
    const ObstacleScaffold& scaffold,
    int attachment_node_id,
    const bdm::Double3& radial_direction) const;
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
  
  // =============================================================================
  // Cell-matrix mechanics: attachment data
  // =============================================================================
  
  // FEM communication state.
  std::string cell_state_ = "attach";

  // Internal ABM attachment lifecycle.
  CellMatrixLifecycleStatus cell_matrix_lifecycle_status_ =
      CellMatrixLifecycleStatus::kNeedsInitialAttachment;

  // True when the current attachment set requires a new FEM mechanics result.
  bool requires_mechanics_recalculation_ = false;

  double k_ce_ = 0.0;
  double contractile_force_ = 0.0;

  // Persistent attachment records.
  //
  // Each record stores the one-based scaffold node ID, current attachment
  // coordinate and associated local mechanics state.
  std::vector<AttachmentRecord> attachment_records_;

  inline
  double CalculateReferenceDetachmentProbability(
      const double k_ecm,
      const double optimum_k_ecm,
      const double low_k_ecm_slope,
      const double high_k_ecm_slope,
      const double min_probability,
      const double max_probability,
      const double k_ecm_weight) const {
    /*
     * Function goal
     * -------------
     * Calculate the reference detachment probability for one attachment from
     * its local scaffold stiffness (k_ecm).
     *
     * Detachment is least likely near the optimum stiffness and increases
     * towards the configured limit at lower or higher stiffness values (biphasic).
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    ASSERT_(
        k_ecm >= 0.0,
        "Detachment probability calculation requires non-negative k_ecm"
    );

    ASSERT_(
        optimum_k_ecm > 0.0,
        "Detachment probability calculation requires positive optimum k_ecm"
    );

    ASSERT_(
        low_k_ecm_slope > 0.0,
        "Detachment probability calculation requires a positive low-k_ecm "
        "slope"
    );

    ASSERT_(
        high_k_ecm_slope > 0.0,
        "Detachment probability calculation requires a positive high-k_ecm "
        "slope"
    );

    ASSERT_(
        min_probability >= 0.0 &&
        min_probability <= 1.0,
        "Minimum detachment probability must be between 0 and 1"
    );

    ASSERT_(
        max_probability >= min_probability &&
        max_probability <= 1.0,
        "Maximum detachment probability must be between the minimum "
        "probability and 1"
    );

    ASSERT_(
        k_ecm_weight >= 0.0 &&
        k_ecm_weight <= 1.0,
        "Detachment k_ecm weight must be between 0 and 1"
    );

    // -------------------------------------------------------------------------
    // Step 2: Select the slope on the appropriate side of the optimum
    // -------------------------------------------------------------------------

    const double selected_slope =
        k_ecm <= optimum_k_ecm
            ? low_k_ecm_slope
            : high_k_ecm_slope;

    const double distance_from_optimum =
        std::abs(k_ecm - optimum_k_ecm);

    // -------------------------------------------------------------------------
    // Step 3: Calculate the bounded stiffness response
    // -------------------------------------------------------------------------

    const double stiffness_response =
        1.0 -
        std::exp(
            -selected_slope * distance_from_optimum
        );

    // -------------------------------------------------------------------------
    // Step 4: Calculate the reference probability
    // -------------------------------------------------------------------------

    const double reference_probability =
        min_probability +
        k_ecm_weight *
        (max_probability - min_probability) *
        stiffness_response;

    ASSERT_(
        reference_probability >= 0.0 &&
        reference_probability <= 1.0,
        "Calculated reference detachment probability is outside [0, 1]"
    );

    return reference_probability;
  }

  inline
  double CalculateTimeAdjustedDetachmentProbability(
      const double reference_probability,
      const double time_step,
      const double reference_detachment_time) const {
    /*
     * Function goal
     * -------------
     * Convert a detachment probability defined over a reference time interval
     * into the equivalent probability for the current ABM timestep.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    ASSERT_(
        reference_probability >= 0.0 &&
        reference_probability <= 1.0,
        "Reference detachment probability must be between 0 and 1"
    );

    ASSERT_(
        time_step > 0.0,
        "Time-adjusted detachment probability requires a positive timestep"
    );

    ASSERT_(
        reference_detachment_time > 0.0,
        "Time-adjusted detachment probability requires a positive reference "
        "time"
    );

    // -------------------------------------------------------------------------
    // Step 2: Scale the probability to the current timestep
    // -------------------------------------------------------------------------

    const double time_adjusted_probability =
        1.0 -
        std::pow(
            1.0 - reference_probability,
            time_step / reference_detachment_time
        );

    ASSERT_(
        time_adjusted_probability >= 0.0 &&
        time_adjusted_probability <= 1.0,
        "Calculated time-adjusted detachment probability is outside [0, 1]"
    );

    return time_adjusted_probability;
  }

  inline
  std::vector<double> CalculateAttachmentDetachmentProbabilities(
      const double optimum_k_ecm,
      const double low_k_ecm_slope,
      const double high_k_ecm_slope,
      const double min_probability,
      const double max_probability,
      const double k_ecm_weight,
      const double time_step,
      const double reference_detachment_time) const {
    /*
     * Function goal
     * -------------
     * Calculate the timestep-adjusted detachment probability for every current
     * attachment using its latest valid local scaffold stiffness (k_ecm).
     *
     * Returned probabilities follow the same order as attachment_records_.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the current cell-matrix state
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Attachment detachment probabilities require an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Attachment detachment probabilities require attachment records"
    );

    // A cell must always retain at least one attachment.
    if (attachment_records_.size() == 1) {
      return {0.0};
    }

    ASSERT_(
        !requires_mechanics_recalculation_,
        "Attachment detachment probabilities require current FEM mechanics"
    );

    // -------------------------------------------------------------------------
    // Step 2: Calculate one probability per attachment
    // -------------------------------------------------------------------------

    std::vector<double> detachment_probabilities;

    detachment_probabilities.reserve(
        attachment_records_.size()
    );

    for (const auto& attachment : attachment_records_) {
      ASSERT_(
          attachment.node_id > 0,
          "Attachment detachment probability encountered an invalid node ID"
      );

      ASSERT_(
          attachment.has_valid_k_ecm,
          "Attachment detachment probability requires valid k_ecm"
      );

      ASSERT_(
          std::isfinite(attachment.k_ecm),
          "Attachment detachment probability encountered non-finite k_ecm"
      );

      ASSERT_(
          attachment.k_ecm >= 0.0,
          "Attachment detachment probability encountered negative k_ecm"
      );

      const double reference_probability =
          this->CalculateReferenceDetachmentProbability(
              attachment.k_ecm,
              optimum_k_ecm,
              low_k_ecm_slope,
              high_k_ecm_slope,
              min_probability,
              max_probability,
              k_ecm_weight
          );

      const double time_adjusted_probability =
          this->CalculateTimeAdjustedDetachmentProbability(
              reference_probability,
              time_step,
              reference_detachment_time
          );

      detachment_probabilities.push_back(
          time_adjusted_probability
      );
    }

    ASSERT_(
        detachment_probabilities.size() ==
            attachment_records_.size(),
        "Attachment detachment probability count does not match attachment "
        "record count"
    );

    return detachment_probabilities;
  }

  inline
  std::vector<std::size_t> SelectAttachmentIndicesForDetachment(
      const std::vector<double>& detachment_probabilities,
      const std::vector<double>& random_draws) const {
    /*
     * Function goal
     * -------------
     * Select attachments for detachment using one random draw per attachment,
     * while ensuring that an established cell always retains at least one.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the current attachment state
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Attachment detachment selection requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Attachment detachment selection requires attachment records"
    );

    ASSERT_(
        detachment_probabilities.size() ==
            attachment_records_.size(),
        "Detachment probability count does not match attachment count"
    );

    ASSERT_(
        random_draws.size() ==
            attachment_records_.size(),
        "Detachment random-draw count does not match attachment count"
    );

    // A previously attached cell must never lose its final attachment.
    if (attachment_records_.size() == 1) {
      ASSERT_(
          detachment_probabilities.front() == 0.0,
          "A single attachment must have zero detachment probability"
      );

      return {};
    }

    // -------------------------------------------------------------------------
    // Step 2: Apply one stochastic draw to each attachment
    // -------------------------------------------------------------------------

    std::vector<std::size_t> selected_indices;

    selected_indices.reserve(
        attachment_records_.size() - 1
    );

    for (std::size_t i = 0;
         i < attachment_records_.size();
         ++i) {
      const double probability =
          detachment_probabilities[i];

      const double random_draw =
          random_draws[i];

      ASSERT_(
          std::isfinite(probability) &&
          probability >= 0.0 &&
          probability <= 1.0,
          "Attachment detachment selection encountered an invalid probability"
      );

      ASSERT_(
          std::isfinite(random_draw) &&
          random_draw >= 0.0 &&
          random_draw <= 1.0,
          "Attachment detachment selection encountered an invalid random draw"
      );

      // Use < so that an attachment with probability zero can never detach.
      if (random_draw < probability) {
        selected_indices.push_back(i);
      }
    }

    // -------------------------------------------------------------------------
    // Step 3: Preserve one attachment if every draw selected detachment
    // -------------------------------------------------------------------------

    if (selected_indices.size() ==
        attachment_records_.size()) {
      std::size_t retained_index =
          selected_indices.front();

      double smallest_selection_margin =
          detachment_probabilities[retained_index] -
          random_draws[retained_index];

      for (const std::size_t index : selected_indices) {
        const double selection_margin =
            detachment_probabilities[index] -
            random_draws[index];

        const bool has_smaller_margin =
            selection_margin < smallest_selection_margin;

        const bool equal_margin_lower_node_id =
            selection_margin == smallest_selection_margin &&
            attachment_records_[index].node_id <
                attachment_records_[retained_index].node_id;

        if (has_smaller_margin ||
            equal_margin_lower_node_id) {
          retained_index = index;
          smallest_selection_margin =
              selection_margin;
        }
      }

      selected_indices.erase(
          std::remove(
              selected_indices.begin(),
              selected_indices.end(),
              retained_index
          ),
          selected_indices.end()
      );
    }

    ASSERT_(
        selected_indices.size() <
            attachment_records_.size(),
        "Attachment detachment selection attempted to remove every attachment"
    );

    return selected_indices;
  }

  inline
  bool ApplySelectedAttachmentDetachments(
      const std::vector<std::size_t>& selected_indices) {
    /*
     * Function goal
     * -------------
     * Remove the attachments selected for detachment while preserving at least
     * one retained attachment and marking the cell mechanics as outdated.
     *
     * Returns true when the attachment set changes.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the current cell state
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Attachment detachment requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Attachment detachment requires attachment records"
    );

    // No selected attachments means that no update is required.
    if (selected_indices.empty()) {
      return false;
    }

    ASSERT_(
        selected_indices.size() <
            attachment_records_.size(),
        "Attachment detachment cannot remove every attachment"
    );

    // -------------------------------------------------------------------------
    // Step 2: Validate and record the selected indices
    // -------------------------------------------------------------------------

    std::unordered_set<std::size_t> selected_index_set;

    for (const std::size_t index : selected_indices) {
      ASSERT_(
          index < attachment_records_.size(),
          "Attachment detachment encountered an out-of-range index"
      );

      const bool inserted =
          selected_index_set.insert(index).second;

      ASSERT_(
          inserted,
          "Attachment detachment received a duplicate selected index"
      );
    }

    // -------------------------------------------------------------------------
    // Step 3: Preserve the retained attachment records
    // -------------------------------------------------------------------------

    std::vector<AttachmentRecord> retained_attachments;

    retained_attachments.reserve(
        attachment_records_.size() -
        selected_indices.size()
    );

    for (std::size_t index = 0;
         index < attachment_records_.size();
         ++index) {
      if (selected_index_set.find(index) !=
          selected_index_set.end()) {
        continue;
      }

      retained_attachments.push_back(
          attachment_records_[index]
      );
    }

    ASSERT_(
        !retained_attachments.empty(),
        "Attachment detachment produced an empty retained attachment set"
    );

    ASSERT_(
        retained_attachments.size() +
            selected_indices.size() ==
            attachment_records_.size(),
        "Attachment detachment produced an inconsistent attachment count"
    );

    // -------------------------------------------------------------------------
    // Step 4: Update the attachment set
    // -------------------------------------------------------------------------

    this->SetAttachmentRecords(
        retained_attachments
    );

    this->MarkMechanicsForRecalculation();

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Attachment detachment unexpectedly changed the cell lifecycle"
    );

    ASSERT_(
        requires_mechanics_recalculation_,
        "Attachment detachment must request mechanics recalculation"
    );

    return true;
  }

  inline
  bool AttemptAttachmentDetachment(
      const double optimum_k_ecm,
      const double low_k_ecm_slope,
      const double high_k_ecm_slope,
      const double min_probability,
      const double max_probability,
      const double k_ecm_weight,
      const double time_step,
      const double reference_detachment_time) {
    /*
     * Function goal
     * -------------
     * Calculate attachment-level detachment probabilities, perform one random
     * draw per attachment and remove the selected attachments.
     *
     * Returns true when at least one attachment is removed.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the current attachment state
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Stochastic attachment detachment requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Stochastic attachment detachment requires attachment records"
    );

    // The final retained attachment cannot detach.
    if (attachment_records_.size() == 1) {
      return false;
    }

    // -------------------------------------------------------------------------
    // Step 2: Calculate one probability per attachment
    // -------------------------------------------------------------------------

    const std::vector<double> detachment_probabilities =
        this->CalculateAttachmentDetachmentProbabilities(
            optimum_k_ecm,
            low_k_ecm_slope,
            high_k_ecm_slope,
            min_probability,
            max_probability,
            k_ecm_weight,
            time_step,
            reference_detachment_time
        );

    ASSERT_(
        detachment_probabilities.size() ==
            attachment_records_.size(),
        "Stochastic detachment probability count does not match attachment "
        "count"
    );

    // -------------------------------------------------------------------------
    // Step 3: Generate one random draw per attachment
    // -------------------------------------------------------------------------

    auto* simulation =
        bdm::Simulation::GetActive();

    ASSERT_(
        simulation != nullptr,
        "Stochastic attachment detachment requires an active simulation"
    );

    auto* random_generator =
        simulation->GetRandom();

    ASSERT_(
        random_generator != nullptr,
        "Stochastic attachment detachment requires a random-number generator"
    );

    std::vector<double> random_draws;

    random_draws.reserve(
        attachment_records_.size()
    );

    for (std::size_t index = 0;
         index < attachment_records_.size();
         ++index) {
      random_draws.push_back(
          random_generator->Uniform(0.0, 1.0)
      );
    }

    // -------------------------------------------------------------------------
    // Step 4: Select the attachments whose draws permit detachment
    // -------------------------------------------------------------------------

    const std::vector<std::size_t> selected_indices =
        this->SelectAttachmentIndicesForDetachment(
            detachment_probabilities,
            random_draws
        );

    // -------------------------------------------------------------------------
    // Step 5: Apply the selected detachments
    // -------------------------------------------------------------------------

    return this->ApplySelectedAttachmentDetachments(
        selected_indices
    );
  }

  inline
  bool ProcessAttachmentDetachment() {
    /*
     * Function goal
     * -------------
     * Confirm that the cell is eligible for attachment turnover, read the
     * phenotype-specific detachment parameters and attempt stochastic
     * detachment.
     *
     * Returns true when at least one attachment is removed.
     */

    // -------------------------------------------------------------------------
    // Step 1: Exclude cells that cannot undergo attachment turnover
    // -------------------------------------------------------------------------

    if (this->GetPhenotype() < 1) {
      return false;
    }

    if (cell_matrix_lifecycle_status_ !=
        CellMatrixLifecycleStatus::kEstablished) {
      return false;
    }

    if (attachment_records_.size() <= 1) {
      return false;
    }

    // Do not use attachment mechanics that are waiting for FEM recalculation.
    if (requires_mechanics_recalculation_) {
      return false;
    }

    // -------------------------------------------------------------------------
    // Step 2: Read the phenotype mechanics configuration
    // -------------------------------------------------------------------------

    const std::string& phenotype_name =
        this->params()->get<std::string>(
            "phenotype_ID/" +
            std::to_string(this->GetPhenotype())
        );

    const std::string mech_base =
        phenotype_name + "/cell_matrix_mechanics";

    const bool mechanics_enabled =
        this->params()->have_parameter<bool>(
            mech_base + "/enabled"
        ) &&
        this->params()->get<bool>(
            mech_base + "/enabled"
        );

    if (!mechanics_enabled) {
      return false;
    }

    // -------------------------------------------------------------------------
    // Step 3: Confirm that current attachment stiffness values are usable
    // -------------------------------------------------------------------------

    for (const auto& attachment : attachment_records_) {
      ASSERT_(
          attachment.node_id > 0,
          "Attachment detachment encountered a non-positive node ID"
      );

      if (!attachment.has_valid_k_ecm) {
        return false;
      }

      ASSERT_(
          std::isfinite(attachment.k_ecm) &&
          attachment.k_ecm >= 0.0,
          "Attachment detachment encountered invalid k_ecm"
      );
    }

    // -------------------------------------------------------------------------
    // Step 4: Read the detachment parameters
    // -------------------------------------------------------------------------

    const double optimum_k_ecm =
        this->params()->get<double>(
            mech_base + "/detachment_optimum_kecm"
        );

    const double low_k_ecm_slope =
        this->params()->get<double>(
            mech_base + "/detachment_low_kecm_slope"
        );

    const double high_k_ecm_slope =
        this->params()->get<double>(
            mech_base + "/detachment_high_kecm_slope"
        );

    const double min_probability =
        this->params()->get<double>(
            mech_base + "/detachment_min_probability"
        );

    const double max_probability =
        this->params()->get<double>(
            mech_base + "/detachment_max_probability"
        );

    const double k_ecm_weight =
        this->params()->get<double>(
            mech_base + "/detachment_kecm_weight"
        );

    const double reference_detachment_time =
        this->params()->get<double>(
            mech_base + "/reference_detachment_time"
        );

    const double time_step =
        this->params()->get<double>(
            "time_step"
        );

    // -------------------------------------------------------------------------
    // Step 5: Attempt stochastic detachment
    // -------------------------------------------------------------------------

    return this->AttemptAttachmentDetachment(
        optimum_k_ecm,
        low_k_ecm_slope,
        high_k_ecm_slope,
        min_probability,
        max_probability,
        k_ecm_weight,
        time_step,
        reference_detachment_time
    );
  }

  inline
  double CalculateTimeAdjustedAttachmentFormationMean(
      const double mean_additions_per_reference_time,
      const double time_step,
      const double reference_attachment_formation_time) const {
    /*
     * Function goal
     * -------------
     * Scale the expected number of attachment additions from its configured
     * reference interval to the current ABM timestep.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    ASSERT_(
        mean_additions_per_reference_time >= 0.0,
        "Mean attachment additions per reference time must be non-negative"
    );

    ASSERT_(
        time_step > 0.0,
        "Attachment formation requires a positive timestep"
    );

    ASSERT_(
        reference_attachment_formation_time > 0.0,
        "Attachment formation requires a positive reference time"
    );

    // -------------------------------------------------------------------------
    // Step 2: Scale the mean to the current timestep
    // -------------------------------------------------------------------------

    const double time_adjusted_mean =
        mean_additions_per_reference_time *
        time_step /
        reference_attachment_formation_time;

    ASSERT_(
        std::isfinite(time_adjusted_mean) &&
        time_adjusted_mean >= 0.0,
        "Calculated attachment formation mean is invalid"
    );

    return time_adjusted_mean;
  }

  inline
  std::size_t SampleAttachmentAdditionCount(
      const double time_adjusted_mean,
      const std::size_t available_attachment_slots) const {
    /*
     * Function goal
     * -------------
     * Sample the number of attachment additions from a Poisson distribution
     * and limit the result to the number of available attachment slots.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    ASSERT_(
        std::isfinite(time_adjusted_mean) &&
        time_adjusted_mean >= 0.0,
        "Attachment addition sampling requires a valid non-negative mean"
    );

    // A cell at its maximum attachment count cannot form more attachments.
    if (available_attachment_slots == 0) {
      return 0;
    }

    // A zero mean disables attachment formation without requiring a draw.
    if (time_adjusted_mean == 0.0) {
      return 0;
    }

    // -------------------------------------------------------------------------
    // Step 2: Access BioDynaMo's random-number generator
    // -------------------------------------------------------------------------

    auto* simulation =
        bdm::Simulation::GetActive();

    ASSERT_(
        simulation != nullptr,
        "Attachment addition sampling requires an active simulation"
    );

    auto* random_generator =
        simulation->GetRandom();

    ASSERT_(
        random_generator != nullptr,
        "Attachment addition sampling requires a random-number generator"
    );

    // -------------------------------------------------------------------------
    // Step 3: Draw the requested number of attachment additions
    // -------------------------------------------------------------------------

    const int sampled_addition_count =
        random_generator->Poisson(
            time_adjusted_mean
        );

    ASSERT_(
        sampled_addition_count >= 0,
        "Poisson sampling returned a negative attachment addition count"
    );

    // -------------------------------------------------------------------------
    // Step 4: Cap the count by the available attachment slots
    // -------------------------------------------------------------------------

    const std::size_t capped_addition_count =
        std::min(
            static_cast<std::size_t>(
                sampled_addition_count
            ),
            available_attachment_slots
        );

    ASSERT_(
        capped_addition_count <= available_attachment_slots,
        "Attachment addition count exceeds the available attachment slots"
    );

    return capped_addition_count;
  }

  inline
  std::size_t DetermineAttachmentAdditionCount() const {
    /*
     * Function goal
     * -------------
     * Determine how many new attachments the cell should attempt to form
     * during the current timestep, limited by its available attachment slots.
     */

    // -------------------------------------------------------------------------
    // Step 1: Confirm that attachment formation is applicable
    // -------------------------------------------------------------------------

    if (this->GetPhenotype() < 1) {
      return 0;
    }

    if (cell_matrix_lifecycle_status_ !=
        CellMatrixLifecycleStatus::kEstablished) {
      return 0;
    }

    ASSERT_(
        !attachment_records_.empty(),
        "Attachment addition count requires at least one retained attachment"
    );

    // -------------------------------------------------------------------------
    // Step 2: Read the phenotype mechanics configuration
    // -------------------------------------------------------------------------

    const std::string& phenotype_name =
        this->params()->get<std::string>(
            "phenotype_ID/" +
            std::to_string(this->GetPhenotype())
        );

    const std::string mech_base =
        phenotype_name + "/cell_matrix_mechanics";

    const bool mechanics_enabled =
        this->params()->have_parameter<bool>(
            mech_base + "/enabled"
        ) &&
        this->params()->get<bool>(
            mech_base + "/enabled"
        );

    if (!mechanics_enabled) {
      return 0;
    }

    const int maximum_attachment_count =
        this->params()->get<int>(
            mech_base + "/num_attachments"
        );

    ASSERT_(
        maximum_attachment_count > 0,
        "Maximum attachment count must be positive"
    );

    ASSERT_(
        attachment_records_.size() <=
            static_cast<std::size_t>(
                maximum_attachment_count
            ),
        "Current attachment count exceeds the configured maximum"
    );

    // -------------------------------------------------------------------------
    // Step 3: Calculate the number of available attachment slots
    // -------------------------------------------------------------------------

    const std::size_t available_attachment_slots =
        static_cast<std::size_t>(
            maximum_attachment_count
        ) -
        attachment_records_.size();

    // Avoid probability calculations and random sampling when the cell is full.
    if (available_attachment_slots == 0) {
      return 0;
    }

    // -------------------------------------------------------------------------
    // Step 4: Calculate the timestep-adjusted Poisson mean
    // -------------------------------------------------------------------------

    const double mean_additions_per_reference_time =
        this->params()->get<double>(
            mech_base +
            "/mean_attachment_additions_per_reference_time"
        );

    const double reference_attachment_formation_time =
        this->params()->get<double>(
            mech_base +
            "/reference_attachment_formation_time"
        );

    const double time_step =
        this->params()->get<double>(
            "time_step"
        );

    const double time_adjusted_mean =
        this->CalculateTimeAdjustedAttachmentFormationMean(
            mean_additions_per_reference_time,
            time_step,
            reference_attachment_formation_time
        );

    // -------------------------------------------------------------------------
    // Step 5: Sample and cap the requested addition count
    // -------------------------------------------------------------------------

    return this->SampleAttachmentAdditionCount(
        time_adjusted_mean,
        available_attachment_slots
    );
  }

  inline
  std::vector<int> GenerateAttachmentCandidateNodeIds(
      const ObstacleScaffold& active_scaffold,
      const std::size_t requested_addition_count,
      const double min_cell_reach_radius,
      const double max_cell_reach_radius) const {
    /*
     * Function goal
     * -------------
     * Return scaffold node IDs that satisfy the pairwise attachment reach
     * constraints for the cell's current retained attachment set.
     *
     * Candidate generation is skipped when the cell has not requested any
     * attachment additions.
     */

    // -------------------------------------------------------------------------
    // Step 1: Skip cells that will not form attachments
    // -------------------------------------------------------------------------

    if (requested_addition_count == 0) {
      return {};
    }

    // -------------------------------------------------------------------------
    // Step 2: Validate the current cell and scaffold state
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Attachment candidate generation requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Attachment candidate generation requires retained attachments"
    );

    ASSERT_(
        min_cell_reach_radius >= 0.0,
        "Minimum cell reach radius must be non-negative"
    );

    ASSERT_(
        max_cell_reach_radius > 0.0,
        "Maximum cell reach radius must be positive"
    );

    ASSERT_(
        min_cell_reach_radius <=
            2.0 * max_cell_reach_radius,
        "Minimum attachment distance exceeds the maximum permitted distance"
    );

    const double maximum_pairwise_distance =
        2.0 * max_cell_reach_radius;

    // -------------------------------------------------------------------------
    // Step 3: Select a deterministic retained attachment as the search centre
    // -------------------------------------------------------------------------

    const auto search_attachment_it =
        std::min_element(
            attachment_records_.begin(),
            attachment_records_.end(),
            [](
                const AttachmentRecord& a,
                const AttachmentRecord& b) {
              return a.node_id < b.node_id;
            }
        );

    ASSERT_(
        search_attachment_it != attachment_records_.end(),
        "Attachment candidate generation could not select a search centre"
    );

    ASSERT_(
        active_scaffold.HasNode(
            search_attachment_it->node_id
        ),
        "Attachment candidate generation could not find retained node ID "
        + std::to_string(search_attachment_it->node_id)
    );

    const bdm::Double3& search_centre =
        active_scaffold.GetNodePosition(
            search_attachment_it->node_id
        );

    // -------------------------------------------------------------------------
    // Step 4: Query only the local scaffold region
    // -------------------------------------------------------------------------

    const std::vector<int> nearby_node_ids =
        active_scaffold.GetNodeIdsWithinRadius(
            search_centre,
            maximum_pairwise_distance
        );

    // -------------------------------------------------------------------------
    // Step 5: Apply exact pairwise reach checks
    // -------------------------------------------------------------------------

    std::vector<int> candidate_node_ids;

    for (const int candidate_node_id : nearby_node_ids) {
      ASSERT_(
          candidate_node_id > 0,
          "Attachment candidate generation encountered a non-positive node ID"
      );

      // Existing attachments cannot be selected again.
      bool already_attached = false;

      for (const auto& attachment : attachment_records_) {
        if (attachment.node_id == candidate_node_id) {
          already_attached = true;
          break;
        }
      }

      if (already_attached) {
        continue;
      }

      // Isolated nodes cannot form valid scaffold attachments.
      if (active_scaffold
              .GetConnectedNodeIds(candidate_node_id)
              .empty()) {
        continue;
      }

      const bdm::Double3& candidate_position =
          active_scaffold.GetNodePosition(
              candidate_node_id
          );

      bool satisfies_reach_constraints = true;

      for (const auto& retained_attachment :
           attachment_records_) {
        const bdm::Double3 difference =
            candidate_position -
            retained_attachment.position;

        const double distance =
            L2norm(difference);

        if (distance < min_cell_reach_radius ||
            distance > maximum_pairwise_distance) {
          satisfies_reach_constraints = false;
          break;
        }
      }

      if (!satisfies_reach_constraints) {
        continue;
      }

      candidate_node_ids.push_back(
          candidate_node_id
      );
    }

    // Preserve deterministic processing and testing.
    std::sort(
        candidate_node_ids.begin(),
        candidate_node_ids.end()
    );

    return candidate_node_ids;
  }

  inline
  double CalculateCandidateProximityWeight(
      const bdm::Double3& candidate_position,
      const double max_cell_reach_radius,
      const double proximity_sensitivity) const {
    /*
     * Function goal
     * -------------
     * Calculate a relative candidate weight based on the distance to the
     * nearest retained attachment.
     *
     * Candidates closer to a retained attachment receive a greater weight.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Candidate proximity weighting requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Candidate proximity weighting requires retained attachments"
    );

    ASSERT_(
        max_cell_reach_radius > 0.0,
        "Candidate proximity weighting requires a positive maximum reach"
    );

    ASSERT_(
        proximity_sensitivity >= 0.0,
        "Candidate proximity sensitivity must be non-negative"
    );

    ASSERT_(
        std::isfinite(candidate_position[0]) &&
        std::isfinite(candidate_position[1]) &&
        std::isfinite(candidate_position[2]),
        "Candidate proximity weighting received a non-finite position"
    );

    // -------------------------------------------------------------------------
    // Step 2: Find the nearest retained attachment
    // -------------------------------------------------------------------------

    const double maximum_pairwise_distance =
        2.0 * max_cell_reach_radius;

    double nearest_attachment_distance =
        maximum_pairwise_distance;

    bool found_retained_attachment = false;

    for (const auto& attachment : attachment_records_) {
      const bdm::Double3 difference =
          candidate_position -
          attachment.position;

      const double distance =
          L2norm(difference);

      ASSERT_(
          std::isfinite(distance) &&
          distance >= 0.0,
          "Candidate proximity weighting calculated an invalid distance"
      );

      if (!found_retained_attachment ||
          distance < nearest_attachment_distance) {
        nearest_attachment_distance = distance;
        found_retained_attachment = true;
      }
    }

    ASSERT_(
        found_retained_attachment,
        "Candidate proximity weighting could not find a retained attachment"
    );

    // -------------------------------------------------------------------------
    // Step 3: Normalise the distance
    // -------------------------------------------------------------------------

    const double normalised_distance =
        std::min(
            nearest_attachment_distance /
                maximum_pairwise_distance,
            1.0
        );

    // -------------------------------------------------------------------------
    // Step 4: Calculate the bounded relative weight
    // -------------------------------------------------------------------------

    const double proximity_weight =
        std::exp(
            -proximity_sensitivity *
            normalised_distance
        );

    ASSERT_(
        std::isfinite(proximity_weight) &&
        proximity_weight > 0.0 &&
        proximity_weight <= 1.0,
        "Calculated candidate proximity weight is outside (0, 1]"
    );

    return proximity_weight;
  }

  inline
  std::vector<double> CalculateAttachmentCandidateWeights(
      const ObstacleScaffold& active_scaffold,
      const std::vector<int>& candidate_node_ids,
      const double max_cell_reach_radius,
      const double proximity_sensitivity) const {
    /*
     * Function goal
     * -------------
     * Calculate one proximity-based selection weight for every geometrically
     * valid candidate attachment node.
     *
     * Returned weights follow the same order as candidate_node_ids.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the candidate list
    // -------------------------------------------------------------------------

    if (candidate_node_ids.empty()) {
      return {};
    }

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Candidate attachment weighting requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Candidate attachment weighting requires retained attachments"
    );

    ASSERT_(
        max_cell_reach_radius > 0.0,
        "Candidate attachment weighting requires a positive maximum reach"
    );

    ASSERT_(
        proximity_sensitivity >= 0.0,
        "Candidate proximity sensitivity must be non-negative"
    );

    // -------------------------------------------------------------------------
    // Step 2: Calculate one weight per candidate
    // -------------------------------------------------------------------------

    std::vector<double> candidate_weights;

    candidate_weights.reserve(
        candidate_node_ids.size()
    );

    std::unordered_set<int> encountered_node_ids;

    for (const int candidate_node_id : candidate_node_ids) {
      ASSERT_(
          candidate_node_id > 0,
          "Candidate attachment weighting encountered a non-positive node ID"
      );

      ASSERT_(
          active_scaffold.HasNode(candidate_node_id),
          "Candidate attachment weighting could not find scaffold node ID "
          + std::to_string(candidate_node_id)
      );

      const bool inserted =
          encountered_node_ids.insert(
              candidate_node_id
          ).second;

      ASSERT_(
          inserted,
          "Candidate attachment weighting received a duplicate node ID"
      );

      const bdm::Double3& candidate_position =
          active_scaffold.GetNodePosition(
              candidate_node_id
          );

      const double candidate_weight =
          this->CalculateCandidateProximityWeight(
              candidate_position,
              max_cell_reach_radius,
              proximity_sensitivity
          );

      candidate_weights.push_back(
          candidate_weight
      );
    }

    // -------------------------------------------------------------------------
    // Step 3: Confirm candidate-weight alignment
    // -------------------------------------------------------------------------

    ASSERT_(
        candidate_weights.size() ==
            candidate_node_ids.size(),
        "Candidate weight count does not match candidate node count"
    );

    return candidate_weights;
  }

  inline
  int SelectAttachmentCandidateNodeId(
      const std::vector<int>& candidate_node_ids,
      const std::vector<double>& candidate_weights) const {
    /*
     * Function goal
     * -------------
     * Select one candidate scaffold node using proximity-based relative
     * weights.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the candidate data
    // -------------------------------------------------------------------------

    ASSERT_(
        !candidate_node_ids.empty(),
        "Weighted candidate selection requires at least one candidate"
    );

    ASSERT_(
        candidate_weights.size() ==
            candidate_node_ids.size(),
        "Candidate node and weight counts do not match"
    );

    double total_weight = 0.0;

    for (std::size_t index = 0;
         index < candidate_node_ids.size();
         ++index) {
      ASSERT_(
          candidate_node_ids[index] > 0,
          "Weighted candidate selection encountered a non-positive node ID"
      );

      ASSERT_(
          std::isfinite(candidate_weights[index]) &&
          candidate_weights[index] > 0.0,
          "Weighted candidate selection encountered an invalid weight"
      );

      total_weight +=
          candidate_weights[index];
    }

    ASSERT_(
        std::isfinite(total_weight) &&
        total_weight > 0.0,
        "Weighted candidate selection requires a positive total weight"
    );

    // -------------------------------------------------------------------------
    // Step 2: Generate a random position within the total weight
    // -------------------------------------------------------------------------

    auto* simulation =
        bdm::Simulation::GetActive();

    ASSERT_(
        simulation != nullptr,
        "Weighted candidate selection requires an active simulation"
    );

    auto* random_generator =
        simulation->GetRandom();

    ASSERT_(
        random_generator != nullptr,
        "Weighted candidate selection requires a random-number generator"
    );

    const double random_draw =
        random_generator->Uniform(
            0.0,
            total_weight
        );

    // -------------------------------------------------------------------------
    // Step 3: Select the candidate containing the random position
    // -------------------------------------------------------------------------

    double cumulative_weight = 0.0;

    for (std::size_t index = 0;
         index < candidate_node_ids.size();
         ++index) {
      cumulative_weight +=
          candidate_weights[index];

      if (random_draw < cumulative_weight) {
        return candidate_node_ids[index];
      }
    }

    /*
     * Floating-point rounding may place a draw at the final cumulative
     * boundary. In that case, select the final candidate.
     */
    return candidate_node_ids.back();
  }

  inline
  void AddNewAttachmentRecord(
      const ObstacleScaffold& active_scaffold,
      const int selected_node_id) {
    /*
     * Function goal
     * -------------
     * Add one selected scaffold node to the cell's authoritative attachment
     * records and mark its mechanics for FEM recalculation.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the cell and selected scaffold node
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "New attachment formation requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "New attachment formation requires retained attachments"
    );

    ASSERT_(
        selected_node_id > 0,
        "New attachment formation received a non-positive scaffold node ID"
    );

    ASSERT_(
        active_scaffold.HasNode(selected_node_id),
        "New attachment formation could not find scaffold node ID "
        + std::to_string(selected_node_id)
    );

    ASSERT_(
        !active_scaffold
             .GetConnectedNodeIds(selected_node_id)
             .empty(),
        "New attachment formation cannot use an isolated scaffold node"
    );

    for (const auto& attachment : attachment_records_) {
      ASSERT_(
          attachment.node_id != selected_node_id,
          "New attachment formation attempted to add an existing attachment"
      );
    }

    const std::size_t previous_attachment_count =
        attachment_records_.size();

    // -------------------------------------------------------------------------
    // Step 2: Construct the new attachment record
    // -------------------------------------------------------------------------

    AttachmentRecord new_attachment;

    new_attachment.node_id =
        selected_node_id;

    new_attachment.position =
        active_scaffold.GetNodePosition(
            selected_node_id
        );

    /*
     * FEM has not yet calculated the local stiffness for this attachment.
     */
    new_attachment.k_ecm = 0.0;
    new_attachment.has_valid_k_ecm = false;
    new_attachment.newly_formed = true;

    // -------------------------------------------------------------------------
    // Step 3: Add the record to the authoritative attachment set
    // -------------------------------------------------------------------------

    std::vector<AttachmentRecord> updated_attachments =
        attachment_records_;

    updated_attachments.push_back(
        new_attachment
    );

    this->SetAttachmentRecords(
        updated_attachments
    );

    this->MarkMechanicsForRecalculation();

    // -------------------------------------------------------------------------
    // Step 4: Validate the updated state
    // -------------------------------------------------------------------------

    ASSERT_(
        attachment_records_.size() ==
            previous_attachment_count + 1,
        "New attachment formation did not increase the attachment count"
    );

    ASSERT_(
        attachment_records_.back().node_id ==
            selected_node_id,
        "New attachment formation stored an incorrect scaffold node ID"
    );

    ASSERT_(
        attachment_records_.back().newly_formed &&
        !attachment_records_.back().has_valid_k_ecm,
        "A newly formed attachment must await FEM mechanics"
    );

    ASSERT_(
        requires_mechanics_recalculation_,
        "New attachment formation must request mechanics recalculation"
    );

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "New attachment formation unexpectedly changed the cell lifecycle"
    );
  }

  inline
  std::size_t AttemptAttachmentFormation(
      const ObstacleScaffold& active_scaffold,
      const std::size_t requested_addition_count,
      const double min_cell_reach_radius,
      const double max_cell_reach_radius,
      const double proximity_sensitivity) {
    /*
     * Function goal
     * -------------
     * Attempt to form the requested number of new attachments using
     * geometrically valid candidates and proximity-weighted selection.
     *
     * The candidate pool is regenerated after every successful addition so
     * later candidates are checked against the updated attachment set.
     *
     * Returns the number of attachments successfully formed.
     */

    // -------------------------------------------------------------------------
    // Step 1: Skip cells that requested no additions
    // -------------------------------------------------------------------------

    if (requested_addition_count == 0) {
      return 0;
    }

    // -------------------------------------------------------------------------
    // Step 2: Validate the current state
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Attachment formation requires an established cell"
    );

    ASSERT_(
        !attachment_records_.empty(),
        "Attachment formation requires retained attachments"
    );

    ASSERT_(
        min_cell_reach_radius >= 0.0,
        "Attachment formation requires a non-negative minimum reach"
    );

    ASSERT_(
        max_cell_reach_radius > 0.0,
        "Attachment formation requires a positive maximum reach"
    );

    ASSERT_(
        proximity_sensitivity >= 0.0,
        "Attachment formation requires a non-negative proximity sensitivity"
    );

    const std::size_t original_attachment_count =
        attachment_records_.size();

    std::size_t formed_attachment_count = 0;

    // -------------------------------------------------------------------------
    // Step 3: Attempt each requested addition sequentially
    // -------------------------------------------------------------------------

    for (std::size_t addition = 0;
         addition < requested_addition_count;
         ++addition) {

      // -----------------------------------------------------------------------
      // Step 3a: Recalculate candidates using the current attachment set
      // -----------------------------------------------------------------------

      const std::vector<int> candidate_node_ids =
          this->GenerateAttachmentCandidateNodeIds(
              active_scaffold,
              1,
              min_cell_reach_radius,
              max_cell_reach_radius
          );

      /*
       * Fewer additions than requested may be formed when no valid candidates
       * remain.
       */
      if (candidate_node_ids.empty()) {
        break;
      }

      // -----------------------------------------------------------------------
      // Step 3b: Calculate the relative candidate weights
      // -----------------------------------------------------------------------

      const std::vector<double> candidate_weights =
          this->CalculateAttachmentCandidateWeights(
              active_scaffold,
              candidate_node_ids,
              max_cell_reach_radius,
              proximity_sensitivity
          );

      ASSERT_(
          candidate_weights.size() ==
              candidate_node_ids.size(),
          "Attachment formation candidate and weight counts do not match"
      );

      // -----------------------------------------------------------------------
      // Step 3c: Select and add one candidate
      // -----------------------------------------------------------------------

      const int selected_node_id =
          this->SelectAttachmentCandidateNodeId(
              candidate_node_ids,
              candidate_weights
          );

      this->AddNewAttachmentRecord(
          active_scaffold,
          selected_node_id
      );

      ++formed_attachment_count;
    }

    // -------------------------------------------------------------------------
    // Step 4: Validate the final attachment state
    // -------------------------------------------------------------------------

    ASSERT_(
        formed_attachment_count <= requested_addition_count,
        "Attachment formation exceeded the requested addition count"
    );

    ASSERT_(
        attachment_records_.size() ==
            original_attachment_count +
            formed_attachment_count,
        "Attachment formation produced an inconsistent attachment count"
    );

    if (formed_attachment_count > 0) {
      ASSERT_(
          requires_mechanics_recalculation_,
          "New attachment formation must request mechanics recalculation"
      );
    }

    return formed_attachment_count;
  }

};
// =============================================================================
} // ...end of namespace
// =============================================================================
#endif // _BIOLOGICAL_CELL_H_
// =============================================================================
