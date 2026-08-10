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

  struct ScaffoldOverlap {
    bool detected = false;

    int element_id = -1;

    bdm::Double3 closest_point = {
        0.0,
        0.0,
        0.0
    };

    double centreline_distance = 0.0;
    double required_separation = 0.0;
    double penetration_depth = 0.0;
  };

  bdm::Double3 SegmentClosestPoint(
    const ObstacleScaffold::Segment& segment,
    const bdm::Double3& point) const;

  ScaffoldOverlap FindDeepestOverlap(
    const ObstacleScaffold& scaffold,
    const bdm::Double3& proposed_position) const;

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
  double GetMinimumCellRadius() const;

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

  bdm::Double3 CalculatePreferredPosition(
    const ObstacleScaffold& scaffold) const;

bool ResolveScaffoldOverlap(
    const ObstacleScaffold& scaffold,
    const bdm::Double3& original_position,
    bdm::Double3* proposed_position) const;
  
  bool RepositionFromAttachments(
    const ObstacleScaffold& scaffold,
    const bdm::Double3& original_position);
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

  // Most recent non-zero cell displacement, used for directional persistence.
  bdm::Double3 last_migration_displacement_ = {0.0, 0.0, 0.0};
  
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
    double SegmentDistance(
        const bdm::Double3& segment_1_start,
        const bdm::Double3& segment_1_end,
        const bdm::Double3& segment_2_start,
        const bdm::Double3& segment_2_end) const {
    /*
    * Function goal
    * -------------
    * Calculate the shortest distance between two finite 3D line segments.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate coordinates
    // -------------------------------------------------------------------------

    for (std::size_t coordinate = 0;
        coordinate < 3;
        ++coordinate) {

        ASSERT_(
            std::isfinite(segment_1_start[coordinate]) &&
            std::isfinite(segment_1_end[coordinate]) &&
            std::isfinite(segment_2_start[coordinate]) &&
            std::isfinite(segment_2_end[coordinate]),
            "Segment-distance calculation received non-finite coordinates"
        );
    }

    const bdm::Double3 u =
        segment_1_end - segment_1_start;

    const bdm::Double3 v =
        segment_2_end - segment_2_start;

    const bdm::Double3 w =
        segment_1_start - segment_2_start;

    const double a = u * u;
    const double b = u * v;
    const double c = v * v;
    const double d = u * w;
    const double e = v * w;

    ASSERT_(
        a > 0.0,
        "Segment-distance calculation received a zero-length first segment"
    );

    ASSERT_(
        c > 0.0,
        "Segment-distance calculation received a zero-length second segment"
    );

    const double denominator =
        a * c - b * b;

    const double tolerance =
        1.0e-12;

    double s_numerator;
    double s_denominator =
        denominator;

    double t_numerator;
    double t_denominator =
        denominator;

    // -------------------------------------------------------------------------
    // Step 2: Find the closest positions along the infinite lines
    // -------------------------------------------------------------------------

    if (denominator < tolerance) {

        // Segments are almost parallel.
        s_numerator = 0.0;
        s_denominator = 1.0;

        t_numerator = e;
        t_denominator = c;

    } else {

        s_numerator =
            b * e - c * d;

        t_numerator =
            a * e - b * d;

        // Clamp the first segment parameter to [0, 1].
        if (s_numerator < 0.0) {
        s_numerator = 0.0;

        t_numerator = e;
        t_denominator = c;

        } else if (s_numerator >
                s_denominator) {

        s_numerator =
            s_denominator;

        t_numerator =
            e + b;

        t_denominator =
            c;
        }
    }

    // -------------------------------------------------------------------------
    // Step 3: Clamp the second segment parameter to [0, 1]
    // -------------------------------------------------------------------------

    if (t_numerator < 0.0) {

        t_numerator = 0.0;

        if (-d < 0.0) {
        s_numerator = 0.0;

        } else if (-d > a) {
        s_numerator =
            s_denominator;

        } else {
        s_numerator = -d;
        s_denominator = a;
        }

    } else if (t_numerator >
                t_denominator) {

        t_numerator =
            t_denominator;

        if ((-d + b) < 0.0) {
        s_numerator = 0.0;

        } else if ((-d + b) > a) {
        s_numerator =
            s_denominator;

        } else {
        s_numerator =
            -d + b;

        s_denominator =
            a;
        }
    }

    // -------------------------------------------------------------------------
    // Step 4: Calculate the closest-point separation
    // -------------------------------------------------------------------------

    const double s =
        std::fabs(s_numerator) < tolerance
            ? 0.0
            : s_numerator / s_denominator;

    const double t =
        std::fabs(t_numerator) < tolerance
            ? 0.0
            : t_numerator / t_denominator;

    const bdm::Double3 separation =
        w + u * s - v * t;

    const double distance =
        L2norm(
            separation
        );

    ASSERT_(
        std::isfinite(distance) &&
        distance >= 0.0,
        "Segment-distance calculation produced an invalid distance"
    );

    return distance;
    }

    inline
    double PointSegmentDistance(
        const bdm::Double3& point,
        const bdm::Double3& segment_start,
        const bdm::Double3& segment_end) const {
    /*
    * Function goal
    * -------------
    * Calculate the shortest distance between a point and a finite 3D segment.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the coordinates
    // -------------------------------------------------------------------------

    for (std::size_t coordinate = 0;
        coordinate < 3;
        ++coordinate) {

        ASSERT_(
            std::isfinite(point[coordinate]) &&
            std::isfinite(segment_start[coordinate]) &&
            std::isfinite(segment_end[coordinate]),
            "Point-segment distance received non-finite coordinates"
        );
    }

    // -------------------------------------------------------------------------
    // Step 2: Calculate the path projection
    // -------------------------------------------------------------------------

    const bdm::Double3 segment_vector =
        segment_end - segment_start;

    const double segment_length_squared =
        segment_vector * segment_vector;

    ASSERT_(
        segment_length_squared > 0.0,
        "Point-segment distance received a zero-length segment"
    );

    double projection =
        (
            (point - segment_start) *
            segment_vector
        ) /
        segment_length_squared;

    if (projection < 0.0) {
        projection = 0.0;
    }

    if (projection > 1.0) {
        projection = 1.0;
    }

    // -------------------------------------------------------------------------
    // Step 3: Calculate the shortest distance
    // -------------------------------------------------------------------------

    const bdm::Double3 closest_point =
        segment_start +
        segment_vector * projection;

    const double distance =
        L2norm(
            point - closest_point
        );

    ASSERT_(
        std::isfinite(distance) &&
        distance >= 0.0,
        "Point-segment distance produced an invalid distance"
    );

    return distance;
    }
    
    inline
    bool ScaffoldPathObstructed(
        const ObstacleScaffold& active_scaffold,
        int candidate_node_id,
        double additional_clearance) const {
    /*
    * Function goal
    * -------------
    * Determine whether an unrelated scaffold segment blocks the straight path
    * from the repositioned cell centre to a candidate attachment node.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the candidate and path
    // -------------------------------------------------------------------------

    ASSERT_(
        candidate_node_id > 0,
        "Scaffold-path obstruction check received a non-positive node ID"
    );

    ASSERT_(
        active_scaffold.HasNode(
            candidate_node_id
        ),
        "Scaffold-path obstruction check could not find candidate node ID "
        + std::to_string(candidate_node_id)
    );

    ASSERT_(
        std::isfinite(additional_clearance) &&
        additional_clearance >= 0.0,
        "Scaffold-path clearance must be finite and non-negative"
    );

    const bdm::Double3 path_start =
        this->GetPosition();

    const bdm::Double3& path_end =
        active_scaffold.GetNodePosition(
            candidate_node_id
        );

    const bdm::Double3 path_vector =
        path_end - path_start;

    const double path_length =
        L2norm(
            path_vector
        );

    ASSERT_(
        std::isfinite(path_length) &&
        path_length > 0.0,
        "Candidate attachment path must have positive length"
    );

    // -------------------------------------------------------------------------
    // Step 2: Define a broad-phase search around the complete candidate path
    // -------------------------------------------------------------------------

    const bdm::Double3 path_midpoint =
        (path_start + path_end) * 0.5;

    /*
    * A box centred on the path midpoint with this half-width contains the
    * complete path. The segment spatial index already accounts for strut radii.
    */
    const double search_radius =
        0.5 * path_length +
        additional_clearance;

    const std::vector<int> nearby_segment_ids =
        active_scaffold.GetNearbySegmentIds(
            path_midpoint,
            search_radius
        );

    // -------------------------------------------------------------------------
    // Step 3: Apply the exact obstruction test
    // -------------------------------------------------------------------------

    for (const int element_id :
        nearby_segment_ids) {

        ASSERT_(
            active_scaffold.HasSegment(
                element_id
            ),
            "Scaffold-path obstruction check could not find element ID "
            + std::to_string(element_id)
        );

        const ObstacleScaffold::Segment& scaffold_segment =
            active_scaffold.GetSegment(
                element_id
            );

        /*
        * Segments terminating at the candidate cannot obstruct the attachment
        * because the path is expected to finish on their centreline.
        */
        if (scaffold_segment.node_id_1 ==
                candidate_node_id ||
            scaffold_segment.node_id_2 ==
                candidate_node_id) {

        continue;
        }

        ASSERT_(
            std::isfinite(scaffold_segment.radius) &&
            scaffold_segment.radius >= 0.0,
            "Scaffold-path obstruction check encountered an invalid strut radius"
        );

        const double path_to_segment_distance =
            this->SegmentDistance(
                path_start,
                path_end,
                scaffold_segment.vertex_0,
                scaffold_segment.vertex_1
            );

        const double required_clearance =
            scaffold_segment.radius +
            additional_clearance;

        const double obstruction_tolerance =
            1.0e-12 *
            std::max(
                1.0,
                required_clearance
            );

        if (path_to_segment_distance <
            required_clearance - obstruction_tolerance) {
        return true;
        }
    }

    return false;
    }
    
    inline
    bool CellPathObstructed(
        const bdm::Double3& candidate_position,
        const double additional_clearance) const {
    /*
    * Function goal
    * -------------
    * Determine whether another BiologicalCell blocks the straight path from
    * this cell's repositioned centre to a candidate attachment node.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the candidate path
    // -------------------------------------------------------------------------

    ASSERT_(
        std::isfinite(additional_clearance) &&
        additional_clearance >= 0.0,
        "Candidate cell clearance must be finite and non-negative"
    );

    const bdm::Double3 path_start =
        this->GetPosition();

    const bdm::Double3 path_end =
        candidate_position;

    const double path_length =
        L2norm(
            path_end - path_start
        );

    ASSERT_(
        std::isfinite(path_length) &&
        path_length > 0.0,
        "Candidate cell-obstruction path must have positive length"
    );

    // -------------------------------------------------------------------------
    // Step 2: Access the current cell population
    // -------------------------------------------------------------------------

    auto* simulation =
        bdm::Simulation::GetActive();

    ASSERT_(
        simulation != nullptr,
        "Cell-path obstruction requires an active simulation"
    );

    auto* resource_manager =
        simulation->GetResourceManager();

    ASSERT_(
        resource_manager != nullptr,
        "Cell-path obstruction requires a resource manager"
    );

    // -------------------------------------------------------------------------
    // Step 3: Test every other BiologicalCell
    // -------------------------------------------------------------------------

    bool obstructed = false;
    std::mutex obstruction_mutex;

    resource_manager->ForEachAgent(
        [&](bdm::Agent* agent) {

            auto* other_cell =
                dynamic_cast<bdm::BiologicalCell*>(
                    agent
                );

            if (other_cell == nullptr ||
                other_cell == this) {
            return;
            }

            const bdm::Double3 other_position =
                other_cell->GetPosition();

            const double other_radius =
                0.5 *
                other_cell->GetDiameter();

            ASSERT_(
                std::isfinite(other_radius) &&
                other_radius > 0.0,
                "Cell-path obstruction encountered an invalid cell radius"
            );

            const double path_distance =
                this->PointSegmentDistance(
                    other_position,
                    path_start,
                    path_end
                );

            const double required_clearance =
                other_radius +
                additional_clearance;

            const double obstruction_tolerance =
                1.0e-12 *
                std::max(
                    1.0,
                    required_clearance
                );

            if (path_distance <
                required_clearance -
                obstruction_tolerance) {

                std::lock_guard<std::mutex> lock(
                    obstruction_mutex
                );

                obstructed = true;
            }
            
        }
    );

    return obstructed;
    }

    inline
    std::vector<int> GenerateAttachmentCandidateNodeIds(
        const ObstacleScaffold& active_scaffold,
        const double min_attachment_separation,
        const double max_attachment_separation,
        const double candidate_scaffold_clearance,
        const double candidate_cell_clearance) const {
    /*
    * Function goal
    * -------------
    * Return valid scaffold-node candidates for new cell attachments.
    *
    * Candidates must satisfy the configured pairwise separation limits relative
    * to every retained attachment and have an unobstructed path from the cell
    * centre through both the scaffold and neighbouring cells.
    *
    * One retained attachment is used as the spatial-query centre because every
    * valid candidate must lie within max_attachment_separation of every retained
    * attachment.
    */
    
    // -------------------------------------------------------------------------
    // Step 1: Validate the current state
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
        std::isfinite(min_attachment_separation) &&
        min_attachment_separation >= 0.0,
        "Minimum attachment separation must be finite and non-negative"
    );

    ASSERT_(
        std::isfinite(max_attachment_separation) &&
        max_attachment_separation > 0.0,
        "Maximum attachment separation must be finite and positive"
    );

    ASSERT_(
        min_attachment_separation <=
            max_attachment_separation,
        "Minimum attachment separation exceeds maximum attachment separation"
    );

    ASSERT_(
        std::isfinite(candidate_scaffold_clearance) &&
        candidate_scaffold_clearance >= 0.0,
        "Candidate scaffold clearance must be finite and non-negative"
    );

    ASSERT_(
        std::isfinite(candidate_cell_clearance) &&
        candidate_cell_clearance >= 0.0,
        "Candidate cell clearance must be finite and non-negative"
    );

    // -------------------------------------------------------------------------
    // Step 2: Select one deterministic retained attachment for broad-phase search
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
        search_attachment_it !=
            attachment_records_.end(),
        "Attachment candidate generation could not select a search attachment"
    );

    ASSERT_(
        active_scaffold.HasNode(
            search_attachment_it->node_id
        ),
        "Attachment candidate generation could not find retained node ID "
        + std::to_string(
            search_attachment_it->node_id
        )
    );

    const bdm::Double3& search_centre =
        active_scaffold.GetNodePosition(
            search_attachment_it->node_id
        );

    // -------------------------------------------------------------------------
    // Step 3: Retrieve scaffold nodes inside the maximum separation
    // -------------------------------------------------------------------------

    const std::vector<int> nearby_node_ids =
        active_scaffold.GetNodeIdsWithinRadius(
            search_centre,
            max_attachment_separation
        );

    // -------------------------------------------------------------------------
    // Step 4: Record existing attachment IDs
    // -------------------------------------------------------------------------

    std::unordered_set<int> attached_node_ids;

    for (const auto& attachment :
        attachment_records_) {

        ASSERT_(
            attachment.node_id > 0,
            "Attachment candidate generation encountered an invalid retained "
            "attachment node ID"
        );

        ASSERT_(
            active_scaffold.HasNode(
                attachment.node_id
            ),
            "Attachment candidate generation could not find retained scaffold "
            "node ID "
            + std::to_string(
                attachment.node_id
            )
        );

        attached_node_ids.insert(
            attachment.node_id
        );
    }

    // -------------------------------------------------------------------------
    // Step 5: Apply exact attachment-separation constraints
    // -------------------------------------------------------------------------

    std::vector<int> candidate_node_ids;

    for (const int candidate_node_id :
        nearby_node_ids) {

        ASSERT_(
            candidate_node_id > 0,
            "Attachment candidate generation encountered a non-positive node ID"
        );

        // Existing attachment nodes cannot be candidates.
        if (attached_node_ids.find(
                candidate_node_id
            ) != attached_node_ids.end()) {
        continue;
        }

        // Isolated scaffold nodes cannot form meaningful attachments.
        if (active_scaffold
                .GetConnectedNodeIds(
                    candidate_node_id
                )
                .empty()) {
        continue;
        }

        const bdm::Double3& candidate_position =
            active_scaffold.GetNodePosition(
                candidate_node_id
            );

        ASSERT_(
            std::isfinite(candidate_position[0]) &&
            std::isfinite(candidate_position[1]) &&
            std::isfinite(candidate_position[2]),
            "Attachment candidate generation encountered a non-finite candidate "
            "position"
        );

        bool valid_separation = true;

        for (const auto& retained_attachment :
            attachment_records_) {

        const bdm::Double3&
            retained_position =
                active_scaffold.GetNodePosition(
                    retained_attachment.node_id
                );

        const double separation =
            L2norm(
                candidate_position -
                retained_position
            );

        ASSERT_(
            std::isfinite(separation) &&
            separation >= 0.0,
            "Attachment candidate generation calculated an invalid separation"
        );

        if (separation <
                min_attachment_separation ||
            separation >
                max_attachment_separation) {

            valid_separation = false;
            break;
        }
        }

        if (!valid_separation) {
        continue;
        }

        // Reject candidates whose attachment path crosses another scaffold strut.
        if (this->ScaffoldPathObstructed(
                active_scaffold,
                candidate_node_id,
                candidate_scaffold_clearance)) {

        continue;
        }

        // Candidate must not lie behind another cell.
        if (this->CellPathObstructed(
                candidate_position,
                candidate_cell_clearance)) {

        continue;
        }

        candidate_node_ids.push_back(
            candidate_node_id
        );
    }

    // -------------------------------------------------------------------------
    // Step 6: Preserve deterministic candidate ordering
    // -------------------------------------------------------------------------

    std::sort(
        candidate_node_ids.begin(),
        candidate_node_ids.end()
    );

    return candidate_node_ids;
    }
    
    inline
    double CalculateCandidateDistanceWeight(
        const double candidate_distance,
        const double minimum_candidate_distance,
        const double maximum_candidate_distance,
        const double distance_sensitivity) const {
    /*
    * Function goal
    * -------------
    * Calculate the candidate weight associated with its distance from the
    * repositioned cell centre.
    *
    * Candidates closer to the cell centre receive greater weight.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    ASSERT_(
        std::isfinite(candidate_distance) &&
        candidate_distance >= 0.0,
        "Candidate distance weighting received an invalid distance"
    );

    ASSERT_(
        std::isfinite(minimum_candidate_distance) &&
        minimum_candidate_distance >= 0.0,
        "Candidate distance weighting received an invalid minimum distance"
    );

    ASSERT_(
        std::isfinite(maximum_candidate_distance) &&
        maximum_candidate_distance >=
            minimum_candidate_distance,
        "Candidate distance weighting received invalid distance limits"
    );

    ASSERT_(
        std::isfinite(distance_sensitivity) &&
        distance_sensitivity >= 0.0,
        "Candidate distance sensitivity must be finite and non-negative"
    );

    // -------------------------------------------------------------------------
    // Step 2: Normalise the distance within the candidate pool
    // -------------------------------------------------------------------------

    double normalised_distance = 0.0;

    const double distance_range =
        maximum_candidate_distance -
        minimum_candidate_distance;

    if (distance_range > 0.0) {
        normalised_distance =
            (
            candidate_distance -
            minimum_candidate_distance
            ) /
            distance_range;
    }

    ASSERT_(
        normalised_distance >= 0.0 &&
        normalised_distance <= 1.0,
        "Normalised candidate distance is outside [0, 1]"
    );

    // -------------------------------------------------------------------------
    // Step 3: Calculate the bounded distance weight
    // -------------------------------------------------------------------------

    const double distance_weight =
        std::exp(
            -distance_sensitivity *
            normalised_distance
        );

    ASSERT_(
        std::isfinite(distance_weight) &&
        distance_weight > 0.0 &&
        distance_weight <= 1.0,
        "Calculated candidate distance weight is outside (0, 1]"
    );

    return distance_weight;
    }
  
    inline
    double CalculateCandidatePersistenceWeight(
        const bdm::Double3& candidate_position,
        const bdm::Double3& recent_movement,
        const double persistence_sensitivity) const {
    /*
    * Function goal
    * -------------
    * Calculate a directional-persistence weight that favours attachment
    * candidates lying ahead of the cell's most recent movement.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    ASSERT_(
        std::isfinite(persistence_sensitivity) &&
        persistence_sensitivity >= 0.0,
        "Candidate directional-persistence sensitivity must be finite and "
        "non-negative"
    );

    for (std::size_t coordinate = 0;
        coordinate < 3;
        ++coordinate) {

        ASSERT_(
            std::isfinite(candidate_position[coordinate]) &&
            std::isfinite(recent_movement[coordinate]),
            "Candidate persistence weighting received non-finite coordinates"
        );
    }

    // -------------------------------------------------------------------------
    // Step 2: Determine whether a meaningful movement direction exists
    // -------------------------------------------------------------------------

    const double movement_magnitude =
        L2norm(
            recent_movement
        );

    const double movement_tolerance =
        this->params()->get<double>(
            "migration_tolerance"
        );

    /*
    * If the cell did not move far enough to define a reliable direction,
    * directional persistence is disabled for this candidate.
    */
    if (movement_magnitude <=
        movement_tolerance) {

        return 1.0;
    }

    // -------------------------------------------------------------------------
    // Step 3: Calculate the direction from the current cell centre to candidate
    // -------------------------------------------------------------------------

    const bdm::Double3 candidate_direction =
        candidate_position -
        this->GetPosition();

    const double candidate_distance =
        L2norm(
            candidate_direction
        );

    ASSERT_(
        std::isfinite(candidate_distance) &&
        candidate_distance > 0.0,
        "Candidate persistence weighting requires a positive candidate distance"
    );

    // -------------------------------------------------------------------------
    // Step 4: Calculate directional alignment
    // -------------------------------------------------------------------------

    double cosine_alignment =
        (
            recent_movement *
            candidate_direction
        ) /
        (
            movement_magnitude *
            candidate_distance
        );

    /*
    * Protect against small floating-point excursions beyond [-1, 1].
    */
    cosine_alignment =
        std::max(
            -1.0,
            std::min(
                1.0,
                cosine_alignment
            )
        );

    // -------------------------------------------------------------------------
    // Step 5: Convert alignment to a penalty
    // -------------------------------------------------------------------------

    const double directional_penalty =
        0.5 *
        (
            1.0 -
            cosine_alignment
        );

    ASSERT_(
        directional_penalty >= 0.0 &&
        directional_penalty <= 1.0,
        "Candidate directional penalty is outside [0, 1]"
    );

    // -------------------------------------------------------------------------
    // Step 6: Calculate the persistence weight
    // -------------------------------------------------------------------------

    const double persistence_weight =
        std::exp(
            -persistence_sensitivity *
            directional_penalty
        );

    ASSERT_(
        std::isfinite(persistence_weight) &&
        persistence_weight > 0.0 &&
        persistence_weight <= 1.0,
        "Calculated candidate persistence weight is outside (0, 1]"
    );

    return persistence_weight;
    }
    
    inline
    std::vector<double> CalculateAttachmentCandidateWeights(
        const ObstacleScaffold& active_scaffold,
        const std::vector<int>& candidate_node_ids,
        const bdm::Double3& recent_movement,
        const double distance_sensitivity,
        const double persistence_sensitivity) const {
    /*
    * Function goal
    * -------------
    * Calculate one combined weight for every valid attachment candidate.
    *
    * Each weight combines the candidate's distance from the repositioned cell
    * centre with its alignment to the cell's most recent movement direction.
    *
    * Returned weights follow the same order as candidate_node_ids.
    */

    if (candidate_node_ids.empty()) {
        return {};
    }

    ASSERT_(
        std::isfinite(distance_sensitivity) &&
        distance_sensitivity >= 0.0,
        "Candidate cell-distance sensitivity must be finite and non-negative"
    );

    ASSERT_(
        std::isfinite(persistence_sensitivity) &&
        persistence_sensitivity >= 0.0,
        "Candidate directional-persistence sensitivity must be finite and "
        "non-negative"
    );

    const bdm::Double3 cell_position =
        this->GetPosition();

    // Calculate candidate distances from the current cell centre.
    std::vector<double> candidate_distances;

    candidate_distances.reserve(
        candidate_node_ids.size()
    );

    for (const int candidate_node_id :
        candidate_node_ids) {

        ASSERT_(
            active_scaffold.HasNode(candidate_node_id),
            "Candidate weighting could not find scaffold node ID "
            + std::to_string(candidate_node_id)
        );

        const bdm::Double3& candidate_position =
            active_scaffold.GetNodePosition(
                candidate_node_id
            );

        const double candidate_distance =
            L2norm(
                candidate_position -
                cell_position
            );

        ASSERT_(
            std::isfinite(candidate_distance) &&
            candidate_distance >= 0.0,
            "Candidate weighting calculated an invalid cell distance"
        );

        candidate_distances.push_back(
            candidate_distance
        );
    }

    // Determine the distance range across the candidate pool.
    const auto distance_limits =
        std::minmax_element(
            candidate_distances.begin(),
            candidate_distances.end()
        );

    const double minimum_candidate_distance =
        *distance_limits.first;

    const double maximum_candidate_distance =
        *distance_limits.second;

    // Calculate the combined weight for every candidate.
    std::vector<double> candidate_weights;

    candidate_weights.reserve(
        candidate_node_ids.size()
    );

    for (std::size_t index = 0;
        index < candidate_node_ids.size();
        ++index) {

        const int candidate_node_id =
            candidate_node_ids[index];

        const bdm::Double3& candidate_position =
            active_scaffold.GetNodePosition(
                candidate_node_id
            );

        const double distance_weight =
            this->CalculateCandidateDistanceWeight(
                candidate_distances[index],
                minimum_candidate_distance,
                maximum_candidate_distance,
                distance_sensitivity
            );

        const double persistence_weight =
            this->CalculateCandidatePersistenceWeight(
                candidate_position,
                recent_movement,
                persistence_sensitivity
            );

        const double combined_weight =
            distance_weight *
            persistence_weight;

        ASSERT_(
            std::isfinite(combined_weight) &&
            combined_weight > 0.0 &&
            combined_weight <= 1.0,
            "Combined candidate weight is outside (0, 1]"
        );

        candidate_weights.push_back(
            combined_weight
        );
    }

    ASSERT_(
        candidate_weights.size() ==
            candidate_node_ids.size(),
        "Candidate weight count does not match candidate node count"
    );

    return candidate_weights;
    }

    inline
    std::vector<double> CalculateCandidateSelectionProbabilities(
        const std::vector<double>& candidate_weights,
        const double random_selection_strength) const {
    /*
    * Function goal
    * -------------
    * Convert candidate bias weights into attachment-selection probabilities.
    *
    * random_selection_strength controls the balance between biased and random
    * selection:
    *
    *   0.0 -> highest-weight candidate is selected deterministically
    *   0.5 -> probabilities are proportional to candidate weights
    *   1.0 -> all candidates have equal probability
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the inputs
    // -------------------------------------------------------------------------

    if (candidate_weights.empty()) {
        return {};
    }

    ASSERT_(
        std::isfinite(random_selection_strength) &&
        random_selection_strength >= 0.0 &&
        random_selection_strength <= 1.0,
        "Random selection strength must be finite and between 0 and 1"
    );

    for (const double weight : candidate_weights) {
        ASSERT_(
            std::isfinite(weight) &&
            weight > 0.0,
            "Candidate selection requires positive finite weights"
        );
    }

    std::vector<double> selection_probabilities(
        candidate_weights.size(),
        0.0
    );

    // -------------------------------------------------------------------------
    // Step 2: Handle fully deterministic selection
    // -------------------------------------------------------------------------

    if (random_selection_strength == 0.0) {

        const auto highest_weight =
            std::max_element(
                candidate_weights.begin(),
                candidate_weights.end()
            );

        const std::size_t selected_index =
            static_cast<std::size_t>(
                std::distance(
                    candidate_weights.begin(),
                    highest_weight
                )
            );

        selection_probabilities[selected_index] =
            1.0;

        return selection_probabilities;
    }

    // -------------------------------------------------------------------------
    // Step 3: Calculate the bias exponent
    // -------------------------------------------------------------------------

    const double gamma =
        (1.0 - random_selection_strength) /
        random_selection_strength;

    // -------------------------------------------------------------------------
    // Step 4: Calculate numerically stable relative scores
    // -------------------------------------------------------------------------

    std::vector<double> log_scores;

    log_scores.reserve(
        candidate_weights.size()
    );

    for (const double weight :
        candidate_weights) {

        log_scores.push_back(
            gamma *
            std::log(weight)
        );
    }

    const double maximum_log_score =
        *std::max_element(
            log_scores.begin(),
            log_scores.end()
        );

    double total_score = 0.0;

    for (std::size_t index = 0;
        index < log_scores.size();
        ++index) {

        const double score =
            std::exp(
                log_scores[index] -
                maximum_log_score
            );

        selection_probabilities[index] =
            score;

        total_score += score;
    }

    ASSERT_(
        std::isfinite(total_score) &&
        total_score > 0.0,
        "Candidate selection produced an invalid probability total"
    );

    // -------------------------------------------------------------------------
    // Step 5: Normalise the scores into probabilities
    // -------------------------------------------------------------------------

    double probability_sum = 0.0;

    for (double& probability :
        selection_probabilities) {

        probability /=
            total_score;

        ASSERT_(
            std::isfinite(probability) &&
            probability >= 0.0 &&
            probability <= 1.0,
            "Candidate selection probability is outside [0, 1]"
        );

        probability_sum += probability;
    }

    ASSERT_(
        std::fabs(probability_sum - 1.0) <= 1.0e-12,
        "Candidate selection probabilities do not sum to 1"
    );

    return selection_probabilities;
    }

    inline
    int SelectAttachmentCandidateNodeId(
        const std::vector<int>& candidate_node_ids,
        const std::vector<double>& selection_probabilities) const {
    /*
    * Function goal
    * -------------
    * Select one attachment candidate using the supplied selection probabilities.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the candidate pool
    // -------------------------------------------------------------------------

    ASSERT_(
        !candidate_node_ids.empty(),
        "Attachment candidate selection requires at least one candidate"
    );

    ASSERT_(
        candidate_node_ids.size() ==
            selection_probabilities.size(),
        "Candidate node and selection-probability counts do not match"
    );

    double probability_sum = 0.0;

    for (std::size_t index = 0;
        index < candidate_node_ids.size();
        ++index) {

        ASSERT_(
            candidate_node_ids[index] > 0,
            "Attachment candidate selection encountered a non-positive node ID"
        );

        ASSERT_(
            std::isfinite(selection_probabilities[index]) &&
            selection_probabilities[index] >= 0.0 &&
            selection_probabilities[index] <= 1.0,
            "Attachment candidate selection encountered an invalid probability"
        );

        probability_sum +=
            selection_probabilities[index];
    }

    ASSERT_(
        std::fabs(probability_sum - 1.0) <= 1.0e-12,
        "Attachment candidate selection probabilities do not sum to 1"
    );

    // -------------------------------------------------------------------------
    // Step 2: Handle deterministic selection
    // -------------------------------------------------------------------------

    for (std::size_t index = 0;
        index < selection_probabilities.size();
        ++index) {

        if (selection_probabilities[index] == 1.0) {
        return candidate_node_ids[index];
        }
    }

    // -------------------------------------------------------------------------
    // Step 3: Generate one random selection draw
    // -------------------------------------------------------------------------

    auto* simulation =
        bdm::Simulation::GetActive();

    ASSERT_(
        simulation != nullptr,
        "Attachment candidate selection requires an active simulation"
    );

    auto* random_generator =
        simulation->GetRandom();

    ASSERT_(
        random_generator != nullptr,
        "Attachment candidate selection requires a random-number generator"
    );

    const double random_draw =
        random_generator->Uniform(
            0.0,
            1.0
        );

    // -------------------------------------------------------------------------
    // Step 4: Select from the cumulative probability distribution
    // -------------------------------------------------------------------------

    double cumulative_probability = 0.0;

    for (std::size_t index = 0;
        index < candidate_node_ids.size();
        ++index) {

        cumulative_probability +=
            selection_probabilities[index];

        if (random_draw <
            cumulative_probability) {

        return candidate_node_ids[index];
        }
    }

    // Protect against floating-point accumulation leaving the final cumulative
    // probability fractionally below 1.
    return candidate_node_ids.back();
    }

    inline
    void AddSelectedAttachment(
        const ObstacleScaffold& active_scaffold,
        const int selected_node_id) {
    /*
    * Function goal
    * -------------
    * Add one ABM-selected scaffold node to the cell's persistent attachment
    * records and mark the cell mechanics for FEM recalculation.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the selected candidate
    // -------------------------------------------------------------------------

    ASSERT_(
        cell_matrix_lifecycle_status_ ==
            CellMatrixLifecycleStatus::kEstablished,
        "Attachment formation requires an established cell"
    );

    ASSERT_(
        selected_node_id > 0,
        "Attachment formation received a non-positive scaffold node ID"
    );

    ASSERT_(
        active_scaffold.HasNode(selected_node_id),
        "Attachment formation could not find scaffold node ID "
        + std::to_string(selected_node_id)
    );

    for (const auto& attachment :
        attachment_records_) {

        ASSERT_(
            attachment.node_id != selected_node_id,
            "Attachment formation attempted to add an already attached node"
        );
    }

    // -------------------------------------------------------------------------
    // Step 2: Create the new persistent attachment record
    // -------------------------------------------------------------------------

    AttachmentRecord new_attachment;

    new_attachment.node_id =
        selected_node_id;

    new_attachment.position =
        active_scaffold.GetNodePosition(
            selected_node_id
        );

    new_attachment.k_ecm =
        0.0;

    new_attachment.has_valid_k_ecm =
        false;

    new_attachment.newly_formed =
        true;

    // -------------------------------------------------------------------------
    // Step 3: Append the attachment to the existing persistent state
    // -------------------------------------------------------------------------

    std::vector<AttachmentRecord> updated_attachments =
        attachment_records_;

    updated_attachments.push_back(
        new_attachment
    );

    this->UpdateAttachmentRecordsFromAbm(
        updated_attachments
    );

    // -------------------------------------------------------------------------
    // Step 4: Validate the updated cell-matrix state
    // -------------------------------------------------------------------------

    ASSERT_(
        this->RequiresMechanicsRecalculation(),
        "New attachment formation must request FEM mechanics recalculation"
    );

    ASSERT_(
        this->GetNumberOfAttachmentRecords() ==
            updated_attachments.size(),
        "Attachment formation produced an inconsistent attachment count"
    );

    this->ValidateCellMatrixState();
    }

    inline
    std::size_t AttemptAttachmentFormation(
        const ObstacleScaffold& active_scaffold,
        const std::size_t requested_addition_count,
        const bdm::Double3& recent_movement,
        const double min_attachment_separation,
        const double max_attachment_separation,
        const double candidate_scaffold_clearance,
        const double candidate_cell_clearance,
        const double cell_distance_sensitivity,
        const double persistence_sensitivity,
        const double random_selection_strength) {
    /*
    * Function goal
    * -------------
    * Attempt the requested number of new attachments using the current
    * attachment geometry, candidate weights and stochastic selection rule.
    *
    * The candidate pool is regenerated after every successful addition.
    *
    * Returns the number of attachments actually formed.
    */

    // -------------------------------------------------------------------------
    // Step 1: Nothing requested
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

    const std::size_t original_attachment_count =
        attachment_records_.size();

    std::size_t formed_attachment_count = 0;

    // -------------------------------------------------------------------------
    // Step 3: Attempt each requested addition
    // -------------------------------------------------------------------------

    for (std::size_t addition = 0;
        addition < requested_addition_count;
        ++addition) {

        // Recalculate candidates because the attachment set may have changed.
        const std::vector<int> candidate_node_ids =
            this->GenerateAttachmentCandidateNodeIds(
                active_scaffold,
                min_attachment_separation,
                max_attachment_separation,
                candidate_scaffold_clearance,
                candidate_cell_clearance
            );

        // It is valid to form fewer attachments than requested.
        if (candidate_node_ids.empty()) {
        break;
        }

        // Calculate Stage 17 candidate weights.
        const std::vector<double> candidate_weights =
            this->CalculateAttachmentCandidateWeights(
                active_scaffold,
                candidate_node_ids,
                recent_movement,
                cell_distance_sensitivity,
                persistence_sensitivity
            );

        ASSERT_(
            candidate_weights.size() ==
                candidate_node_ids.size(),
            "Attachment formation candidate and weight counts do not match"
        );

        // Convert relative weights into Stage 18 selection probabilities.
        const std::vector<double> selection_probabilities =
            this->CalculateCandidateSelectionProbabilities(
                candidate_weights,
                random_selection_strength
            );

        ASSERT_(
            selection_probabilities.size() ==
                candidate_node_ids.size(),
            "Attachment formation candidate and probability counts do not match"
        );

        // Select and form one attachment.
        const int selected_node_id =
            this->SelectAttachmentCandidateNodeId(
                candidate_node_ids,
                selection_probabilities
            );

        this->AddSelectedAttachment(
            active_scaffold,
            selected_node_id
        );

        ++formed_attachment_count;
    }

    // -------------------------------------------------------------------------
    // Step 4: Validate the final attachment state
    // -------------------------------------------------------------------------

    ASSERT_(
        formed_attachment_count <=
            requested_addition_count,
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
            "New attachment formation must request FEM mechanics recalculation"
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
