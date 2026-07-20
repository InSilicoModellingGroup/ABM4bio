#ifndef CELL_MATRIX_INTERACTIONS_H_
#define CELL_MATRIX_INTERACTIONS_H_

// =============================================================================
// Project dependencies
// =============================================================================

#include "./global.h"
#include "./biological_cell.h"
#include "./obstacles.h"

// =============================================================================
// Forward declarations
// =============================================================================

/*
 * This function remains defined in ABM4bio.h for now.
 *
 * It can later be moved into the cell behaviour code because it represents a
 * biological contractility law rather than FEM communication logic.
 */
inline double calculate_contractile_force_from_kce(
    double k_ce,
    double min_contractile_force,
    double max_contractile_force,
    double contractility_kce_sensitivity,
    double contractility_kce_midpoint);


// =============================================================================
// Cell-matrix FEM–ABM interaction interface
// =============================================================================

class CellMatrixInteraction {
 public:
  
    // ---------------------------------------------------------------------------
    // FEM export action
    //
    // Defines how a BiologicalCell should be handled during the current
    // ABM-to-FEM export.
    // ---------------------------------------------------------------------------

    enum class FemExportAction {
    kSkip,
    kAttach,
    kContract
    };
 
    // ---------------------------------------------------------------------------
    // Construction
    // ---------------------------------------------------------------------------

    explicit CellMatrixInteraction(Parameters* model_params)
        : params_(model_params)
    {
        ASSERT_(
        params_ != nullptr,
        "CellMatrixInteraction requires a valid Parameters pointer"
        );
    }

    ~CellMatrixInteraction() = default;

    // Prevent construction without access to the model parameters.
    CellMatrixInteraction() = delete;

    // ---------------------------------------------------------------------------
    // Parameter access
    // ---------------------------------------------------------------------------

    inline Parameters* params() {
        return params_;
    }

    inline const Parameters* params() const {
        return params_;
    }

    // ---------------------------------------------------------------------------
    // Cell-matrix mechanics checks
    // ---------------------------------------------------------------------------

    inline
    bool AnyCellMatrixMechanicsEnabled(
        const std::map<int, std::string>& cells) const
        {
            /*
            * Check whether cell-matrix mechanics is enabled for any phenotype.
            */

            for (const auto& cell_type : cells) {
                const std::string& CP_name = cell_type.second;

                const std::string mech_base =
                CP_name + "/cell_matrix_mechanics";

                if (this->params()->have_parameter<bool>(mech_base + "/enabled") &&
                    this->params()->get<bool>(mech_base + "/enabled")) {
                return true;
                }
            }

        return false;
        }

    // ---------------------------------------------------------------------------
    // Attachment formatting
    // ---------------------------------------------------------------------------

    inline
    std::string FormatAttachmentNodeIdsJsonLike(
        const std::vector<int>& node_ids) const
    {
    
    /*
        * Format attachment node IDs as JSON-like list text.
        *
        * Function goal
        * -------------
        * Create a stable string representation that Python can parse robustly for FEM solver.
        *
        * Inputs
        * ------
        * node_ids : std::vector<int>
        *   Attachment node IDs stored on the BiologicalCell.
        *
        * Returns
        * -------
        * std::string
        *   JSON-like list string, for example "[5316, 33603]".
        *
        * Example:
        *
        *   [5316, 33603]
        */

    std::ostringstream oss;

    oss << "[";

    for (std::size_t i = 0; i < node_ids.size(); ++i) {
        oss << node_ids[i];

        if (i + 1 < node_ids.size()) {
        oss << ", ";
        }
    }

    oss << "]";

    return oss.str();
    }
    
    // ---------------------------------------------------------------------------
    // Determine FEM export action
    // ---------------------------------------------------------------------------

    inline
    FemExportAction DetermineFemExportAction(
        const bdm::BiologicalCell* cell,
        const std::map<int, std::string>& cells) const
    {
    /*
    * Function goal
    * -------------
    * Determine whether a cell should be omitted from the FEM input, sent for
    * initial attachment, or sent for contraction.
    *
    * Current rules
    * ---------------------
    * - Necrotic cells are omitted.
    * - Cells from phenotypes without cell-matrix mechanics are omitted.
    * - Cells with exactly one attachment remain entirely within the ABM.
    * - Cells with at least two valid attachments are sent as "contract",
    *   unless the existing movement flag requires LEGACY FEM reattachment.
    * - Cells with no attachments are sent as "attach".
    *
    * The explicit internal cell-matrix status to be introduced will later --------EDIT THIS COMMENT WHEN RELEVANT----------------------------------------------
    * replace some of these temporary attachment-count-based decisions.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the supplied cell pointer
    // -------------------------------------------------------------------------

    ASSERT_(
        cell != nullptr,
        "DetermineFemExportAction received a null BiologicalCell pointer"
    );

    // -------------------------------------------------------------------------
    // Step 2: Exclude necrotic cells
    // -------------------------------------------------------------------------

    const int phenotype_id =
        cell->GetPhenotype();

    if (phenotype_id < 0) {
        return FemExportAction::kSkip;
    }

    // -------------------------------------------------------------------------
    // Step 3: Find the cell phenotype definition
    // -------------------------------------------------------------------------

    const auto phenotype_it =
        cells.find(phenotype_id);

    ASSERT_(
        phenotype_it != cells.end(),
        "DetermineFemExportAction could not find phenotype ID "
        + std::to_string(phenotype_id)
    );

    const std::string& phenotype_name =
        phenotype_it->second;

    const std::string mech_base =
        phenotype_name + "/cell_matrix_mechanics";

    // -------------------------------------------------------------------------
    // Step 4: Exclude cells without enabled cell-matrix mechanics
    // -------------------------------------------------------------------------

    const bool mechanics_enabled =
        this->params()->have_parameter<bool>(mech_base + "/enabled") &&
        this->params()->get<bool>(mech_base + "/enabled");

    if (!mechanics_enabled) {
        return FemExportAction::kSkip;
    }

    // -------------------------------------------------------------------------
    // Step 5: Read the current attachment state
    // -------------------------------------------------------------------------

    const std::size_t attachment_count =
        cell->GetAttachmentNodeIds().size();

    const bool has_valid_mechanics_attachments =
        cell->HasValidMechanicsAttachments();

    // -------------------------------------------------------------------------
    // Step 6: Keep single-attachment cells entirely within the ABM
    // 
    // Assumption: single-attachment cells are not expected to be able to 
    // contract.
    // -------------------------------------------------------------------------

    if (attachment_count == 1) {
        return FemExportAction::kSkip;
    }

    // -------------------------------------------------------------------------
    // Step 7: Contract using two or more valid retained attachments
    //
    // The movement-flag condition is retained temporarily change remains 
    // compatible with the existing migration workflow, attached cells did not 
    // move.
    // -------------------------------------------------------------------------

    if (attachment_count >= 2 &&
        has_valid_mechanics_attachments &&
        !cell->GetMovedDueToMechanics()) {
        return FemExportAction::kContract;
    }

    // -------------------------------------------------------------------------
    // Step 8: Cells without usable contract attachments require FEM attachment
    //
    // At this stage this preserves the existing handling of:
    // - newly introduced cells with no attachments;
    // - cells moved by the legacy mechanics-migration behaviour (floating cells);
    // - cells whose attachment data are not currently valid.
    // -------------------------------------------------------------------------

    return FemExportAction::kAttach;
    }
    
  
    // ---------------------------------------------------------------------------
    // ABM-to-FEM export
    // ---------------------------------------------------------------------------

    inline
    void ExportCellPositions(
        bdm::Simulation& sim,
        const std::map<int, std::string>& cells,
        int time)
    {
    /*
    * Export ABM cell positions and FEM state instructions.
    *
    * Function goal
    * -------------
    * Write the ABM-to-FEM csv file containing each cell position, the requested
    * FEM handling state, reusable attachment node IDs, and the contractile force.
    *
    * File format
    * -----------
    * time,abm_cell_id,x,y,z,cell_state,attachment_node_ids,contractile_force
    */

    // ---------------------------------------------------------------------------
    // Step 1: Generate filename
    // ---------------------------------------------------------------------------

    std::ostringstream filename;

    filename << this->params()->get<std::string>("output_directory")
            << "/cell_positions/cells_t"
            << std::setw(4) << std::setfill('0') << time-1
            << ".csv";

    std::ofstream fpos(filename.str());

    if (!fpos.is_open()) {
        std::cerr << "Could not open file " << filename.str() << "\n";
        return;
    }

    // ---------------------------------------------------------------------------
    // Step 2: Reset the exported-cell record
    //
    // This set must describe only the cells included in the current FEM input.
    // We do not include cells that were skipped because they are necrotic,
    // have only 1 attachment point or from phenotypes where mechanics is 
    // disabled. 
    // ---------------------------------------------------------------------------

    exported_cell_ids_.clear();

    // ---------------------------------------------------------------------------
    // Step 3: Write header
    // ---------------------------------------------------------------------------

    fpos << "time,abm_cell_id,x,y,z,cell_state,attachment_node_ids,contractile_force\n";

    // ---------------------------------------------------------------------------
    // Step 4: Collect BiologicalCell agents safely
    // ---------------------------------------------------------------------------

    // ForEachAgent may be parallel depending on the BioDynaMo execution backend.
    // Therefore, only collect pointers here, and protect the shared vector with a
    // mutex. All CSV line construction, file writing, and cell state updates are
    // performed serially afterwards.

    std::vector<bdm::BiologicalCell*> biological_cells;
    std::mutex biological_cells_mutex;

    auto* rm = sim.GetResourceManager();

    rm->ForEachAgent([&](bdm::Agent* agent) {

        auto* cell = dynamic_cast<bdm::BiologicalCell*>(agent);

        if (!cell) {
        return;
        }

        std::lock_guard<std::mutex> lock(biological_cells_mutex);
        biological_cells.push_back(cell);
    });

    // Optional but useful: make output order deterministic.
    std::sort(
        biological_cells.begin(),
        biological_cells.end(),
        [](bdm::BiologicalCell* a, bdm::BiologicalCell* b) {
        std::ostringstream uid_a;
        std::ostringstream uid_b;

        uid_a << a->GetUid();
        uid_b << b->GetUid();

        return uid_a.str() < uid_b.str();
        }
    );

    // ---------------------------------------------------------------------------
    // Step 5: Build CSV lines serially
    // ---------------------------------------------------------------------------

    std::vector<std::string> csv_lines;
    csv_lines.reserve(biological_cells.size());

    for (auto* cell : biological_cells) {

        const auto& pos = cell->GetPosition();

        std::ostringstream uid_stream;
        uid_stream << cell->GetUid();

        const std::string abm_cell_id = uid_stream.str();

        // -------------------------------------------------------------------------
        // Step 5a: Determine whether and how this cell should be exported
        // -------------------------------------------------------------------------

        const FemExportAction export_action =
        this->DetermineFemExportAction(cell, cells);

        // Cells retained entirely within the ABM are not written to the FEM input.
        if (export_action == FemExportAction::kSkip) {
        continue;
        }

        // -------------------------------------------------------------------------
        // Step 5b: Prepare the FEM instruction
        // -------------------------------------------------------------------------

        std::string cell_state;
        std::string attachment_node_ids_text;
        double contractile_force = 0.0;

        if (export_action == FemExportAction::kContract) {

        cell_state = "contract";

        attachment_node_ids_text =
            this->FormatAttachmentNodeIdsJsonLike(
            cell->GetAttachmentNodeIds()
            );

        contractile_force =
            cell->GetContractileForce();

        // A negative contractile force is not expected by the model.
        // Report the issue and prevent it from reaching the FEM input.
        if (contractile_force < 0.0) {
            std::cerr
            << "[CELL-MATRIX WARNING] Negative contractile force for cell "
            << abm_cell_id
            << ": "
            << contractile_force
            << ". Force has been set to 0.0 before FEM export.\n";

            contractile_force = 0.0;
        }

        } else {

        // Sanity check: only the "attach" action is supported here.
        ASSERT_(
            export_action == FemExportAction::kAttach,
            "ExportCellPositions encountered an unsupported FEM export action"
        );

        cell_state = "attach";
        attachment_node_ids_text = "[]";
        contractile_force = 0.0;
        }

        // -------------------------------------------------------------------------
        // Step 5c: Record the FEM instruction on the cell
        // -------------------------------------------------------------------------

        cell->SetCellState(cell_state);

        // The movement flag has now been communicated to the FEM.
        // Skipped cells do not reach this line, so their flag remains unchanged.
        cell->ClearMovedDueToMechanics();

        // -------------------------------------------------------------------------
        // Step 5d: Record the exported ABM cell ID
        // -------------------------------------------------------------------------

        const bool inserted =
        exported_cell_ids_.insert(abm_cell_id).second;

        ASSERT_(
        inserted,
        "ExportCellPositions encountered duplicate ABM cell UID "
        + abm_cell_id
        );

        std::ostringstream line;

        line << time << ","
            << abm_cell_id << ","
            << pos[0] << ","
            << pos[1] << ","
            << pos[2] << ","
            << cell_state << ","
            << "\"" << attachment_node_ids_text << "\"" << ","
            << contractile_force << "\n";

        csv_lines.push_back(line.str());
    }

    // ---------------------------------------------------------------------------
    // Step 5: Write CSV lines serially and close the file
    // ---------------------------------------------------------------------------

    for (const auto& line : csv_lines) {
        fpos << line;
    }

    fpos.close();

    }

  // ---------------------------------------------------------------------------
  // FEM execution
  // ---------------------------------------------------------------------------

  inline
  int RunFemSolver(
      bdm::Simulation& sim,
      const std::map<int, std::string>& cells,
      int time,
      bool* ran_any)
    {
    /*
    * Run the FEM solver.
    *
    * Function goal
    * -------------
        * Run the FEM solver for the mechanics-enabled phenotype using only the cells
        * written during the current ABM-to-FEM export.
        *
        * The exported_cell_ids_ set is the authoritative source for cell_count.
        * Cells retained entirely within the ABM, including single-attachment cells,
        * are therefore excluded from the FEM calculation.
    *
    * Command format
    * --------------
    * python3 -u FEM_solver_interface.py \
            --step_num <time - 1> \
            --cell_count <cell_count> \
            --min_cell_radius <min_cell_radius> \
            --max_cell_radius <max_cell_radius> \
            --strut_radius <strut_radius> \
            --delta_F <delta_F> \
            --perturbance_dist <perturbance_dist> \
            --random_state <random_state> \
            --num_attachments <num_attachments> \
            --lattice_mesh_path "<lattice_mesh_path>" \
            --verbose \ <= IF TRUE 
            --private_key_path "<private_key_path>" \
            --user_name "<user_name>" \
            --host_name "<host_name>"
    */
    
    ASSERT_(
    ran_any != nullptr,
    "RunFemSolver received a null ran_any pointer"
    );
    
    // Initialise flag.
    *ran_any = false;

    // Temporary limitation: FEM solver supports only one phenotype at a time.
    int mech_enabled_phenotypes = 0;
    std::string enabled_list;

    for (const auto& kv : cells) {
        const int CP_ID = kv.first;
        const std::string& CP_name = kv.second;

        if (CP_ID < 1) {
        continue;
        }

        const std::string mech_base =
        CP_name + "/cell_matrix_mechanics";

        if (this->params()->have_parameter<bool>(mech_base + "/enabled") &&
            this->params()->get<bool>(mech_base + "/enabled")) {

        ++mech_enabled_phenotypes;

        if (!enabled_list.empty()) {
            enabled_list += ", ";
        }

        enabled_list +=
            CP_name + " (ID " + std::to_string(CP_ID) + ")";
        }
    }

    if (mech_enabled_phenotypes > 1) {
        ABORT_(
        "Temporary limitation: FEM solver currently supports only ONE cell "
        "phenotype with cell-matrix mechanics enabled. Enabled phenotypes: "
        + enabled_list
        + ". Please enable mechanics for only one phenotype. Future versions "
            "will support multiple phenotypes."
        );
    }

    // Check whether any phenotype has mechanics enabled.
    bool any_mech_enabled = false;

    for (auto ci = cells.begin(); ci != cells.end(); ++ci) {
        const int CP_ID = ci->first;
        const std::string& CP_name = ci->second;

        if (CP_ID < 1) {
        continue;
        }

        const std::string mech_base =
        CP_name + "/cell_matrix_mechanics";

        if (this->params()->have_parameter<bool>(mech_base + "/enabled") &&
            this->params()->get<bool>(mech_base + "/enabled")) {

        any_mech_enabled = true;
        break;
        }
    }

    // Read global HPC settings once.
    std::string private_key_path;
    std::string user_name;
    std::string host_name;

    if (any_mech_enabled) {
        private_key_path =
        this->params()->get<std::string>(
            "cell_matrix_mechanics/HPC/private_key_path"
        );

        user_name =
        this->params()->get<std::string>(
            "cell_matrix_mechanics/HPC/user_name"
        );

        host_name =
        this->params()->get<std::string>(
            "cell_matrix_mechanics/HPC/host_name"
        );
    }

    for (auto ci = cells.begin(); ci != cells.end(); ++ci) {
        const int CP_ID = ci->first;
        const std::string& CP_name = ci->second;

        // Ignore the necrotic phenotype.
        if (CP_ID < 1) {
        continue;
        }

        const std::string mech_base =
        CP_name + "/cell_matrix_mechanics";

        // Skip phenotypes without cell-matrix mechanics.
        if (!this->params()->have_parameter<bool>(mech_base + "/enabled") ||
            !this->params()->get<bool>(mech_base + "/enabled")) {
        continue;
        }

        // -------------------------------------------------------------------------
        // Phenotype-specific mechanics parameters
        // -------------------------------------------------------------------------

        const double min_cell_radius =
            this->params()->get<double>(
                mech_base + "/min_cell_reach_radius"
            );

        if ((2.0 * min_cell_radius) <
            this->params()->get<double>(CP_name + "/diameter/min")) {
            ABORT_(
                "Model parameter '"
                + mech_base
                + "/min_cell_reach_radius' cannot be less than '"
                + CP_name
                + "/diameter/min'"
            );
            }

        const double max_cell_radius =
            this->params()->get<double>(
                mech_base + "/max_cell_reach_radius"
            );

        const double strut_radius =
            this->params()->get<double>(
                mech_base + "/strut_radius"
            );

        const double delta_F =
            this->params()->get<double>(
                mech_base + "/delta_F"
            );

        const double perturbance_dist =
            this->params()->get<double>(
                mech_base + "/perturbance_dist"
            );

        const int random_state =
            this->params()->get<int>(
                mech_base + "/random_state"
            );

        const int num_attachments =
            this->params()->get<int>(
                mech_base + "/num_attachments"
            );

        const std::string lattice_mesh_path =
            this->params()->get<std::string>(
                mech_base + "/lattice_mesh_path"
            );

        const bool verbose =
            this->params()->get<bool>(
                mech_base + "/verbose"
            );

        // -------------------------------------------------------------------------
        // Determine the number of cells included in the current FEM input
        //
        // ExportCellPositions() has already excluded cells that must remain entirely
        // within the ABM. The exported UID set therefore defines the FEM cell count.
        // -------------------------------------------------------------------------

        const int cell_count =
            static_cast<int>(exported_cell_ids_.size());

        // Do not launch the FEM solver when no eligible cells were exported.
        if (cell_count == 0) {
            this->PrepareNoInteractionStep(time); // Take undeformed scaffold
            continue;
        }
        
        // FEM solver will run for this phenotype.
        *ran_any = true;

        // -------------------------------------------------------------------------
        // Build Python command
        // -------------------------------------------------------------------------

        std::string cmd =
        "python3 -u FEM_solver_interface.py ";

        cmd +=
        "--step_num "
        + std::to_string(time - 1)
        + " ";

        cmd +=
        "--cell_count "
        + std::to_string(cell_count)
        + " ";

        cmd +=
        "--min_cell_radius "
        + std::to_string(min_cell_radius)
        + " ";

        cmd +=
        "--max_cell_radius "
        + std::to_string(max_cell_radius)
        + " ";

        cmd +=
        "--strut_radius "
        + std::to_string(strut_radius)
        + " ";

        cmd +=
        "--delta_F "
        + std::to_string(delta_F)
        + " ";

        cmd +=
        "--perturbance_dist "
        + std::to_string(perturbance_dist)
        + " ";

        cmd +=
        "--random_state "
        + std::to_string(random_state)
        + " ";

        cmd +=
        "--num_attachments "
        + std::to_string(num_attachments)
        + " ";

        cmd +=
        "--lattice_mesh_path \""
        + lattice_mesh_path
        + "\" ";

        if (verbose) {
        cmd += "--verbose ";
        }

        // Global HPC arguments.
        cmd +=
        "--private_key_path \""
        + private_key_path
        + "\" ";

        cmd +=
        "--user_name \""
        + user_name
        + "\" ";

        cmd +=
        "--host_name \""
        + host_name
        + "\" ";

        std::cout
        << "Running command: "
        << cmd
        << std::endl;

        FILE* pipe =
        popen(cmd.c_str(), "r");

        if (!pipe) {
        std::cerr
            << "Failed to start Python script for phenotype "
            << CP_name
            << "\n";

        return 1;
        }

        char buffer[256];

        while (fgets(buffer, sizeof(buffer), pipe)) {
        std::cout
            << buffer
            << std::flush;
        }

        const int return_code =
        pclose(pipe);

        if (return_code != 0) {
        std::cerr
            << "Python script failed for phenotype "
            << CP_name
            << " with code "
            << return_code
            << "\n";

        return return_code;
        }
    }

    return 0;
    }

    // ---------------------------------------------------------------------------
    // Write scaffold using lattice format version 2.0
    // ---------------------------------------------------------------------------

    inline
    void WriteScaffoldVersion2(
        const ObstacleScaffold& scaffold,
        const std::string& output_path) const
    {
    /*
    * Function goal
    * -------------
    * Write an in-memory scaffold using persistent one-based node and element
    * identifiers in lattice format version 2.0.
    *
    * This method is also used to translate a legacy version 1.0 scaffold after
    * it has been read by ObstacleScaffold.
    */

    // -------------------------------------------------------------------------
    // Step 1: Validate the scaffold data
    // -------------------------------------------------------------------------

    ASSERT_(
        !scaffold.nodes_by_id.empty(),
        "WriteScaffoldVersion2 received a scaffold with no nodes"
    );

    ASSERT_(
        !scaffold.segment.empty(),
        "WriteScaffoldVersion2 received a scaffold with no segments"
    );

    // -------------------------------------------------------------------------
    // Step 2: Sort node IDs for deterministic output
    // -------------------------------------------------------------------------

    std::vector<int> node_ids;
    node_ids.reserve(scaffold.nodes_by_id.size());

    for (const auto& item : scaffold.nodes_by_id) {
        node_ids.push_back(item.first);
    }

    std::sort(node_ids.begin(), node_ids.end());

    // -------------------------------------------------------------------------
    // Step 3: Sort segments by persistent element ID
    // -------------------------------------------------------------------------

    std::vector<const ObstacleScaffold::Segment*> sorted_segments;
    sorted_segments.reserve(scaffold.segment.size());

    for (const auto& segment : scaffold.segment) {
        sorted_segments.push_back(&segment);
    }

    std::sort(
        sorted_segments.begin(),
        sorted_segments.end(),
        [](const ObstacleScaffold::Segment* a,
        const ObstacleScaffold::Segment* b) {
        return a->element_id < b->element_id;
        }
    );

    // -------------------------------------------------------------------------
    // Step 4: Open the output file
    // -------------------------------------------------------------------------

    std::ofstream fout(output_path);

    ASSERT_(
        fout.good(),
        "WriteScaffoldVersion2 could not create file: "
        + output_path
    );

    fout << std::fixed << std::setprecision(6);

    // -------------------------------------------------------------------------
    // Step 5: Write the version and node section
    // -------------------------------------------------------------------------

    fout << "2.0\n";
    fout << node_ids.size() << "\n";

    for (const int node_id : node_ids) {
        const auto& node =
        scaffold.GetNode(node_id);

        fout << node.id << " "
            << node.position[0] << " "
            << node.position[1] << " "
            << node.position[2] << " "
            << node.radius << "\n";
    }

    // -------------------------------------------------------------------------
    // Step 6: Write the element section
    // -------------------------------------------------------------------------

    fout << sorted_segments.size() << "\n";

    for (const auto* segment : sorted_segments) {
        fout << segment->element_id << " "
            << segment->node_id_1 << " "
            << segment->node_id_2 << "\n";
    }

    ASSERT_(
        fout.good(),
        "WriteScaffoldVersion2 failed while writing file: "
        + output_path
    );

    fout.close();
    }
  
    // ---------------------------------------------------------------------------
    // Prepare scaffold for a timestep with no FEM-eligible cells
    // ---------------------------------------------------------------------------

    inline
    void PrepareNoInteractionStep(const int time) const
    {
    /*
    * Function goal
    * -------------
    * Prepare the unloaded scaffold for a timestep in which no cells are sent
    * to the FEM.
    *
    * The saved initial scaffold is read through ObstacleScaffold, allowing
    * legacy version 1.0 input to be translated into the persistent-ID model.
    * The timestep lattice is then written using version 2.0.
    * 
    * Version 2.0 format:
    *
    *   2.0
    *   <number_of_nodes>
    *   <node_id> <x> <y> <z> <radius>
    *   ...
    *   <number_of_elements>
    *   <element_id> <node_id_1> <node_id_2>
    *
    * No cell-mechanics file is produced because no cells were interacted with
    * the scaffold at this timestep.
    */

    // -------------------------------------------------------------------------
    // Step 1: Define the initial and timestep scaffold paths
    // -------------------------------------------------------------------------

    const std::string initial_lattice_path =
        this->params()->get<std::string>("output_directory")
        + "/in/simulation_obstacle.1.scaffold";

    const std::string step_directory =
        this->params()->get<std::string>("output_directory")
        + "/FEM/step_"
        + std::to_string(time - 1);

    const std::string target_lattice_path =
        step_directory + "/lattice.1d";

    const std::string mechanics_file_path =
        step_directory
        + "/cell_mechanics_step_"
        + std::to_string(time - 1)
        + ".dat";

    // -------------------------------------------------------------------------
    // Step 2: Confirm that the saved initial scaffold exists
    // -------------------------------------------------------------------------

    std::ifstream initial_file(initial_lattice_path);

    ASSERT_(
        initial_file.good(),
        "PrepareNoInteractionStep could not find the initial scaffold: "
        + initial_lattice_path
    );

    initial_file.close();

    // -------------------------------------------------------------------------
    // Step 3: Create the timestep output directory
    // -------------------------------------------------------------------------

    const std::string mkdir_command =
        "mkdir -p \"" + step_directory + "\"";

    ASSERT_(
        0 == std::system(mkdir_command.c_str()),
        "PrepareNoInteractionStep could not create directory: "
        + step_directory
    );

    // -------------------------------------------------------------------------
    // Step 4: Read and validate the initial scaffold
    //
    // ObstacleScaffold translates legacy version 1.0 indexing into persistent
    // one-based node and element identifiers.
    // -------------------------------------------------------------------------

    ObstacleScaffold initial_scaffold;

    initial_scaffold.init(
        "scaffold",
        initial_lattice_path
    );

    // -------------------------------------------------------------------------
    // Step 5: Write the unloaded scaffold using version 2.0
    // -------------------------------------------------------------------------

    this->WriteScaffoldVersion2(
        initial_scaffold,
        target_lattice_path
    );

    // -------------------------------------------------------------------------
    // Step 6: Remove any stale mechanics file
    //
    // A previous run may have left a file in the same timestep directory.
    // No mechanics file should exist when no cells were processed.
    // -------------------------------------------------------------------------

    std::remove(mechanics_file_path.c_str());
    }
  
    // ---------------------------------------------------------------------------
    // FEM-to-ABM import
    // ---------------------------------------------------------------------------

    inline
    void ImportFemCells(
    bdm::Simulation& sim,
    const std::map<int, std::string>& cells,
    int time)
    {
    /*
    * Function goal
    * -------------
    * Read FEM-calculated cell mechanics data and assign it to the corresponding
    * ABM cells using their BioDynaMo UID strings.
    *
    * Only cells included in the current ABM-to-FEM export are expected in the FEM
    * result file. Cells retained entirely within the ABM are left unchanged.
    *
    * Input file format
    * -----------------
    * The first line contains the number of cells processed by the FEM:
    *
    *   <number_of_cells>
    *
    * Each subsequent row contains:
    *
    *   <abm_cell_id>
    *   <cell_x> <cell_y> <cell_z>
    *   <k_ce>
    *   <number_of_attachments>
    *   <attachment_node_ids>
    *   <attachment_coordinates>
    *   <attachment_k_ecm_values>
    *
    * For a cell with N attachments:
    *
    *   abm_cell_id
    *   cell_x cell_y cell_z
    *   k_ce
    *   N
    *   node_id_1 ... node_id_N
    *   attach_x_1 attach_y_1 attach_z_1
    *   ...
    *   attach_x_N attach_y_N attach_z_N
    *   k_ecm_1 ... k_ecm_N
    *
    * Attachment node IDs use persistent one-based scaffold/FEM IDs.
    * Example for a cell with two attachments:
    *
    *   0-0 249.2765 504.1778 595.3332 41.6827 2
    *   55495 78082
    *   270.2821 518.9576 592.0374
    *   230.3962 489.4993 593.8962
    *   42.3303 23.3542
    *
    * Notes
    * -----
    * - `abm_cell_id` is the BioDynaMo UID string, for example "0-0".
    * - Attachment node IDs use the persistent one-based FEM/scaffold IDs.
    * - Attachment coordinates represent the current deformed scaffold geometry.
    * - One k_ecm value is provided for each attachment.
    * - The number of node IDs, coordinates and k_ecm values must all equal N.
    */

    // ---------------------------------------------------------------------------
    // Step 1: Build FEM result filename
    // ---------------------------------------------------------------------------

    char buf[2048];

    const int time_inc = time - 1;

    std::snprintf(buf, sizeof(buf),
                    "./results/FEM/step_%d/cell_mechanics_step_%d.dat",
                    time_inc, time_inc);

    const std::string fname(buf);

    // ---------------------------------------------------------------------------
    // Step 2: Open FEM result file
    // ---------------------------------------------------------------------------

    std::ifstream fin(fname);

    ASSERT_(fin, "FEM import: could not open file " + fname);

    // ---------------------------------------------------------------------------
    // Step 3: Read number of cells
    // ---------------------------------------------------------------------------

    int n_cells = -1;

    fin >> n_cells;

    ASSERT_(fin && n_cells >= 0,
            "FEM import: invalid n_cells in " + fname);

    // ---------------------------------------------------------------------------
    // Step 4: Build ABM cell lookup using BioDynaMo UID string
    // ---------------------------------------------------------------------------

    // Collect BiologicalCell pointers safely first.
    std::vector<bdm::BiologicalCell*> biological_cells;
    std::mutex biological_cells_mutex;

    auto* rm = sim.GetResourceManager();

    rm->ForEachAgent([&](bdm::Agent* agent) {
        auto* cell = dynamic_cast<bdm::BiologicalCell*>(agent);

        if (!cell) {
        return;
        }

        std::lock_guard<std::mutex> lock(biological_cells_mutex);
        biological_cells.push_back(cell);
    });

    // Build the lookup serially after the agent collection step.
    std::unordered_map<std::string, bdm::BiologicalCell*> cell_lookup;

    for (auto* cell : biological_cells) {
        std::ostringstream uid_stream;
        uid_stream << cell->GetUid();

        const std::string abm_cell_id = uid_stream.str();

        ASSERT_(cell_lookup.find(abm_cell_id) == cell_lookup.end(),
                "FEM import: duplicate ABM cell UID found: " + abm_cell_id);

        cell_lookup[abm_cell_id] = cell;
    }

    // ---------------------------------------------------------------------------
    // Step 5: Validate the expected FEM result count
    //
    // The FEM result file should contain one row for every cell exported during
    // the current coupling step, not one row for every ABM cell.
    // ---------------------------------------------------------------------------

    ASSERT_(
        n_cells == static_cast<int>(exported_cell_ids_.size()),
        "FEM import: file n_cells does not match the number of cells exported "
        "to the FEM"
    );

    // Confirm that every exported UID still identifies a BiologicalCell.
    for (const auto& exported_id : exported_cell_ids_) {
        ASSERT_(
            cell_lookup.find(exported_id) != cell_lookup.end(),
            "FEM import: exported ABM cell UID no longer exists: "
            + exported_id
        );
    }

    // ---------------------------------------------------------------------------
    // Step 5: Track only the cells expected from the FEM
    //
    // Single-attachment cells are not exported to the FEM, so they are not expected 
    // to be present in the FEM result file. 
    // ---------------------------------------------------------------------------

    std::unordered_map<std::string, bool> imported_flags;

    for (const auto& exported_id : exported_cell_ids_) {
        imported_flags.emplace(exported_id, false);
    }

    const int debug_cells_to_print = 3;

    for (int i = 0; i < n_cells; ++i) {

        std::string abm_cell_id;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double k_ce = 0.0;
        int n_attach = 0;

        fin >> abm_cell_id >> x >> y >> z >> k_ce >> n_attach;

        ASSERT_(fin && n_attach >= 0,
                "FEM import: failed parsing header for row " + std::to_string(i)
                + " in " + fname);

        auto cell_it = cell_lookup.find(abm_cell_id);

        ASSERT_(cell_it != cell_lookup.end(),
                "FEM import: could not find ABM cell with abm_cell_id = "
                + abm_cell_id);

        // The FEM must not return cells that were omitted from the current export.
        auto imported_flag_it =
        imported_flags.find(abm_cell_id);

        ASSERT_(
            imported_flag_it != imported_flags.end(),
            "FEM import: FEM returned unexpected or non-exported abm_cell_id = "
            + abm_cell_id
        );

        ASSERT_(
            !imported_flag_it->second,
            "FEM import: duplicate FEM row for abm_cell_id = "
            + abm_cell_id
        );

        bdm::BiologicalCell* cell = cell_it->second;

        // -------------------------------------------------------------------------
        // Step 6a: Read attachment node IDs
        // -------------------------------------------------------------------------

        std::vector<int> attachment_node_ids(n_attach);

        for (int a = 0; a < n_attach; ++a) {
        fin >> attachment_node_ids[a];

        ASSERT_(fin,
                "FEM import: not enough attachment node IDs for abm_cell_id = "
                + abm_cell_id + " in " + fname);
        }

        // -------------------------------------------------------------------------
        // Step 6b: Read attachment coordinates
        // -------------------------------------------------------------------------

        std::vector<bdm::Double3> attachment_points(n_attach);

        for (int a = 0; a < n_attach; ++a) {
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;

        fin >> ax >> ay >> az;

        ASSERT_(fin,
                "FEM import: not enough attachment xyz values for abm_cell_id = "
                + abm_cell_id + " in " + fname);

        attachment_points[a] = bdm::Double3{ax, ay, az};
        }

        // -------------------------------------------------------------------------
        // Step 6c: Read attachment stiffness values
        // -------------------------------------------------------------------------

        std::vector<double> k_values(n_attach);

        for (int a = 0; a < n_attach; ++a) {
        fin >> k_values[a];

        ASSERT_(fin,
                "FEM import: not enough stiffness values for abm_cell_id = "
                + abm_cell_id + " in " + fname);
        }

        // -------------------------------------------------------------------------
        // Step 6d: Update the matched ABM cell
        // -------------------------------------------------------------------------

        cell->SetPosition(bdm::Double3{x, y, z});

        cell->ClearAttachmentNodeIds();
        cell->ClearAttachmentPoints();
        cell->ClearAttachmentStiffness();

        // Set points and stiffness before node IDs, or IDs before both, both are
        // valid because the vectors are empty after clearing. This order keeps the
        // geometry and stiffness assignment similar to the original implementation.
        cell->SetAttachmentPoints(attachment_points);
        cell->SetAttachmentStiffness(k_values);
        cell->SetAttachmentNodeIds(attachment_node_ids);

        // -------------------------------------------------------------------------
        // Temporary Stage 6 test: force one selected cell to retain one attachment
        // -------------------------------------------------------------------------

        const std::string single_attachment_test_uid = "0-0";

        if (abm_cell_id == single_attachment_test_uid &&
            cell->GetAttachmentNodeIds().size() > 1) {

        const std::vector<int> test_node_ids = {
            cell->GetAttachmentNodeIds().front()
        };

        const std::vector<bdm::Double3> test_points = {
            cell->GetAttachmentPoints().front()
        };

        const std::vector<double> test_stiffness = {
            cell->GetAttachmentStiffness().front()
        };

        cell->ClearAttachmentNodeIds();
        cell->ClearAttachmentPoints();
        cell->ClearAttachmentStiffness();

        cell->SetAttachmentPoints(test_points);
        cell->SetAttachmentStiffness(test_stiffness);
        cell->SetAttachmentNodeIds(test_node_ids);

        std::cout
            << "[STAGE 6 TEST] Forced cell "
            << abm_cell_id
            << " to retain one attachment at node "
            << test_node_ids.front()
            << "\n";
        }

        cell->SetKce(k_ce);

        // -------------------------------------------------------------------------
        // Step 6e: Update adaptive cell-specific contractile force
        // -------------------------------------------------------------------------

        double adaptive_contractile_force = 0.0;

        const int phenotype_id = cell->GetPhenotype();

        auto phenotype_it = cells.find(phenotype_id);

        if (phenotype_it != cells.end() && phenotype_id >= 1) {

        const std::string& CP_name = phenotype_it->second;
        const std::string mech_base = CP_name + "/cell_matrix_mechanics";

        const bool mechanics_enabled =
            this->params()->have_parameter<bool>(mech_base + "/enabled") &&
            this->params()->get<bool>(mech_base + "/enabled");

        if (mechanics_enabled) {

            const double min_contractile_force =
                this->params()->get<double>(mech_base + "/min_contractile_force");

            const double max_contractile_force =
                this->params()->get<double>(mech_base + "/max_contractile_force");

            const double contractility_kce_sensitivity =
                this->params()->get<double>(mech_base + "/contractility_kce_sensitivity");

            const double contractility_kce_midpoint =
                this->params()->get<double>(mech_base + "/contractility_kce_midpoint");

            adaptive_contractile_force =
                calculate_contractile_force_from_kce(
                    k_ce,
                    min_contractile_force,
                    max_contractile_force,
                    contractility_kce_sensitivity,
                    contractility_kce_midpoint
                );
        }
        }

        cell->SetContractileForce(adaptive_contractile_force);


        // Do not update cell_state here.
        // cell_state_ represents the last ABM-to-FEM state exported by ABM.
        // This prevents cells that were just attached by FEM from migrating before
        // they have completed a dedicated "contract" step.

        imported_flag_it->second = true;

        if (i < debug_cells_to_print) {
        std::cout << "[FEM IMPORT] abm_cell_id=" << abm_cell_id
                    << " pos=(" << x << ", " << y << ", " << z << ")"
                    << " k_ce=" << k_ce
                    << " contractile_force=" << adaptive_contractile_force
                    << " n_attach=" << n_attach
                    << "\n";
        }
    }

    // ---------------------------------------------------------------------------
    // Step 7: Confirm every ABM cell was imported exactly once
    // ---------------------------------------------------------------------------

    for (const auto& item : imported_flags) {
        ASSERT_(
            item.second,
            "FEM import: no FEM row was returned for exported abm_cell_id = "
            + item.first
        );
    }

    fin.close();

    }

    private:
    // Model parameters are owned by ABM4bio, not by this class.
    Parameters* params_;

    /*
    * This will be populated during Stage 6.
    *
    * It will allow ImportFemCells() to validate exactly which ABM cells were
    * exported to and returned by the FEM.
    */
    std::unordered_set<std::string> exported_cell_ids_;
};


// =============================================================================
// Method definitions
// =============================================================================

#endif  // CELL_MATRIX_INTERACTION_H_