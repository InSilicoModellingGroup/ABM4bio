#ifndef CELL_MATRIX_INTERACTIONS_H_
#define CELL_MATRIX_INTERACTIONS_H_

// =============================================================================
// Standard-library dependencies
// =============================================================================

/*

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

*/

// =============================================================================
// Project dependencies
// =============================================================================

#include "./global.h"
#include "./biological_cell.h"

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
  // Step 2: Write header
  // ---------------------------------------------------------------------------

  fpos << "time,abm_cell_id,x,y,z,cell_state,attachment_node_ids,contractile_force\n";

  // ---------------------------------------------------------------------------
  // Step 3: Collect BiologicalCell agents safely
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
  // Step 4: Build CSV lines serially
  // ---------------------------------------------------------------------------

  std::vector<std::string> csv_lines;
  csv_lines.reserve(biological_cells.size());

  for (auto* cell : biological_cells) {

    const auto& pos = cell->GetPosition();

    std::ostringstream uid_stream;
    uid_stream << cell->GetUid();

    const std::string abm_cell_id = uid_stream.str();

    const bool has_valid_mechanics_attachments =
      cell->HasValidMechanicsAttachments();

    const bool has_enough_contract_attachments =
      has_valid_mechanics_attachments &&
      cell->GetAttachmentNodeIds().size() >= 2;

    std::string cell_state = "attach";
    std::string attachment_node_ids_text = "[]";
    double contractile_force = 0.0;

    if (has_enough_contract_attachments && !cell->GetMovedDueToMechanics()) {
      cell_state = "contract";
      attachment_node_ids_text =
        this->FormatAttachmentNodeIdsJsonLike(cell->GetAttachmentNodeIds());

      contractile_force = cell->GetContractileForce();
    } else {
      cell_state = "attach";
      attachment_node_ids_text = "[]";
      contractile_force = 0.0;
    }

    if (contractile_force < 0.0) {
      contractile_force = 0.0;
    }

    // Store the state that ABM is sending to FEM.
    // CheckMigration() will later use this as the previous FEM-request state.
    cell->SetCellState(cell_state);

    // The movement flag has now been communicated to FEM.
    cell->ClearMovedDueToMechanics();

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
   * Call the pyton code (FEM_solver_interface.py) to transfer files, run and monitor FEM.
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
    // Count cells belonging to this phenotype
    // -------------------------------------------------------------------------

    int cell_count = 0;

    auto* rm = sim.GetResourceManager();

    rm->ForEachAgent([&](bdm::Agent* agent) {
      auto* cell =
        dynamic_cast<bdm::BiologicalCell*>(agent);

      if (!cell) {
        return;
      }

      if (cell->GetPhenotype() == CP_ID) {
        ++cell_count;
      }
    });

    if (cell_count == 0) {
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
    * Read information related to the cells mechanical state and assign it to the ABM cells.
    *
    * Input file format
    * -----------------
    * The first line contains the number of cells included in the FEM results:
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
    * More explicitly, for a cell with N attachments:
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

        ASSERT_(n_cells == static_cast<int>(cell_lookup.size()),
                "FEM import: file n_cells does not match number of BiologicalCell agents");

        // ---------------------------------------------------------------------------
        // Step 5: Read each FEM result row and update the matching ABM cell
        // ---------------------------------------------------------------------------

        std::unordered_map<std::string, bool> imported_flags;

        for (const auto& item : cell_lookup) {
            imported_flags[item.first] = false;
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

            ASSERT_(!imported_flags[abm_cell_id],
                    "FEM import: duplicate FEM row for abm_cell_id = " + abm_cell_id);

            bdm::BiologicalCell* cell = cell_it->second;

            // -------------------------------------------------------------------------
            // Step 5a: Read attachment node IDs
            // -------------------------------------------------------------------------

            std::vector<int> attachment_node_ids(n_attach);

            for (int a = 0; a < n_attach; ++a) {
            fin >> attachment_node_ids[a];

            ASSERT_(fin,
                    "FEM import: not enough attachment node IDs for abm_cell_id = "
                    + abm_cell_id + " in " + fname);
            }

            // -------------------------------------------------------------------------
            // Step 5b: Read attachment coordinates
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
            // Step 5c: Read attachment stiffness values
            // -------------------------------------------------------------------------

            std::vector<double> k_values(n_attach);

            for (int a = 0; a < n_attach; ++a) {
            fin >> k_values[a];

            ASSERT_(fin,
                    "FEM import: not enough stiffness values for abm_cell_id = "
                    + abm_cell_id + " in " + fname);
            }

            // -------------------------------------------------------------------------
            // Step 5d: Update the matched ABM cell
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

            cell->SetKce(k_ce);

            // -------------------------------------------------------------------------
            // Step 5e: Update adaptive cell-specific contractile force
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

            imported_flags[abm_cell_id] = true;

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
        // Step 6: Confirm every ABM cell was imported exactly once
        // ---------------------------------------------------------------------------

        for (const auto& item : imported_flags) {
            ASSERT_(item.second,
                    "FEM import: no FEM row was imported for abm_cell_id = "
                    + item.first);
        }

        // ---------------------------------------------------------------------------
        // Step 7: Check for trailing unread data issues
        // ---------------------------------------------------------------------------

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

// Method implementations will follow here.

#endif  // CELL_MATRIX_INTERACTION_H_