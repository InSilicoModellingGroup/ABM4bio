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

    // ---------------------------------------------------------------------------
    // Initial-attachment cell-path spatial cache
    // ---------------------------------------------------------------------------

    struct InitialAttachmentPathRecord {
        unsigned int cell_uid = 0;

        bdm::Double3 start;
        bdm::Double3 end;
        bdm::Double3 midpoint;

        double length = 0.0;
    };

    struct InitialAttachmentRelocationCandidate {
        bool found = false;

        int segment_id = -1;

        bdm::Double3 closest_point;
        bdm::Double3 proposed_position;

        double original_centreline_distance = 0.0;
        double target_centreline_distance = 0.0;
        double movement_distance = 0.0;
    };

    struct InitialAttachmentScaffoldOverlapResult {
        bool overlaps = false;

        int blocking_segment_id = -1;

        double centreline_distance = 0.0;
        double required_distance = 0.0;
        double clearance_margin = 0.0;
    };

    struct InitialAttachmentGenerationResult {
        bool success = false;

        std::size_t candidate_node_count = 0;
        std::size_t candidate_pair_count = 0;
        std::size_t scaffold_valid_pair_count = 0;
        std::size_t cell_valid_pair_count = 0;

        std::pair<int, int> selected_pair = {
            -1,
            -1
        };
    };

    struct InitialAttachmentPathBucketKey {
        int x;
        int y;
        int z;

        bool operator==(
            const InitialAttachmentPathBucketKey& other) const {

            return x == other.x &&
                y == other.y &&
                z == other.z;
        }
    };

    struct InitialAttachmentPathBucketKeyHash {
        std::size_t operator()(
            const InitialAttachmentPathBucketKey& key) const {

            std::size_t seed = 0;

            seed ^= std::hash<int>()(key.x)
                    + 0x9e3779b9
                    + (seed << 6)
                    + (seed >> 2);

            seed ^= std::hash<int>()(key.y)
                    + 0x9e3779b9
                    + (seed << 6)
                    + (seed >> 2);

            seed ^= std::hash<int>()(key.z)
                    + 0x9e3779b9
                    + (seed << 6)
                    + (seed >> 2);

            return seed;
        }
    };

    struct InitialAttachmentCellRecord {
        unsigned int cell_uid = 0;
        bdm::Double3 position;
        double radius = 0.0;
    };

    struct InitialAttachmentCellBucketKey {
        int x;
        int y;
        int z;

        bool operator==(
            const InitialAttachmentCellBucketKey& other) const {

            return
                x == other.x &&
                y == other.y &&
                z == other.z;
        }
    };

    struct InitialAttachmentCellBucketKeyHash {

        std::size_t operator()(
            const InitialAttachmentCellBucketKey& key) const {

            std::size_t seed = 0;

            seed ^=
                std::hash<int>()(key.x) +
                0x9e3779b9 +
                (seed << 6) +
                (seed >> 2);

            seed ^=
                std::hash<int>()(key.y) +
                0x9e3779b9 +
                (seed << 6) +
                (seed >> 2);

            seed ^=
                std::hash<int>()(key.z) +
                0x9e3779b9 +
                (seed << 6) +
                (seed >> 2);

            return seed;
        }
    };

    struct InitialAttachmentCellOverlapResult {
        bool overlaps = false;

        unsigned int blocking_cell_uid = 0;

        double centre_distance = 0.0;
        double required_distance = 0.0;
        double clearance_margin = 0.0;
    };

    struct InitialAttachmentRelocationCheckResult {
        bool accepted = false;

        InitialAttachmentRelocationCandidate relocation;

        InitialAttachmentScaffoldOverlapResult scaffold_overlap;

        InitialAttachmentCellOverlapResult cell_overlap;
    };

    struct InitialAttachmentRelocationSelectionResult {
        bool accepted = false;

        std::size_t selected_position_index = 0;

        bdm::Double3 selected_position;

        InitialAttachmentScaffoldOverlapResult scaffold_overlap;

        InitialAttachmentCellOverlapResult cell_overlap;
    };

    struct InitialAttachmentRelocationFallbackResult {
        bool accepted = false;

        std::size_t attempted_segment_count = 0;

        std::size_t selected_segment_index = 0;
        int selected_segment_id = -1;

        std::size_t selected_position_index = 0;

        bdm::Double3 selected_position;
    };

    struct InitialAttachmentWorkflowResult {

        bool success = false;
        bool attached_at_original_position = false;
        bool original_position_blocked_by_cell = false;
        bool relocation_attempted = false;
        bool marked_for_death = false;
        bool removal_requested = false;

        std::size_t attempted_segment_count = 0;
        std::size_t attempted_position_count = 0;
        std::size_t valid_relocation_position_count = 0;

        int selected_segment_id = -1;

        std::size_t selected_position_index = 0;

        bdm::Double3 selected_relocation_position;

        InitialAttachmentGenerationResult original_generation_result;
        InitialAttachmentGenerationResult final_generation_result;
    };
    
    struct InitialAttachmentPopulationResult {

        std::size_t processed_cell_count = 0;
        std::size_t attached_cell_count = 0;
        std::size_t attached_at_original_count = 0;
        std::size_t attached_after_relocation_count = 0;
        std::size_t original_position_blocked_count = 0;
        std::size_t removed_cell_count = 0;
        std::size_t total_segment_attempt_count = 0;
        std::size_t total_position_attempt_count = 0;
        std::size_t seeded_established_cell_count = 0;

        unsigned int first_processed_uid = 0;

        unsigned int last_processed_uid = 0;
    };

    struct InitialAttachmentCellPathObstructionResult {

        bool obstructed = false;

        unsigned int blocking_cell_uid = 0;

        double path_distance = 0.0;

        double required_distance = 0.0;
    };
    
    struct InitialAdditionalAttachmentSelectionResult {

        bool success = false;

        std::vector<double> candidate_distances;
        std::vector<double> candidate_weights;
        std::vector<double> selection_probabilities;

        int selected_node_id = -1;
    };
    
    struct InitialAdditionalAttachmentCompletionResult {

        std::size_t initial_attachment_count = 0;

        std::size_t target_attachment_count = 0;

        std::size_t formed_attachment_count = 0;

        std::size_t final_attachment_count = 0;

        bool reached_target = false;

        bool stopped_no_candidates = false;
    };
    
    inline
    void ResetInitialAttachmentCellCache(
        const double bucket_size) {
        ASSERT_(
            std::isfinite(bucket_size) && bucket_size > 0.0,
            "Initial attachment cell-cache bucket size "
            "must be finite and positive"
        );

        initial_attachment_cell_records_.clear();
        initial_attachment_cell_spatial_index_.clear();
        initial_attachment_cell_bucket_size_ = bucket_size;
        initial_attachment_cell_max_radius_ = 0.0;
    }

    inline
    void AddInitialAttachmentCellToCache(
        const unsigned int cell_uid,
        const bdm::Double3& position,
        const double radius) {

        ASSERT_(
            std::isfinite(initial_attachment_cell_bucket_size_) && initial_attachment_cell_bucket_size_ > 0.0,
            "Initial attachment cell cache must be reset before adding cells"
        );

        ASSERT_(
            std::isfinite(radius) && radius > 0.0,
            "Initial attachment cell cache requires a positive cell radius"
        );

        for (int coordinate = 0; coordinate < 3; ++coordinate) {

            ASSERT_(
                std::isfinite(position[coordinate]),
                "Initial attachment cell cache received a non-finite position"
            );
        }

        InitialAttachmentCellRecord record;

        record.cell_uid = cell_uid;
        record.position = position;
        record.radius = radius;

        const std::size_t record_index = initial_attachment_cell_records_.size();

        initial_attachment_cell_records_.push_back(record);

        initial_attachment_cell_max_radius_ =
            std::max(initial_attachment_cell_max_radius_,radius);

        const InitialAttachmentCellBucketKey bucket_key = {
            static_cast<int>(
                std::floor(
                    position[0] / initial_attachment_cell_bucket_size_
                )
            ),
            static_cast<int>(
                std::floor(
                    position[1] / initial_attachment_cell_bucket_size_
                )
            ),
            static_cast<int>(
                std::floor(
                    position[2] / initial_attachment_cell_bucket_size_
                )
            )
        };

        initial_attachment_cell_spatial_index_[bucket_key].push_back(
            record_index
        );
    }

    inline
    InitialAttachmentCellOverlapResult
    CheckInitialAttachmentCellOverlap(
        const unsigned int candidate_cell_uid,
        const bdm::Double3& proposed_position,
        const double candidate_radius,
        const double additional_clearance) const {

        /*
        * Determine if cells are overlapping with nearby cells after repostitioning for 
        * scaffold attachment.
        */

        InitialAttachmentCellOverlapResult result;

        ASSERT_(
            std::isfinite(candidate_radius) && candidate_radius > 0.0,
            "Initial attachment cell-overlap check requires a positive candidate radius"
        );

        ASSERT_(
            std::isfinite(additional_clearance) && additional_clearance >= 0.0,
            "Initial attachment cell-overlap check requires non-negative clearance"
        );

        if (initial_attachment_cell_records_.empty()) {
            return result;
        }

        ASSERT_(
            initial_attachment_cell_bucket_size_ > 0.0,
            "Initial attachment cell cache has invalid bucket size"
        );

        const double search_radius =
            candidate_radius + initial_attachment_cell_max_radius_ + additional_clearance;

        const int min_bucket_x =
            static_cast<int>(
                std::floor(
                    (proposed_position[0] - search_radius) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int max_bucket_x =
            static_cast<int>(
                std::floor(
                    (proposed_position[0] + search_radius) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int min_bucket_y =
            static_cast<int>(
                std::floor(
                    (proposed_position[1] - search_radius) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int max_bucket_y =
            static_cast<int>(
                std::floor(
                    (proposed_position[1] + search_radius) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int min_bucket_z =
            static_cast<int>(
                std::floor(
                    (proposed_position[2] - search_radius) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int max_bucket_z =
            static_cast<int>(
                std::floor(
                    (proposed_position[2] + search_radius) /
                    initial_attachment_cell_bucket_size_
                )
            );

        double smallest_margin = std::numeric_limits<double>::max();

        for (int bucket_x = min_bucket_x; bucket_x <= max_bucket_x; ++bucket_x) {
            for (int bucket_y = min_bucket_y; bucket_y <= max_bucket_y; ++bucket_y) {
                for (int bucket_z = min_bucket_z; bucket_z <= max_bucket_z; ++bucket_z) {

                    const InitialAttachmentCellBucketKey key = {
                        bucket_x,
                        bucket_y,
                        bucket_z
                    };

                    const auto bucket_it = initial_attachment_cell_spatial_index_.find(key);

                    if (bucket_it == initial_attachment_cell_spatial_index_.end()) {
                        continue;
                    }

                    for (const std::size_t record_index : bucket_it->second) {

                        ASSERT_(record_index < initial_attachment_cell_records_.size(),
                            "Initial attachment cell cache contains an invalid record index"
                        );

                        const auto& existing_cell =
                            initial_attachment_cell_records_[
                                record_index
                            ];

                        if (existing_cell.cell_uid == candidate_cell_uid) {
                            continue;
                        }

                        const double centre_distance = L2norm(
                            proposed_position - existing_cell.position
                        );

                        const double required_distance =
                            candidate_radius + existing_cell.radius + additional_clearance;

                        const double clearance_margin =
                            centre_distance - required_distance;

                        if (clearance_margin < smallest_margin) {

                            smallest_margin = clearance_margin;

                            result.blocking_cell_uid = existing_cell.cell_uid;

                            result.centre_distance = centre_distance;

                            result.required_distance = required_distance;

                            result.clearance_margin = clearance_margin;
                        }
                    }
                }
            }
        }

        if (smallest_margin == std::numeric_limits<double>::max()) {
            return result;
        }

        const double tolerance = 1.0e-12 * std::max(1.0,result.required_distance);
        result.overlaps = result.clearance_margin < -tolerance;

        return result;
    }

    inline
    InitialAttachmentCellPathObstructionResult
    CheckInitialAttachmentCellPathObstruction(
        const unsigned int candidate_cell_uid,
        const bdm::Double3& path_start,
        const bdm::Double3& path_end,
        const double additional_clearance) const
    {
        InitialAttachmentCellPathObstructionResult result;

        ASSERT_(
            initial_attachment_cell_bucket_size_ > 0.0,
            "Initial attachment cell-path obstruction requires "
            "an initialised cell cache"
        );

        ASSERT_(
            std::isfinite(additional_clearance) &&
            additional_clearance >= 0.0,
            "Initial attachment cell-path obstruction requires "
            "non-negative finite clearance"
        );

        const bdm::Double3 path_vector =
            path_end - path_start;

        const double path_length =
            L2norm(path_vector);

        ASSERT_(
            std::isfinite(path_length) &&
            path_length > 0.0,
            "Initial attachment cell-path obstruction requires "
            "a positive path length"
        );

        // ---------------------------------------------------------------------
        // Step 1: Define a bounding box around the complete candidate path
        // ---------------------------------------------------------------------

        const double search_padding =
            initial_attachment_cell_max_radius_ +
            additional_clearance;

        const int min_x =
            static_cast<int>(
                std::floor(
                    (
                        std::min(path_start[0], path_end[0]) -
                        search_padding
                    ) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int max_x =
            static_cast<int>(
                std::floor(
                    (
                        std::max(path_start[0], path_end[0]) +
                        search_padding
                    ) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int min_y =
            static_cast<int>(
                std::floor(
                    (
                        std::min(path_start[1], path_end[1]) -
                        search_padding
                    ) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int max_y =
            static_cast<int>(
                std::floor(
                    (
                        std::max(path_start[1], path_end[1]) +
                        search_padding
                    ) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int min_z =
            static_cast<int>(
                std::floor(
                    (
                        std::min(path_start[2], path_end[2]) -
                        search_padding
                    ) /
                    initial_attachment_cell_bucket_size_
                )
            );

        const int max_z =
            static_cast<int>(
                std::floor(
                    (
                        std::max(path_start[2], path_end[2]) +
                        search_padding
                    ) /
                    initial_attachment_cell_bucket_size_
                )
            );

        // ---------------------------------------------------------------------
        // Step 2: Check only cached cells in nearby spatial buckets
        // ---------------------------------------------------------------------

        for (int x = min_x; x <= max_x; ++x) {
            for (int y = min_y; y <= max_y; ++y) {
                for (int z = min_z; z <= max_z; ++z) {

                    const InitialAttachmentCellBucketKey bucket_key = {
                        x,
                        y,
                        z
                    };

                    const auto bucket_it =
                        initial_attachment_cell_spatial_index_.find(
                            bucket_key
                        );

                    if (bucket_it ==
                        initial_attachment_cell_spatial_index_.end()) {

                        continue;
                    }

                    for (const std::size_t record_index :
                        bucket_it->second) {

                        ASSERT_(
                            record_index <
                                initial_attachment_cell_records_.size(),
                            "Initial attachment cell cache contains "
                            "an invalid record index"
                        );

                        const InitialAttachmentCellRecord& existing_cell =
                            initial_attachment_cell_records_[
                                record_index
                            ];

                        if (existing_cell.cell_uid ==
                            candidate_cell_uid) {

                            continue;
                        }

                        // ---------------------------------------------------------
                        // Exact point-to-segment distance
                        // ---------------------------------------------------------

                        const bdm::Double3 start_to_cell =
                            existing_cell.position -
                            path_start;

                        const double path_length_squared =
                            path_vector[0] * path_vector[0] +
                            path_vector[1] * path_vector[1] +
                            path_vector[2] * path_vector[2];

                        ASSERT_(
                            path_length_squared > 0.0,
                            "Initial attachment candidate path "
                            "has zero squared length"
                        );

                        double path_fraction =
                            (
                                start_to_cell[0] * path_vector[0] +
                                start_to_cell[1] * path_vector[1] +
                                start_to_cell[2] * path_vector[2]
                            ) /
                            path_length_squared;

                        path_fraction =
                            std::max(
                                0.0,
                                std::min(
                                    1.0,
                                    path_fraction
                                )
                            );

                        const bdm::Double3 closest_point =
                            path_start +
                            path_vector *
                            path_fraction;

                        const double path_distance =
                            L2norm(
                                existing_cell.position -
                                closest_point
                            );

                        const double required_distance =
                            existing_cell.radius +
                            additional_clearance;

                        const double tolerance =
                            1.0e-12 *
                            std::max(
                                1.0,
                                required_distance
                            );

                        if (path_distance <
                            required_distance - tolerance) {

                            result.obstructed =
                                true;

                            result.blocking_cell_uid =
                                existing_cell.cell_uid;

                            result.path_distance =
                                path_distance;

                            result.required_distance =
                                required_distance;

                            return result;
                        }
                    }
                }
            }
        }

        return result;
    }
    
    inline
    double CalculateInitialAttachmentSegmentDistance(
        const bdm::Double3& segment_1_start,
        const bdm::Double3& segment_1_end,
        const bdm::Double3& segment_2_start,
        const bdm::Double3& segment_2_end) const
    {
        const bdm::Double3 u =
            segment_1_end - segment_1_start;

        const bdm::Double3 v =
            segment_2_end - segment_2_start;

        const bdm::Double3 w =
            segment_1_start - segment_2_start;

        const auto dot =
            [](const bdm::Double3& a,
            const bdm::Double3& b) {

                return
                    a[0] * b[0] +
                    a[1] * b[1] +
                    a[2] * b[2];
            };

        const double a = dot(u, u);
        const double b = dot(u, v);
        const double c = dot(v, v);
        const double d = dot(u, w);
        const double e = dot(v, w);

        const double epsilon = 1.0e-12;

        // Both segments collapse to points.
        if (a <= epsilon &&
            c <= epsilon) {

            return L2norm(
                segment_1_start -
                segment_2_start
            );
        }

        double s = 0.0;
        double t = 0.0;

        // First segment collapses to a point.
        if (a <= epsilon) {

            t =
                std::max(
                    0.0,
                    std::min(
                        1.0,
                        e / c
                    )
                );

        // Second segment collapses to a point.
        } else if (c <= epsilon) {

            s =
                std::max(
                    0.0,
                    std::min(
                        1.0,
                        -d / a
                    )
                );

        } else {

            const double denominator =
                a * c -
                b * b;

            if (std::abs(denominator) >
                epsilon) {

                s =
                    std::max(
                        0.0,
                        std::min(
                            1.0,
                            (
                                b * e -
                                c * d
                            ) /
                            denominator
                        )
                    );
            }

            t =
                (
                    b * s +
                    e
                ) /
                c;

            if (t < 0.0) {

                t = 0.0;

                s =
                    std::max(
                        0.0,
                        std::min(
                            1.0,
                            -d / a
                        )
                    );

            } else if (t > 1.0) {

                t = 1.0;

                s =
                    std::max(
                        0.0,
                        std::min(
                            1.0,
                            (
                                b - d
                            ) /
                            a
                        )
                    );
            }
        }

        const bdm::Double3 closest_1 =
            segment_1_start +
            u * s;

        const bdm::Double3 closest_2 =
            segment_2_start +
            v * t;

        const double distance =
            L2norm(
                closest_1 -
                closest_2
            );

        ASSERT_(
            std::isfinite(distance) &&
            distance >= 0.0,
            "Initial attachment segment-distance "
            "calculation produced an invalid distance"
        );

        return distance;
    }
    
    inline
    bool CheckInitialAdditionalAttachmentScaffoldPathObstruction(
        const bdm::Double3& path_start,
        const int candidate_node_id,
        const ObstacleScaffold& active_scaffold,
        const double additional_clearance) const
    {
        ASSERT_(
            candidate_node_id > 0,
            "Additional initial attachment scaffold-path check "
            "received a non-positive node ID"
        );

        ASSERT_(
            active_scaffold.HasNode(candidate_node_id),
            "Additional initial attachment scaffold-path check "
            "could not find candidate node"
        );

        ASSERT_(
            std::isfinite(additional_clearance) &&
            additional_clearance >= 0.0,
            "Additional initial attachment scaffold clearance "
            "must be finite and non-negative"
        );

        const bdm::Double3& path_end =
            active_scaffold.GetNodePosition(
                candidate_node_id
            );

        const double path_length =
            L2norm(
                path_end -
                path_start
            );

        ASSERT_(
            std::isfinite(path_length) &&
            path_length > 0.0,
            "Additional initial attachment path "
            "must have positive length"
        );

        const bdm::Double3 path_midpoint =
            (path_start + path_end) * 0.5;

        const double search_radius =
            0.5 * path_length +
            additional_clearance;

        const std::vector<int> nearby_segment_ids =
            active_scaffold.GetNearbySegmentIds(
                path_midpoint,
                search_radius
            );

        for (const int element_id :
            nearby_segment_ids) {

            ASSERT_(
                active_scaffold.HasSegment(element_id),
                "Additional initial attachment scaffold-path check "
                "could not find nearby segment"
            );

            const auto& scaffold_segment =
                active_scaffold.GetSegment(
                    element_id
                );

            // The candidate lies on these segments by definition.
            if (scaffold_segment.node_id_1 ==
                    candidate_node_id ||
                scaffold_segment.node_id_2 ==
                    candidate_node_id) {

                continue;
            }

            const double path_to_segment_distance =
                this->CalculateInitialAttachmentSegmentDistance(
                    path_start,
                    path_end,
                    scaffold_segment.vertex_0,
                    scaffold_segment.vertex_1
                );

            const double required_clearance =
                scaffold_segment.radius +
                additional_clearance;

            const double tolerance =
                1.0e-12 *
                std::max(
                    1.0,
                    required_clearance
                );

            if (path_to_segment_distance <
                required_clearance - tolerance) {

                return true;
            }
        }

        return false;
    }
    
    inline
    InitialAttachmentRelocationCheckResult
    GenerateValidInitialAttachmentRelocationCandidate(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const double initial_search_radius,
        const double scaffold_clearance,
        const double cell_clearance) const
    {
        InitialAttachmentRelocationCheckResult result;

        // ---------------------------------------------------------------------
        // Step 1: Generate the direct minimum-movement relocation candidate
        // ---------------------------------------------------------------------

        result.relocation =
            this->GenerateInitialAttachmentRelocationCandidate(
                cell,
                active_scaffold,
                initial_search_radius,
                scaffold_clearance
            );

        if (!result.relocation.found) {
            return result;
        }

        // ---------------------------------------------------------------------
        // Step 2: Reject positions that overlap any scaffold segment
        // ---------------------------------------------------------------------

        result.scaffold_overlap =
            this->CheckInitialAttachmentScaffoldOverlap(
                cell,
                active_scaffold,
                result.relocation.proposed_position,
                scaffold_clearance
            );

        if (result.scaffold_overlap.overlaps) {
            return result;
        }

        // ---------------------------------------------------------------------
        // Step 3: Reject positions that overlap an already accepted cell
        // ---------------------------------------------------------------------

        const double cell_radius =
            0.5 *
            cell.GetDiameter();

        ASSERT_(
            std::isfinite(cell_radius) &&
            cell_radius > 0.0,
            "Initial attachment relocation requires a positive cell radius"
        );

        result.cell_overlap =
            this->CheckInitialAttachmentCellOverlap(
                static_cast<unsigned int>(
                    cell.GetUid().GetIndex()
                ),
                result.relocation.proposed_position,
                cell_radius,
                cell_clearance
            );

        if (result.cell_overlap.overlaps) {
            return result;
        }

        // ---------------------------------------------------------------------
        // Step 4: Candidate clears both tests
        // ---------------------------------------------------------------------

        result.accepted =
            true;

        return result;
    }

    inline
    std::vector<bdm::Double3>
    GenerateInitialAttachmentRelocationPositions(
        const ObstacleScaffold& active_scaffold,
        const InitialAttachmentRelocationCandidate& relocation,
        const std::size_t number_of_positions) const
    {
        ASSERT_(
            relocation.found,
            "Alternative relocation generation requires "
            "a valid relocation candidate"
        );

        ASSERT_(
            relocation.segment_id > 0 &&
            active_scaffold.HasSegment(
                relocation.segment_id
            ),
            "Alternative relocation generation could not find "
            "the selected scaffold segment"
        );

        ASSERT_(
            number_of_positions > 0,
            "Alternative relocation generation requires "
            "at least one position"
        );

        ASSERT_(
            std::isfinite(
                relocation.target_centreline_distance
            ) &&
            relocation.target_centreline_distance > 0.0,
            "Alternative relocation generation requires "
            "a positive target centreline distance"
        );

        std::vector<bdm::Double3> positions;

        positions.reserve(
            number_of_positions
        );

        // ---------------------------------------------------------------------
        // Step 1: Always try the minimum-movement position first
        // ---------------------------------------------------------------------

        positions.push_back(
            relocation.proposed_position
        );

        if (number_of_positions == 1) {
            return positions;
        }

        // ---------------------------------------------------------------------
        // Step 2: Construct the selected strut axis
        // ---------------------------------------------------------------------

        const auto& segment =
            active_scaffold.GetSegment(
                relocation.segment_id
            );

        bdm::Double3 axis =
            segment.vertex_1 -
            segment.vertex_0;

        const double axis_length =
            L2norm(
                axis
            );

        ASSERT_(
            std::isfinite(axis_length) &&
            axis_length > 0.0,
            "Alternative relocation generation encountered "
            "a zero-length scaffold segment"
        );

        axis *=
            1.0 /
            axis_length;

        // ---------------------------------------------------------------------
        // Step 3: Orient the ring towards the direct relocation direction
        // ---------------------------------------------------------------------

        const bdm::Double3 direct_radial =
            relocation.proposed_position -
            relocation.closest_point;

        const double direct_axis_component =
            direct_radial[0] * axis[0] +
            direct_radial[1] * axis[1] +
            direct_radial[2] * axis[2];

        bdm::Double3 basis_1 =
            direct_radial -
            axis *
            direct_axis_component;

        double basis_1_length =
            L2norm(
                basis_1
            );

        /*
        * If the direct direction is parallel to the strut axis, construct a
        * deterministic perpendicular direction instead.
        */
        if (basis_1_length <= 1.0e-12) {

            const bdm::Double3 reference =
                std::fabs(axis[0]) < 0.9
                    ? bdm::Double3{1.0, 0.0, 0.0}
                    : bdm::Double3{0.0, 1.0, 0.0};

            basis_1 = {
                axis[1] * reference[2] -
                    axis[2] * reference[1],

                axis[2] * reference[0] -
                    axis[0] * reference[2],

                axis[0] * reference[1] -
                    axis[1] * reference[0]
            };

            basis_1_length =
                L2norm(
                    basis_1
                );

            ASSERT_(
                basis_1_length > 0.0,
                "Alternative relocation generation could not "
                "construct a perpendicular basis"
            );
        }

        basis_1 *=
            1.0 /
            basis_1_length;

        bdm::Double3 basis_2 = {
            axis[1] * basis_1[2] -
                axis[2] * basis_1[1],

            axis[2] * basis_1[0] -
                axis[0] * basis_1[2],

            axis[0] * basis_1[1] -
                axis[1] * basis_1[0]
        };

        const double basis_2_length =
            L2norm(
                basis_2
            );

        ASSERT_(
            basis_2_length > 0.0,
            "Alternative relocation generation produced "
            "an invalid second basis vector"
        );

        basis_2 *=
            1.0 /
            basis_2_length;

        // ---------------------------------------------------------------------
        // Step 4: Generate evenly distributed alternatives around the strut
        // ---------------------------------------------------------------------

        const std::size_t alternative_count =
            number_of_positions - 1;

        const double pi =
            std::acos(
                -1.0
            );

        for (std::size_t alternative = 0;
            alternative < alternative_count;
            ++alternative) {

            /*
            * Half-step offset prevents the first ring position from duplicating
            * the direct minimum-movement position.
            */
            const double angle =
                2.0 *
                pi *
                (
                    static_cast<double>(alternative) +
                    0.5
                ) /
                static_cast<double>(
                    alternative_count
                );

            const bdm::Double3 radial_direction =
                basis_1 *
                    std::cos(angle) +
                basis_2 *
                    std::sin(angle);

            const bdm::Double3 candidate_position =
                relocation.closest_point +
                radial_direction *
                    relocation.target_centreline_distance;

            positions.push_back(
                candidate_position
            );
        }

        ASSERT_(
            positions.size() ==
                number_of_positions,
            "Alternative relocation generation produced "
            "an incorrect number of positions"
        );

        return positions;
    }
    
    inline
    InitialAttachmentRelocationSelectionResult
    SelectValidInitialAttachmentRelocationPosition(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const std::vector<bdm::Double3>& candidate_positions,
        const double scaffold_clearance,
        const double cell_clearance) const
    {
        InitialAttachmentRelocationSelectionResult result;

        ASSERT_(
            !candidate_positions.empty(),
            "Initial attachment relocation selection requires "
            "at least one candidate position"
        );

        const double cell_radius =
            0.5 *
            cell.GetDiameter();

        ASSERT_(
            std::isfinite(cell_radius) &&
            cell_radius > 0.0,
            "Initial attachment relocation selection requires "
            "a positive cell radius"
        );

        // ---------------------------------------------------------------------
        // Test candidates in deterministic order:
        // direct position first, then rotated alternatives.
        // ---------------------------------------------------------------------

        for (std::size_t position_index = 0;
            position_index < candidate_positions.size();
            ++position_index) {

            const bdm::Double3& candidate_position =
                candidate_positions[position_index];

            const auto scaffold_overlap =
                this->CheckInitialAttachmentScaffoldOverlap(
                    cell,
                    active_scaffold,
                    candidate_position,
                    scaffold_clearance
                );

            if (scaffold_overlap.overlaps) {
                continue;
            }

            const auto cell_overlap =
                this->CheckInitialAttachmentCellOverlap(
                    static_cast<unsigned int>(
                        cell.GetUid().GetIndex()
                    ),
                    candidate_position,
                    cell_radius,
                    cell_clearance
                );

            if (cell_overlap.overlaps) {
                continue;
            }

            result.accepted =
                true;

            result.selected_position_index =
                position_index;

            result.selected_position =
                candidate_position;

            result.scaffold_overlap =
                scaffold_overlap;

            result.cell_overlap =
                cell_overlap;

            return result;
        }

        return result;
    }
    
    inline
    InitialAttachmentRelocationFallbackResult
    SelectValidInitialAttachmentRelocationAcrossSegments(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const std::vector<int>& ordered_segment_ids,
        const std::size_t positions_per_segment,
        const double scaffold_clearance,
        const double cell_clearance) const
    {
        InitialAttachmentRelocationFallbackResult result;

        ASSERT_(
            !ordered_segment_ids.empty(),
            "Initial attachment relocation fallback "
            "requires at least one scaffold segment"
        );

        ASSERT_(
            positions_per_segment > 0,
            "Initial attachment relocation fallback "
            "requires at least one position per segment"
        );

        // ---------------------------------------------------------------------
        // Try scaffold segments from nearest to furthest
        // ---------------------------------------------------------------------

        for (std::size_t segment_index = 0;
            segment_index < ordered_segment_ids.size();
            ++segment_index) {

            const int segment_id =
                ordered_segment_ids[
                    segment_index
                ];

            ++result.attempted_segment_count;

            const auto relocation =
                this->GenerateInitialAttachmentRelocationCandidateForSegment(
                    cell,
                    active_scaffold,
                    segment_id,
                    scaffold_clearance
                );

            ASSERT_(
                relocation.found,
                "Initial attachment relocation fallback "
                "could not generate a segment candidate"
            );

            const auto candidate_positions =
                this->GenerateInitialAttachmentRelocationPositions(
                    active_scaffold,
                    relocation,
                    positions_per_segment
                );

            const auto selection =
                this->SelectValidInitialAttachmentRelocationPosition(
                    cell,
                    active_scaffold,
                    candidate_positions,
                    scaffold_clearance,
                    cell_clearance
                );

            if (!selection.accepted) {
                continue;
            }

            result.accepted =
                true;

            result.selected_segment_index =
                segment_index;

            result.selected_segment_id =
                segment_id;

            result.selected_position_index =
                selection.selected_position_index;

            result.selected_position =
                selection.selected_position;

            return result;
        }

        return result;
    }
    
    inline
    std::vector<int>
    GenerateInitialAdditionalAttachmentCandidateNodeIds(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const double min_attachment_separation,
        const double max_attachment_separation,
        const double candidate_scaffold_clearance,
        const double candidate_cell_clearance) const
    {
        ASSERT_(
            cell.GetCellMatrixLifecycleStatus() ==
                bdm::BiologicalCell::
                    CellMatrixLifecycleStatus::kEstablished,
            "Initial additional-attachment candidate generation "
            "requires an established cell"
        );

        ASSERT_(
            cell.HasAttachmentRecords(),
            "Initial additional-attachment candidate generation "
            "requires existing attachments"
        );

        ASSERT_(
            std::isfinite(min_attachment_separation) &&
            min_attachment_separation >= 0.0,
            "Minimum attachment separation must be finite "
            "and non-negative"
        );

        ASSERT_(
            std::isfinite(max_attachment_separation) &&
            max_attachment_separation > 0.0 &&
            min_attachment_separation <=
                max_attachment_separation,
            "Initial additional-attachment candidate generation "
            "received invalid separation limits"
        );

        // ---------------------------------------------------------------------
        // Step 1: Select one deterministic retained attachment for the
        //         scaffold spatial query
        // ---------------------------------------------------------------------

        const auto& attachment_records =
            cell.GetAttachmentRecords();

        const auto search_attachment_it =
            std::min_element(
                attachment_records.begin(),
                attachment_records.end(),
                [](
                    const bdm::BiologicalCell::AttachmentRecord& a,
                    const bdm::BiologicalCell::AttachmentRecord& b) {

                    return a.node_id < b.node_id;
                }
            );

        ASSERT_(
            search_attachment_it !=
                attachment_records.end(),
            "Could not select a retained attachment "
            "for additional-attachment search"
        );

        const int search_node_id =
            search_attachment_it->node_id;

        ASSERT_(
            active_scaffold.HasNode(
                search_node_id
            ),
            "Additional-attachment search could not find "
            "retained scaffold node ID "
            + std::to_string(search_node_id)
        );

        const bdm::Double3& search_centre =
            active_scaffold.GetNodePosition(
                search_node_id
            );

        // ---------------------------------------------------------------------
        // Step 2: Retrieve only nearby scaffold nodes
        // ---------------------------------------------------------------------

        const std::vector<int> nearby_node_ids =
            active_scaffold.GetNodeIdsWithinRadius(
                search_centre,
                max_attachment_separation
            );

        std::unordered_set<int> attached_node_ids;

        for (const auto& attachment :
            attachment_records) {

            attached_node_ids.insert(
                attachment.node_id
            );
        }

        std::vector<int> candidate_node_ids;

        // ---------------------------------------------------------------------
        // Step 3: Apply exact candidate filters
        // ---------------------------------------------------------------------

        for (const int candidate_node_id :
            nearby_node_ids) {

            if (attached_node_ids.count(
                    candidate_node_id
                ) > 0) {

                continue;
            }

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

            // -------------------------------------------------------------
            // Pairwise separation from every retained attachment
            // -------------------------------------------------------------

            bool separation_valid =
                true;

            for (const auto& retained_attachment :
                attachment_records) {

                const bdm::Double3& retained_position =
                    active_scaffold.GetNodePosition(
                        retained_attachment.node_id
                    );

                const double separation =
                    L2norm(
                        candidate_position -
                        retained_position
                    );

                if (separation <
                        min_attachment_separation ||
                    separation >
                        max_attachment_separation) {

                    separation_valid =
                        false;

                    break;
                }
            }

            if (!separation_valid) {
                continue;
            }

            // -------------------------------------------------------------
            // Scaffold obstruction
            // -------------------------------------------------------------

            if (this->CheckInitialAdditionalAttachmentScaffoldPathObstruction(
                    cell.GetPosition(),
                    candidate_node_id,
                    active_scaffold,
                    candidate_scaffold_clearance)) {

                continue;
            }

            // -------------------------------------------------------------
            // Cell obstruction using the scalable initialisation cache
            // -------------------------------------------------------------

            const auto cell_obstruction =
                this->CheckInitialAttachmentCellPathObstruction(
                    static_cast<unsigned int>(
                        cell.GetUid().GetIndex()
                    ),
                    cell.GetPosition(),
                    candidate_position,
                    candidate_cell_clearance
                );

            if (cell_obstruction.obstructed) {
                continue;
            }

            candidate_node_ids.push_back(
                candidate_node_id
            );
        }

        // Preserve deterministic ordering before weighting/selection.
        std::sort(
            candidate_node_ids.begin(),
            candidate_node_ids.end()
        );

        return candidate_node_ids;
    }
    
    inline
    InitialAdditionalAttachmentSelectionResult
    SelectInitialAdditionalAttachment(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const std::vector<int>& candidate_node_ids,
        const double distance_sensitivity,
        const double random_selection_strength,
        const unsigned int draw_index) const
    {
        InitialAdditionalAttachmentSelectionResult result;

        if (candidate_node_ids.empty()) {
            return result;
        }

        ASSERT_(
            std::isfinite(distance_sensitivity) &&
            distance_sensitivity >= 0.0,
            "Additional initial attachment distance sensitivity "
            "must be finite and non-negative"
        );

        ASSERT_(
            std::isfinite(random_selection_strength) &&
            random_selection_strength >= 0.0 &&
            random_selection_strength <= 1.0,
            "Additional initial attachment random-selection strength "
            "must lie within [0, 1]"
        );

        // ---------------------------------------------------------------------
        // Step 1: Calculate candidate distances
        // ---------------------------------------------------------------------

        result.candidate_distances.reserve(
            candidate_node_ids.size()
        );

        for (const int node_id :
            candidate_node_ids) {

            ASSERT_(
                active_scaffold.HasNode(node_id),
                "Additional attachment selection could not find candidate node"
            );

            const double distance =
                L2norm(
                    active_scaffold.GetNodePosition(node_id) -
                    cell.GetPosition()
                );

            ASSERT_(
                std::isfinite(distance) &&
                distance >= 0.0,
                "Additional attachment selection calculated "
                "an invalid candidate distance"
            );

            result.candidate_distances.push_back(
                distance
            );
        }

        const auto distance_limits =
            std::minmax_element(
                result.candidate_distances.begin(),
                result.candidate_distances.end()
            );

        const double minimum_distance =
            *distance_limits.first;

        const double maximum_distance =
            *distance_limits.second;

        const double distance_range =
            maximum_distance -
            minimum_distance;

        // ---------------------------------------------------------------------
        // Step 2: Convert distance into bias weights
        // ---------------------------------------------------------------------

        result.candidate_weights.reserve(
            candidate_node_ids.size()
        );

        for (const double distance :
            result.candidate_distances) {

            double normalised_distance =
                0.0;

            if (distance_range > 0.0) {

                normalised_distance =
                    (
                        distance -
                        minimum_distance
                    ) /
                    distance_range;
            }

            const double weight =
                std::exp(
                    -distance_sensitivity *
                    normalised_distance
                );

            ASSERT_(
                std::isfinite(weight) &&
                weight > 0.0 &&
                weight <= 1.0,
                "Additional initial attachment produced "
                "an invalid distance weight"
            );

            result.candidate_weights.push_back(
                weight
            );
        }

        // ---------------------------------------------------------------------
        // Step 3: Convert weights into selection probabilities
        // ---------------------------------------------------------------------

        result.selection_probabilities.assign(
            candidate_node_ids.size(),
            0.0
        );

        if (random_selection_strength == 0.0) {

            const auto highest_weight =
                std::max_element(
                    result.candidate_weights.begin(),
                    result.candidate_weights.end()
                );

            const std::size_t selected_index =
                static_cast<std::size_t>(
                    std::distance(
                        result.candidate_weights.begin(),
                        highest_weight
                    )
                );

            result.selection_probabilities[
                selected_index
            ] = 1.0;

        } else {

            const double gamma =
                (
                    1.0 -
                    random_selection_strength
                ) /
                random_selection_strength;

            std::vector<double> log_scores;

            log_scores.reserve(
                candidate_node_ids.size()
            );

            for (const double weight :
                result.candidate_weights) {

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

            double total_score =
                0.0;

            for (std::size_t index = 0;
                index < log_scores.size();
                ++index) {

                const double score =
                    std::exp(
                        log_scores[index] -
                        maximum_log_score
                    );

                result.selection_probabilities[
                    index
                ] = score;

                total_score += score;
            }

            ASSERT_(
                std::isfinite(total_score) &&
                total_score > 0.0,
                "Additional initial attachment produced "
                "an invalid probability total"
            );

            for (double& probability :
                result.selection_probabilities) {

                probability /=
                    total_score;
            }
        }

        double probability_sum =
            0.0;

        for (const double probability :
            result.selection_probabilities) {

            ASSERT_(
                std::isfinite(probability) &&
                probability >= 0.0 &&
                probability <= 1.0,
                "Additional initial attachment probability "
                "is outside [0, 1]"
            );

            probability_sum +=
                probability;
        }

        ASSERT_(
            std::fabs(
                probability_sum -
                1.0
            ) <= 1.0e-12,
            "Additional initial attachment probabilities "
            "do not sum to 1"
        );

        // ---------------------------------------------------------------------
        // Step 4: Reproducible stochastic selection
        // ---------------------------------------------------------------------

        const std::string& phenotype_name =
            this->params()->get<std::string>(
                "phenotype_ID/" +
                std::to_string(
                    cell.GetPhenotype()
                )
            );

        const int random_state =
            this->params()->get<int>(
                phenotype_name +
                "/cell_matrix_mechanics/random_state"
            );

        std::seed_seq seed {
            static_cast<unsigned int>(
                random_state
            ),
            static_cast<unsigned int>(
                cell.GetUid().GetIndex()
            ),
            draw_index,
            7303u
        };

        std::mt19937 generator(seed);

        std::uniform_real_distribution<double>
            distribution(
                0.0,
                1.0
            );

        const double random_draw =
            distribution(generator);

        double cumulative_probability =
            0.0;

        std::size_t selected_index =
            candidate_node_ids.size() - 1;

        for (std::size_t index = 0;
            index <
                result.selection_probabilities.size();
            ++index) {

            cumulative_probability +=
                result.selection_probabilities[
                    index
                ];

            if (random_draw <
                cumulative_probability) {

                selected_index =
                    index;

                break;
            }
        }

        result.selected_node_id =
            candidate_node_ids[
                selected_index
            ];

        result.success =
            true;

        return result;
    }
    
    inline
    bool CommitInitialAdditionalAttachment(
        bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const int selected_node_id) 
    {
        ASSERT_(
            selected_node_id > 0,
            "Initial additional attachment commit "
            "received a non-positive node ID"
        );

        ASSERT_(
            active_scaffold.HasNode(
                selected_node_id
            ),
            "Initial additional attachment commit "
            "could not find selected scaffold node"
        );

        const std::size_t attachment_count_before =
            cell.GetNumberOfAttachmentRecords();

        ASSERT_(
            attachment_count_before >= 2,
            "Initial additional attachment commit "
            "requires an established initial pair"
        );

        for (const auto& attachment :
            cell.GetAttachmentRecords()) {

            ASSERT_(
                attachment.node_id !=
                    selected_node_id,
                "Initial additional attachment commit "
                "received an already attached node"
            );
        }

        // Store the existing attachment node IDs before the new node is committed.
        // Each of these will form one new pair with selected_node_id.
        std::vector<int> previous_attachment_node_ids;

        previous_attachment_node_ids.reserve(
            cell.GetNumberOfAttachmentRecords()
        );

        for (const auto& attachment :
            cell.GetAttachmentRecords()) {

            previous_attachment_node_ids.push_back(
                attachment.node_id
            );
        }

        cell.AddInitialSelectedAttachment(
            active_scaffold,
            selected_node_id
        );

        // ---------------------------------------------------------------------
        // Add only the newly created attachment-pair paths to the cache
        // ---------------------------------------------------------------------

        for (const int previous_node_id :
            previous_attachment_node_ids) {

            ASSERT_(
                active_scaffold.HasNode(
                    previous_node_id
                ),
                "Initial additional attachment cache update "
                "could not find previous scaffold node"
            );

            this->AddInitialAttachmentPathToCache(
                static_cast<unsigned int>(
                    cell.GetUid().GetIndex()
                ),
                active_scaffold.GetNodePosition(
                    previous_node_id
                ),
                active_scaffold.GetNodePosition(
                    selected_node_id
                )
            );
        }

        ASSERT_(
            cell.GetNumberOfAttachmentRecords() ==
                attachment_count_before + 1,
            "Initial additional attachment commit "
            "did not add exactly one attachment"
        );

        ASSERT_(
            cell.RequiresMechanicsRecalculation(),
            "Initial additional attachment must "
            "request FEM mechanics recalculation"
        );

        cell.ValidateCellMatrixState();

        return true;
    }
    
    inline
    InitialAdditionalAttachmentCompletionResult
    CompleteInitialAdditionalAttachmentsForCell(
        bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold) 
    {
        InitialAdditionalAttachmentCompletionResult result;

        ASSERT_(
            cell.GetCellMatrixLifecycleStatus() ==
                bdm::BiologicalCell::
                    CellMatrixLifecycleStatus::kEstablished,
            "Initial additional attachment completion "
            "requires an established cell"
        );

        ASSERT_(
            cell.GetNumberOfAttachmentRecords() >= 2,
            "Initial additional attachment completion "
            "requires an initial attachment pair"
        );

        // ---------------------------------------------------------------------
        // Step 1: Read phenotype parameters
        // ---------------------------------------------------------------------

        const std::string& phenotype_name =
            this->params()->get<std::string>(
                "phenotype_ID/" +
                std::to_string(
                    cell.GetPhenotype()
                )
            );

        const std::string mech_base =
            phenotype_name +
            "/cell_matrix_mechanics";

        const int target_attachment_count_input =
            this->params()->get<int>(
                mech_base +
                "/num_attachments"
            );

        ASSERT_(
            target_attachment_count_input >= 2,
            "Initial attachment generation requires "
            "num_attachments >= 2"
        );

        const std::size_t target_attachment_count =
            static_cast<std::size_t>(
                target_attachment_count_input
            );

        const double min_attachment_separation =
            this->params()->get<double>(
                mech_base +
                "/min_attachment_separation"
            );

        const double max_attachment_separation =
            this->params()->get<double>(
                mech_base +
                "/max_attachment_separation"
            );

        const double scaffold_clearance =
            this->params()->get<double>(
                mech_base +
                "/candidate_scaffold_clearance"
            );

        const double cell_clearance =
            this->params()->get<double>(
                mech_base +
                "/candidate_cell_clearance"
            );

        const double distance_sensitivity =
            this->params()->get<double>(
                mech_base +
                "/candidate_cell_distance_sensitivity"
            );

        const double random_selection_strength =
            this->params()->get<double>(
                mech_base +
                "/random_selection_strength"
            );

        result.initial_attachment_count =
            cell.GetNumberOfAttachmentRecords();

        result.target_attachment_count =
            target_attachment_count;

        ASSERT_(
            result.initial_attachment_count <=
                target_attachment_count,
            "Initial cell already exceeds configured "
            "num_attachments"
        );

        // ---------------------------------------------------------------------
        // Step 2: Add attachments sequentially
        // ---------------------------------------------------------------------

        unsigned int draw_index = 0;

        while (
            cell.GetNumberOfAttachmentRecords() <
            target_attachment_count
        ) {

            const std::vector<int> candidate_node_ids =
                this->GenerateInitialAdditionalAttachmentCandidateNodeIds(
                    cell,
                    active_scaffold,
                    min_attachment_separation,
                    max_attachment_separation,
                    scaffold_clearance,
                    cell_clearance
                );

            if (candidate_node_ids.empty()) {

                result.stopped_no_candidates =
                    true;

                break;
            }

            const InitialAdditionalAttachmentSelectionResult
                selection =
                    this->SelectInitialAdditionalAttachment(
                        cell,
                        active_scaffold,
                        candidate_node_ids,
                        distance_sensitivity,
                        random_selection_strength,
                        draw_index
                    );

            ASSERT_(
                selection.success,
                "Initial additional attachment selection failed"
            );

            ASSERT_(
                selection.selected_node_id > 0,
                "Initial additional attachment selection "
                "returned an invalid node"
            );

            const bool committed =
                this->CommitInitialAdditionalAttachment(
                    cell,
                    active_scaffold,
                    selection.selected_node_id
                );

            ASSERT_(
                committed,
                "Initial additional attachment commit failed"
            );

            ++result.formed_attachment_count;

            ++draw_index;
        }

        // ---------------------------------------------------------------------
        // Step 3: Validate final state
        // ---------------------------------------------------------------------

        result.final_attachment_count =
            cell.GetNumberOfAttachmentRecords();

        result.reached_target =
            (
                result.final_attachment_count ==
                target_attachment_count
            );

        ASSERT_(
            result.final_attachment_count <=
                target_attachment_count,
            "Initial attachment completion exceeded "
            "the configured target"
        );

        ASSERT_(
            result.final_attachment_count ==
                result.initial_attachment_count +
                result.formed_attachment_count,
            "Initial additional attachment accounting "
            "is inconsistent"
        );

        ASSERT_(
            result.reached_target ||
            result.stopped_no_candidates,
            "Initial additional attachment completion "
            "stopped unexpectedly"
        );

        cell.ValidateCellMatrixState();

        return result;
    }
    
    inline
    std::string FormatAttachmentNodeIdsJsonLike(
        const std::vector<bdm::BiologicalCell::AttachmentRecord>& records) const
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

        for (std::size_t i = 0; i < records.size(); ++i) {
            ASSERT_(
            records[i].node_id > 0,
            "Cannot export an attachment record with a non-positive node ID"
            );

            oss << records[i].node_id;

            if (i + 1 < records.size()) {
            oss << ", ";
            }
        }

        oss << "]";

        return oss.str();
    }

    inline
    std::string FormatPairContractileForcesJsonLike(
        const std::vector<bdm::BiologicalCell::PairMechanicsRecord>& records) const {
        /*
        * Function goal
        * -------------
        * Format pair node IDs and their contractile forces as a JSON-like nested
        * list that can be parsed directly by Python. [node_a, node_b, contractile_force]
        *
        * Example: 
        * [[723, 1100, 47.2], [723, 730, 50.0]]
        */

        std::ostringstream oss;

        oss << "[";

        for (std::size_t i = 0; i < records.size(); ++i) {

            const auto& record = records[i];

            ASSERT_(
                record.node_id_a > 0 &&
                record.node_id_b > 0 &&
                record.node_id_a != record.node_id_b,
                "Cannot export invalid pair-mechanics node IDs"
            );

            ASSERT_(
                std::isfinite(record.contractile_force) &&
                record.contractile_force >= 0.0,
                "Cannot export an invalid pair contractile force"
            );

            oss
                << "[" << record.node_id_a << ", " << record.node_id_b
                << ", " << record.contractile_force << "]";

            if (i + 1 < records.size()) {
                oss << ", ";
            }
        }

        oss << "]";

        return oss.str();
    }

    inline
    FemExportAction DetermineFemExportAction(
        const bdm::BiologicalCell* cell,
        const std::map<int, std::string>& cells) const {

        /*
        * Determine whether a cell should be omitted from the FEM input, sent for
        * initial attachment, or sent for contraction.
        *
        * The decision uses the internal ABM attachment lifecycle and the persistent
        * attachment-record count. The FEM communication strings "attach" and
        * "contract" are assigned later during export.
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

        const int phenotype_id = cell->GetPhenotype();

        if (phenotype_id < 0) {
            return FemExportAction::kSkip;
        }

        // -------------------------------------------------------------------------
        // Step 3: Find the cell phenotype definition
        // -------------------------------------------------------------------------

        const auto phenotype_it = cells.find(phenotype_id);

        ASSERT_(
            phenotype_it != cells.end(),
            "DetermineFemExportAction could not find phenotype ID "
            + std::to_string(phenotype_id)
        );

        const std::string& phenotype_name = phenotype_it->second;
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

        // Confirm that lifecycle and attachment data agree before export.
        cell->ValidateCellMatrixState();

        // -------------------------------------------------------------------------
        // Step 5: Read the internal attachment lifecycle
        // -------------------------------------------------------------------------

        using LifecycleStatus =
            bdm::BiologicalCell::CellMatrixLifecycleStatus;

        const LifecycleStatus lifecycle_status =
            cell->GetCellMatrixLifecycleStatus();

        const std::size_t attachment_count =
            cell->GetNumberOfAttachmentRecords();

        // -------------------------------------------------------------------------
        // Step 6: Exclude cells marked for death
        // -------------------------------------------------------------------------

        if (lifecycle_status == LifecycleStatus::kMarkedForDeath) {
            return FemExportAction::kSkip;
        }

        // -------------------------------------------------------------------------
        // Step 7: Send new or retrying cells for initial FEM attachment
        // -------------------------------------------------------------------------

        if (lifecycle_status == LifecycleStatus::kNeedsInitialAttachment ||
            lifecycle_status == LifecycleStatus::kInitialAttachmentRetry) {

            ASSERT_(
                attachment_count == 0,
                "A cell requiring initial FEM attachment already contains persistent "
                "attachment records"
            );

            return FemExportAction::kAttach;
        }

        // -------------------------------------------------------------------------
        // Step 8: Validate the established-cell attachment state
        // -------------------------------------------------------------------------

        ASSERT_(
            lifecycle_status == LifecycleStatus::kEstablished,
            "DetermineFemExportAction encountered an unsupported cell-matrix "
            "lifecycle status"
        );

        ASSERT_(
            attachment_count > 0,
            "An established cell must retain at least one persistent attachment"
        );

        // -------------------------------------------------------------------------
        // Step 9: Keep single-attachment cells entirely within the ABM
        // -------------------------------------------------------------------------

        if (attachment_count == 1) {
            return FemExportAction::kSkip;
        }

        // -------------------------------------------------------------------------
        // Step 10: Contract using the current persistent attachment set
        //
        // Stored k_ecm values may be invalid or awaiting recalculation. The FEM only
        // requires the persistent node IDs and the carried-forward contractile force.
        // -------------------------------------------------------------------------

        return FemExportAction::kContract;
    }

    inline
    void ExportInitialAttachmentConfiguration(bdm::Simulation& sim) const{
        
        /*
        * Export the initial attachment configuration of all BiologicalCell agents.
        * This csv will later allow same cell populations to be tested with different 
        * model parameters.
        * 
        * Example file format:
        * -------------------
        * abm_cell_id,phenotype_id,x,y,z,attachment_node_ids
        * 0-0,1,-44.89,-10.35,19.76,"[723,727,1100]"
        */
        
        const std::string filename =
            this->params()->get<std::string>("output_directory")
            + "/initial_attachments.csv";

        std::ofstream fout(filename);

        ASSERT_(
            fout.good(),
            "Could not create initial attachment file \""
            + filename + "\""
        );

        fout << std::fixed << std::setprecision(6);

        fout
            << "abm_cell_id,"
            << "phenotype_id,"
            << "x,y,z,"
            << "attachment_node_ids"
            << std::endl;

        std::vector<bdm::BiologicalCell*> biological_cells;
        std::mutex biological_cells_mutex;

        auto* rm = sim.GetResourceManager();

        rm->ForEachAgent(
            [&](bdm::Agent* agent) {

                auto* cell =
                    dynamic_cast<bdm::BiologicalCell*>(
                        agent
                    );

                if (cell == nullptr) {
                    return;
                }

                if (!cell->HasAttachmentRecords()) {
                    return;
                }

                std::lock_guard<std::mutex> lock(
                    biological_cells_mutex
                );

                biological_cells.push_back(cell);
            }
        );

        std::sort(
            biological_cells.begin(),
            biological_cells.end(),
            [](const bdm::BiologicalCell* a,
            const bdm::BiologicalCell* b) {

                return a->GetUid().GetIndex()
                    < b->GetUid().GetIndex();
            }
        );

        for (const auto* cell : biological_cells) {

            std::ostringstream uid;
            uid << cell->GetUid();

            const auto& position = cell->GetPosition();

            fout
                << uid.str() << ","
                << cell->GetPhenotype() << ","
                << position[0] << ","
                << position[1] << ","
                << position[2] << ",\""
                << this->FormatAttachmentNodeIdsJsonLike(
                    cell->GetAttachmentRecords()
                )
                << "\""
                << std::endl;
        }
    }

    inline
    std::vector<int>
    GenerateInitialAttachmentRelocationSegmentIds(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const double initial_search_radius,
        const std::size_t maximum_segment_count) const
    {
        ASSERT_(
            std::isfinite(initial_search_radius) &&
            initial_search_radius > 0.0,
            "Initial attachment relocation segment search requires "
            "a positive search radius"
        );

        ASSERT_(
            maximum_segment_count > 0,
            "Initial attachment relocation segment search requires "
            "at least one requested segment"
        );

        ASSERT_(
            !active_scaffold.segment.empty(),
            "Initial attachment relocation segment search requires "
            "scaffold segments"
        );

        const bdm::Double3 cell_position =
            cell.GetPosition();

        const std::size_t target_segment_count =
            std::min(
                maximum_segment_count,
                active_scaffold.segment.size()
            );

        struct SegmentDistanceRecord {
            int segment_id = -1;
            double distance = 0.0;
        };

        std::vector<SegmentDistanceRecord>
            distance_records;

        double search_radius =
            initial_search_radius;

        constexpr std::size_t max_search_expansions =
            32;

        // ---------------------------------------------------------------------
        // Expand the local search until the requested nearest segments are
        // guaranteed to be contained inside the broad-phase search.
        // ---------------------------------------------------------------------

        for (std::size_t expansion = 0;
            expansion < max_search_expansions;
            ++expansion) {

            const std::vector<int> nearby_segment_ids =
                active_scaffold.GetNearbySegmentIds(
                    cell_position,
                    search_radius
                );

            distance_records.clear();

            distance_records.reserve(
                nearby_segment_ids.size()
            );

            for (const int segment_id :
                nearby_segment_ids) {

                ASSERT_(
                    active_scaffold.HasSegment(
                        segment_id
                    ),
                    "Initial attachment relocation segment search "
                    "could not find scaffold segment"
                );

                const auto& segment =
                    active_scaffold.GetSegment(
                        segment_id
                    );

                const bdm::Double3 segment_vector =
                    segment.vertex_1 -
                    segment.vertex_0;

                const double segment_length_squared =
                    segment_vector *
                    segment_vector;

                ASSERT_(
                    std::isfinite(segment_length_squared) &&
                    segment_length_squared > 0.0,
                    "Initial attachment relocation segment search "
                    "encountered a zero-length segment"
                );

                const bdm::Double3 point_offset =
                    cell_position -
                    segment.vertex_0;

                double projected_fraction =
                    (
                        point_offset *
                        segment_vector
                    ) /
                    segment_length_squared;

                projected_fraction =
                    std::max(
                        0.0,
                        std::min(
                            1.0,
                            projected_fraction
                        )
                    );

                const bdm::Double3 closest_point =
                    segment.vertex_0 +
                    segment_vector *
                    projected_fraction;

                const double distance =
                    L2norm(
                        cell_position -
                        closest_point
                    );

                ASSERT_(
                    std::isfinite(distance) &&
                    distance >= 0.0,
                    "Initial attachment relocation segment search "
                    "calculated an invalid distance"
                );

                distance_records.push_back({
                    segment_id,
                    distance
                });
            }

            std::sort(
                distance_records.begin(),
                distance_records.end(),
                [](
                    const SegmentDistanceRecord& a,
                    const SegmentDistanceRecord& b) {

                    if (a.distance != b.distance) {
                        return
                            a.distance <
                            b.distance;
                    }

                    return
                        a.segment_id <
                        b.segment_id;
                }
            );

            /*
            * Once the requested kth segment lies inside the current search
            * radius, no segment outside the search region can be closer.
            */
            if (distance_records.size() >=
                target_segment_count) {

                const double furthest_selected_distance =
                    distance_records[
                        target_segment_count - 1
                    ].distance;

                if (furthest_selected_distance <=
                    search_radius) {

                    break;
                }
            }

            search_radius *=
                2.0;
        }

        ASSERT_(
            distance_records.size() >=
                target_segment_count,
            "Initial attachment relocation segment search "
            "could not retrieve the requested number of segments"
        );

        std::vector<int> ordered_segment_ids;

        ordered_segment_ids.reserve(
            target_segment_count
        );

        for (std::size_t i = 0;
            i < target_segment_count;
            ++i) {

            ordered_segment_ids.push_back(
                distance_records[i].segment_id
            );
        }

        return ordered_segment_ids;
    }
    
    inline
    InitialAttachmentRelocationCandidate
    GenerateInitialAttachmentRelocationCandidateForSegment(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const int segment_id,
        const double scaffold_clearance) const
    {
        InitialAttachmentRelocationCandidate result;

        ASSERT_(
            segment_id > 0 &&
            active_scaffold.HasSegment(segment_id),
            "Initial attachment relocation could not find "
            "the requested scaffold segment"
        );

        ASSERT_(
            std::isfinite(scaffold_clearance) &&
            scaffold_clearance >= 0.0,
            "Initial attachment relocation requires "
            "non-negative scaffold clearance"
        );

        const bdm::Double3 original_position =
            cell.GetPosition();

        const double cell_radius =
            0.5 *
            cell.GetDiameter();

        ASSERT_(
            std::isfinite(cell_radius) &&
            cell_radius > 0.0,
            "Initial attachment relocation requires "
            "a positive cell radius"
        );

        const auto& segment =
            active_scaffold.GetSegment(
                segment_id
            );

        ASSERT_(
            std::isfinite(segment.radius) &&
            segment.radius >= 0.0,
            "Initial attachment relocation encountered "
            "an invalid strut radius"
        );

        const bdm::Double3 segment_vector =
            segment.vertex_1 -
            segment.vertex_0;

        const double segment_length_squared =
            segment_vector *
            segment_vector;

        ASSERT_(
            std::isfinite(segment_length_squared) &&
            segment_length_squared > 0.0,
            "Initial attachment relocation encountered "
            "a zero-length segment"
        );

        // ---------------------------------------------------------------------
        // Find the closest point on this specific segment
        // ---------------------------------------------------------------------

        const bdm::Double3 point_offset =
            original_position -
            segment.vertex_0;

        double projected_fraction =
            (
                point_offset *
                segment_vector
            ) /
            segment_length_squared;

        projected_fraction =
            std::max(
                0.0,
                std::min(
                    1.0,
                    projected_fraction
                )
            );

        const bdm::Double3 closest_point =
            segment.vertex_0 +
            segment_vector *
            projected_fraction;

        bdm::Double3 radial_direction =
            original_position -
            closest_point;

        double radial_norm =
            L2norm(
                radial_direction
            );

        const double target_centreline_distance =
            cell_radius +
            segment.radius +
            scaffold_clearance;

        // ---------------------------------------------------------------------
        // Construct a deterministic perpendicular direction if required
        // ---------------------------------------------------------------------

        if (radial_norm <= 1.0e-12) {

            bdm::Double3 axis =
                segment_vector;

            const double axis_length =
                L2norm(
                    axis
                );

            ASSERT_(
                axis_length > 0.0,
                "Initial attachment relocation encountered "
                "a zero-length segment"
            );

            axis *=
                1.0 /
                axis_length;

            const bdm::Double3 reference =
                std::fabs(axis[0]) < 0.9
                    ? bdm::Double3{1.0, 0.0, 0.0}
                    : bdm::Double3{0.0, 1.0, 0.0};

            radial_direction = {
                axis[1] * reference[2] -
                    axis[2] * reference[1],

                axis[2] * reference[0] -
                    axis[0] * reference[2],

                axis[0] * reference[1] -
                    axis[1] * reference[0]
            };

            radial_norm =
                L2norm(
                    radial_direction
                );

            ASSERT_(
                radial_norm > 0.0,
                "Initial attachment relocation could not "
                "construct a perpendicular direction"
            );
        }

        radial_direction *=
            1.0 /
            radial_norm;

        const bdm::Double3 proposed_position =
            closest_point +
            radial_direction *
            target_centreline_distance;

        result.found =
            true;

        result.segment_id =
            segment_id;

        result.closest_point =
            closest_point;

        result.proposed_position =
            proposed_position;

        result.original_centreline_distance =
            radial_norm;

        result.target_centreline_distance =
            target_centreline_distance;

        result.movement_distance =
            L2norm(
                proposed_position -
                original_position
            );

        return result;
    }
    
    inline
    InitialAttachmentRelocationCandidate
    GenerateInitialAttachmentRelocationCandidate(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const double initial_search_radius,
        const double scaffold_clearance) const {

        /*
        * Find the smallest radial movement to get a cell within a defined distance to
        * a scaffold obstacle.
        */

        InitialAttachmentRelocationCandidate result;

        ASSERT_(
            std::isfinite(initial_search_radius) && initial_search_radius > 0.0,
            "Initial attachment relocation requires a positive search radius"
        );

        ASSERT_(
            std::isfinite(scaffold_clearance) && scaffold_clearance >= 0.0,
            "Initial attachment relocation requires non-negative scaffold clearance"
        );

        ASSERT_(
            !active_scaffold.segment.empty(),
            "Initial attachment relocation requires scaffold segments"
        );

        const bdm::Double3 original_position = cell.GetPosition();

        const double cell_radius = 0.5 * cell.GetDiameter();

        ASSERT_(
            std::isfinite(cell_radius) && cell_radius > 0.0,
            "Initial attachment relocation requires a positive cell radius"
        );

        // ---------------------------------------------------------------------
        // Retrieve a local set of scaffold segments.
        // Expand the search only when no scaffold segment is found locally.
        // ---------------------------------------------------------------------

        double search_radius = initial_search_radius;

        std::vector<int> nearby_segment_ids;

        constexpr std::size_t max_search_expansions = 32;

        for (std::size_t expansion = 0; expansion < max_search_expansions; ++expansion) {

            nearby_segment_ids =
                active_scaffold.GetNearbySegmentIds(
                    original_position,
                    search_radius
                );

            if (!nearby_segment_ids.empty()) {
                break;
            }

            search_radius *= 2.0;
        }

        ASSERT_(
            !nearby_segment_ids.empty(),
            "Initial attachment relocation could not locate a nearby scaffold segment"
        );

        // ---------------------------------------------------------------------
        // Find the exact nearest segment among the broad-phase results.
        // ---------------------------------------------------------------------

        int nearest_segment_id = -1;

        bdm::Double3 nearest_closest_point;

        double nearest_distance = std::numeric_limits<double>::max();

        auto evaluate_segments = [&](const std::vector<int>& segment_ids) {
                for (const int segment_id : segment_ids) {

                    ASSERT_(
                        active_scaffold.HasSegment(segment_id),
                        "Initial attachment relocation could not find scaffold segment"
                    );

                    const auto& segment = active_scaffold.GetSegment(segment_id);

                    const bdm::Double3 segment_vector = segment.vertex_1 - segment.vertex_0;
                    const double segment_length_squared = segment_vector * segment_vector;

                    ASSERT_(
                        std::isfinite(segment_length_squared) && segment_length_squared > 0.0,
                        "Initial attachment relocation encountered a zero-length segment"
                    );

                    const bdm::Double3 point_offset = original_position - segment.vertex_0;

                    double projected_fraction = (point_offset * segment_vector) / segment_length_squared;
                    projected_fraction = std::max(0.0,std::min(1.0,projected_fraction));
                    const bdm::Double3 closest_point = segment.vertex_0 + segment_vector * projected_fraction;
                    const double distance = L2norm(original_position - closest_point);

                    if (distance < nearest_distance) {

                        nearest_distance = distance;
                        nearest_segment_id = segment_id;
                        nearest_closest_point = closest_point;
                    }
                }
            };

        evaluate_segments(nearby_segment_ids);

        ASSERT_(
            nearest_segment_id > 0 && std::isfinite(nearest_distance),
            "Initial attachment relocation could not identify the nearest segment"
        );

        /*
        * The first broad-phase query can contain coarse bucket matches.
        * Re-query using the exact best distance so any genuinely closer
        * segment is also considered.
        */
        const double verification_radius = std::max(search_radius, nearest_distance);

        const std::vector<int> verification_segment_ids =
            active_scaffold.GetNearbySegmentIds(original_position, verification_radius);

        evaluate_segments(verification_segment_ids);

        // ---------------------------------------------------------------------
        // Place the cell immediately outside the nearest strut.
        // ---------------------------------------------------------------------

        const auto& nearest_segment = active_scaffold.GetSegment(nearest_segment_id);

        ASSERT_(
            std::isfinite(nearest_segment.radius) && nearest_segment.radius >= 0.0,
            "Initial attachment relocation encountered an invalid strut radius"
        );

        const double target_centreline_distance = cell_radius + nearest_segment.radius + scaffold_clearance;
        bdm::Double3 radial_direction = original_position - nearest_closest_point;
        double radial_norm = L2norm(radial_direction);

        // If the cell lies exactly on the segment centreline, choose one
        // deterministic perpendicular direction.
        if (radial_norm <= 1.0e-12) {

            const bdm::Double3 segment_vector = nearest_segment.vertex_1 - nearest_segment.vertex_0;

            const double segment_length = L2norm(segment_vector);

            ASSERT_(
                segment_length > 0.0,
                "Initial attachment relocation encountered a zero-length segment"
            );

            const bdm::Double3 axis = segment_vector * (1.0 / segment_length);

            const bdm::Double3 reference =
                std::fabs(axis[0]) < 0.9
                    ? bdm::Double3{1.0, 0.0, 0.0}
                    : bdm::Double3{0.0, 1.0, 0.0};

            radial_direction = {
                axis[1] * reference[2] - axis[2] * reference[1],
                axis[2] * reference[0] - axis[0] * reference[2],
                axis[0] * reference[1] - axis[1] * reference[0]
            };

            radial_norm = L2norm(radial_direction);

            ASSERT_(
                radial_norm > 0.0,
                "Initial attachment relocation could not construct "
                "a perpendicular direction"
            );
        }

        radial_direction *= 1.0 / radial_norm;

        const bdm::Double3 proposed_position = nearest_closest_point + radial_direction *
            target_centreline_distance;

        result.found = true;
        result.segment_id = nearest_segment_id;
        result.closest_point = nearest_closest_point;
        result.proposed_position = proposed_position;
        result.original_centreline_distance = nearest_distance;
        result.target_centreline_distance = target_centreline_distance;
        result.movement_distance = L2norm(proposed_position - original_position);

        return result;
    }

    inline
    InitialAttachmentScaffoldOverlapResult
    CheckInitialAttachmentScaffoldOverlap(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const bdm::Double3& proposed_position,
        const double scaffold_clearance) const {

        /*
        * Uses existing scaffold spatial index for local broad-phase query, then performs
        * exact point-to-segment distance checks to detect penetration of any nearby struts.
        */

        InitialAttachmentScaffoldOverlapResult result;

        ASSERT_(
            std::isfinite(scaffold_clearance) && scaffold_clearance >= 0.0,
            "Initial attachment scaffold-overlap check requires "
            "non-negative clearance"
        );

        const double cell_radius = 0.5 * cell.GetDiameter();

        ASSERT_(
            std::isfinite(cell_radius) && cell_radius > 0.0,
            "Initial attachment scaffold-overlap check requires "
            "a positive cell radius"
        );

        for (int coordinate = 0; coordinate < 3; ++coordinate) {

            ASSERT_(
                std::isfinite( proposed_position[coordinate] ),
                "Initial attachment scaffold-overlap check received "
                "a non-finite proposed position"
            );
        }

        // ---------------------------------------------------------------------
        // Retrieve only scaffold segments close enough to overlap.
        //
        // Segment spatial-index entries already include each strut radius,
        // therefore the query radius only needs the cell radius + clearance.
        // ---------------------------------------------------------------------

        const double search_radius = cell_radius + scaffold_clearance;

        const std::vector<int> nearby_segment_ids =
            active_scaffold.GetNearbySegmentIds(
                proposed_position,
                search_radius
            );

        double smallest_margin = std::numeric_limits<double>::max();

        // ---------------------------------------------------------------------
        // Calculate the exact distance to every nearby segment.
        // ---------------------------------------------------------------------

        for (const int segment_id : nearby_segment_ids) {

            ASSERT_(
                active_scaffold.HasSegment(segment_id),
                "Initial attachment scaffold-overlap check could not "
                "find scaffold segment"
            );

            const auto& segment = active_scaffold.GetSegment(segment_id);

            ASSERT_(
                std::isfinite(segment.radius) && segment.radius >= 0.0,
                "Initial attachment scaffold-overlap check encountered "
                "an invalid strut radius"
            );

            const bdm::Double3 segment_vector = segment.vertex_1 - segment.vertex_0;

            const double segment_length_squared = segment_vector * segment_vector;

            ASSERT_(
                std::isfinite(segment_length_squared) && segment_length_squared > 0.0,
                "Initial attachment scaffold-overlap check encountered "
                "a zero-length segment"
            );

            const bdm::Double3 point_offset = proposed_position - segment.vertex_0;

            double projected_fraction = (point_offset * segment_vector) / segment_length_squared;

            projected_fraction = std::max(0.0,std::min(1.0,projected_fraction));

            const bdm::Double3 closest_point = segment.vertex_0 + segment_vector * projected_fraction;
            const double centreline_distance = L2norm(proposed_position - closest_point);
            const double required_distance = cell_radius + segment.radius + scaffold_clearance;
            const double clearance_margin = centreline_distance - required_distance;

            if (clearance_margin < smallest_margin) {

                smallest_margin = clearance_margin;
                result.blocking_segment_id = segment_id;
                result.centreline_distance = centreline_distance;
                result.required_distance = required_distance;
                result.clearance_margin = clearance_margin;
            }
        }

        // No nearby scaffold means no possible overlap.
        if (nearby_segment_ids.empty()) {
            return result;
        }

        const double tolerance = 1.0e-12 * std::max(1.0,result.required_distance);
        result.overlaps = result.clearance_margin < -tolerance;

        return result;
    }

    inline
    void ImportInitialAttachmentConfiguration(
        bdm::Simulation& sim,
        const std::vector<ObstacleScaffold>& active_scaffolds) const{

        /*
        * Import the initial attachment configuration of all BiologicalCell agents.
        * This csv has been previously exported by ExportInitialAttachmentConfiguration().
        */

        ASSERT_(
            active_scaffolds.size() == 1,
            "Initial attachment import requires exactly one active scaffold"
        );

        const ObstacleScaffold& active_scaffold = active_scaffolds.front();

        ASSERT_(
            this->params()->have_parameter<std::string>(
                "initial_attachment_file"
            ),
            "Initial attachment import requires parameter "
            "\"initial_attachment_file\""
        );

        const std::string filename =
            this->params()->get<std::string>(
                "initial_attachment_file"
            );

        ASSERT_(
            !filename.empty(),
            "Initial attachment file path cannot be empty"
        );

        std::ifstream fin(filename);

        ASSERT_(
            fin.good(),
            "Could not open initial attachment file \""
            + filename + "\""
        );

        // -------------------------------------------------------------------------
        // Collect BiologicalCell pointers safely
        // -------------------------------------------------------------------------

        std::vector<bdm::BiologicalCell*>
            biological_cells;

        std::mutex biological_cells_mutex;

        auto* rm =
            sim.GetResourceManager();

        rm->ForEachAgent(
            [&](bdm::Agent* agent) {

                auto* cell =
                    dynamic_cast<bdm::BiologicalCell*>(
                        agent
                    );

                if (cell == nullptr) {
                    return;
                }

                std::lock_guard<std::mutex> lock(
                    biological_cells_mutex
                );

                biological_cells.push_back(
                    cell
                );
            }
        );

        // Deterministic order is useful for validation/debugging.
        std::sort(
            biological_cells.begin(),
            biological_cells.end(),
            [](const bdm::BiologicalCell* a,
            const bdm::BiologicalCell* b) {

                return
                    a->GetUid().GetIndex() <
                    b->GetUid().GetIndex();
            }
        );

        // -------------------------------------------------------------------------
        // Build UID lookup serially
        // -------------------------------------------------------------------------

        std::unordered_map<
            std::string,
            bdm::BiologicalCell*
        > cell_lookup;

        cell_lookup.reserve(
            biological_cells.size()
        );

        for (auto* cell :
            biological_cells) {

            std::ostringstream uid_stream;

            uid_stream <<
                cell->GetUid();

            const std::string abm_cell_id =
                uid_stream.str();

            ASSERT_(
                cell_lookup.find(
                    abm_cell_id
                ) ==
                cell_lookup.end(),
                "Initial attachment import found "
                "duplicate ABM cell UID: " +
                abm_cell_id
            );

            cell_lookup.emplace(
                abm_cell_id,
                cell
            );
        }

        ASSERT_(
            cell_lookup.size() ==
                biological_cells.size(),
            "Initial attachment import cell lookup "
            "has an inconsistent size"
        );



        // -------------------------------------------------------------------------
        // Read header
        // -------------------------------------------------------------------------

        std::string line;

        ASSERT_(
            static_cast<bool>(std::getline(fin, line)),
            "Initial attachment file is empty"
        );

        ASSERT_(
            line == "abm_cell_id,phenotype_id,x,y,z,attachment_node_ids",
            "Initial attachment file has an unexpected header"
        );

        // -------------------------------------------------------------------------
        // Read rows
        // -------------------------------------------------------------------------

        while (std::getline(fin, line)) {

            if (line.empty()) {
                continue;
            }

            std::stringstream ss(line);

            std::string abm_cell_id;
            std::string phenotype_text;
            std::string x_text;
            std::string y_text;
            std::string z_text;
            std::string attachment_text;

            std::getline(ss, abm_cell_id, ',');
            std::getline(ss, phenotype_text, ',');
            std::getline(ss, x_text, ',');
            std::getline(ss, y_text, ',');
            std::getline(ss, z_text, ',');

            // The attachment list is the final field and contains commas itself.
            std::getline(ss, attachment_text);

            ASSERT_(
                !abm_cell_id.empty() &&
                !phenotype_text.empty() &&
                !x_text.empty() &&
                !y_text.empty() &&
                !z_text.empty() &&
                !attachment_text.empty(),
                "Initial attachment import encountered an incomplete CSV row"
            );

            auto cell_it =
                cell_lookup.find(
                    abm_cell_id
                );

            ASSERT_(
                cell_it != cell_lookup.end(),
                "Initial attachment import could not find ABM cell "
                + abm_cell_id
            );

            bdm::BiologicalCell* cell =
                cell_it->second;

            const int phenotype_id =
                std::stoi(
                    phenotype_text
                );

            ASSERT_(
                cell->GetPhenotype() == phenotype_id,
                "Initial attachment import found a phenotype mismatch for cell "
                + abm_cell_id
            );

            const double x =
                std::stod(
                    x_text
                );

            const double y =
                std::stod(
                    y_text
                );

            const double z =
                std::stod(
                    z_text
                );

            // Remove surrounding quotes.
            if (!attachment_text.empty() &&
                attachment_text.front() == '"') {

                attachment_text.erase(
                    attachment_text.begin()
                );
            }

            if (!attachment_text.empty() &&
                attachment_text.back() == '"') {

                attachment_text.pop_back();
            }

            // Remove [ ].
            ASSERT_(
                attachment_text.size() >= 2 &&
                attachment_text.front() == '[' &&
                attachment_text.back() == ']',
                "Initial attachment import found invalid attachment-node format"
            );

            attachment_text =
                attachment_text.substr(
                    1,
                    attachment_text.size() - 2
                );

            std::vector<int> attachment_node_ids;

            std::stringstream attachment_stream(
                attachment_text
            );

            std::string node_text;

            while (std::getline(
                
                attachment_stream,
                node_text,
                ',')) {

                attachment_node_ids.push_back(
                    std::stoi(
                        node_text
                    )
                );
            }

            ASSERT_(
                !attachment_node_ids.empty(),
                "Initial attachment import found no attachment node IDs for cell "
                + abm_cell_id
            );

            // -------------------------------------------------------------------------
            // Warn about attachment geometry outside the configured separation range
            // -------------------------------------------------------------------------

            if (attachment_node_ids.size() >= 2) {

                const std::string& phenotype_name =
                    this->params()->get<std::string>(
                        "phenotype_ID/" +
                        std::to_string(
                            cell->GetPhenotype()
                        )
                    );

                const std::string mech_base =
                    phenotype_name +
                    "/cell_matrix_mechanics";

                const double min_attachment_separation =
                    this->params()->get<double>(
                        mech_base +
                        "/min_attachment_separation"
                    );

                const double max_attachment_separation =
                    this->params()->get<double>(
                        mech_base +
                        "/max_attachment_separation"
                    );

                double imported_min_separation = std::numeric_limits<double>::max();
                double imported_max_separation = 0.0;

                for (std::size_t i = 0; i < attachment_node_ids.size(); ++i) {

                    ASSERT_(
                        active_scaffold.HasNode(
                            attachment_node_ids[i]
                        ),
                        "Initial attachment import could not find scaffold node ID "
                        + std::to_string(
                            attachment_node_ids[i]
                        )
                    );

                    const bdm::Double3& position_i =
                        active_scaffold.GetNodePosition(
                            attachment_node_ids[i]
                        );

                    for (std::size_t j = i + 1; j < attachment_node_ids.size(); ++j) {

                        ASSERT_(
                            active_scaffold.HasNode(
                                attachment_node_ids[j]
                            ),
                            "Initial attachment import could not find scaffold node ID "
                            + std::to_string(
                                attachment_node_ids[j]
                            )
                        );

                        const bdm::Double3& position_j =
                            active_scaffold.GetNodePosition(
                                attachment_node_ids[j]
                            );

                        const bdm::Double3 difference =
                            position_i - position_j;

                        const double separation =
                            std::sqrt(
                                difference[0] * difference[0] +
                                difference[1] * difference[1] +
                                difference[2] * difference[2]
                            );

                        imported_min_separation =
                            std::min(
                                imported_min_separation,
                                separation
                            );

                        imported_max_separation =
                            std::max(
                                imported_max_separation,
                                separation
                            );
                    }
                }

                if (imported_min_separation <
                        min_attachment_separation ||
                    imported_max_separation >
                        max_attachment_separation) {

                    std::cerr
                        << "[INITIAL ATTACHMENT WARNING] Cell "
                        << abm_cell_id
                        << " contains imported attachment separations outside "
                        << "the configured range ["
                        << min_attachment_separation
                        << ", "
                        << max_attachment_separation
                        << "]. Observed range = ["
                        << imported_min_separation
                        << ", "
                        << imported_max_separation
                        << "]. Configuration will still be imported."
                        << std::endl;
                }
            }

            cell->SetPosition(
                bdm::Double3{x,y,z}
            );

            cell->InitialiseAttachmentRecordsFromNodeIds(
                active_scaffold,
                attachment_node_ids
            );
        }
    }

    inline
    void ExportCellPositions(
        bdm::Simulation& sim,
        const std::map<int, std::string>& cells,
        int time){
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
        * time,abm_cell_id,x,y,z,cell_state,attachment_node_ids,
        * contractile_force,pairwise_contractile_forces
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

        fpos
            << "time,abm_cell_id,x,y,z,cell_state,"
            << "attachment_node_ids,contractile_force,"
            << "pair_contractile_forces\n";

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

        const int debug_cells_to_print = 3;
        int debug_exported_cells = 0;

        for (auto* cell : biological_cells) {

            const auto& pos = cell->GetPosition();

            std::ostringstream uid_stream;
            uid_stream << cell->GetUid();

            const std::string abm_cell_id = uid_stream.str();

            // Check if extra printing is enabled for this cell type. This is useful for debugging
            bool mechanics_verbose = false;
            const int phenotype_id = cell->GetPhenotype();
            const auto phenotype_it = cells.find(phenotype_id);

            if (phenotype_it != cells.end() && phenotype_id >= 1) {

                const std::string& CP_name = phenotype_it->second;

                const std::string mech_base = CP_name + "/cell_matrix_mechanics";

                mechanics_verbose =
                    this->params()->have_parameter<bool>(
                        mech_base + "/verbose"
                    ) &&
                    this->params()->get<bool>(
                        mech_base + "/verbose"
                    );
            }

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
            std::string pair_contractile_forces_text = "[]";

            double contractile_force = 0.0;

            if (export_action == FemExportAction::kContract) {

                cell_state = "contract";

                attachment_node_ids_text =
                    this->FormatAttachmentNodeIdsJsonLike(
                        cell->GetAttachmentRecords()
                    );

                contractile_force = cell->GetContractileForce();

                // Pair forces are valid only when they correspond to the current
                // attachment set.
                const std::size_t n_attachments = cell->GetNumberOfAttachmentRecords();

                const std::size_t expected_pairs = n_attachments * (n_attachments - 1) / 2;

                ASSERT_(
                    cell->GetNumberOfPairMechanicsRecords() == expected_pairs,
                    "ABM-to-FEM export: pair-mechanics count does not match "
                    "the current attachment set for cell " + abm_cell_id
                );

                pair_contractile_forces_text =
                    this->FormatPairContractileForcesJsonLike(
                        cell->GetPairMechanicsRecords()
                    );
                

                // A negative contractile force is not expected by the model.
                // Report the issue and prevent it from reaching the FEM input.
                if (contractile_force < 0.0) {
                    std::cerr
                    << "[CELL-MATRIX WARNING] Negative contractile force for cell "
                    << abm_cell_id 
                    << ": " << contractile_force 
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
                pair_contractile_forces_text = "[]";
            }

            // -------------------------------------------------------------------------
            // Step 5c: Record the FEM communication state on the cell
            // -------------------------------------------------------------------------

            cell->SetCellState(cell_state);

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
                << contractile_force <<  ","
                << "\"" << pair_contractile_forces_text << "\""
                << "\n";

            // Optional debug output for the first few exported cells.
            if (mechanics_verbose && debug_exported_cells < debug_cells_to_print) {

                std::cout
                    << "[FEM EXPORT] abm_cell_id=" << abm_cell_id
                    << " pos=(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")"
                    << " state=" << cell_state
                    << " k_ce=" << cell->GetKce()
                    << " contractile_force=" << contractile_force
                    << " n_attach=" << cell->GetNumberOfAttachmentRecords()
                    << "\n";

                for (const auto& pair_record :
                    cell->GetPairMechanicsRecords()) {

                    std::cout
                        << "    [PAIR EXPORT] pair=("
                        << pair_record.node_id_a << ", "
                        << pair_record.node_id_b << ")"
                        << " k_ce=" << pair_record.k_ce
                        << " k_ce_valid=" << pair_record.has_valid_k_ce
                        << " force=" << pair_record.contractile_force
                        << "\n";
                }
            }

            csv_lines.push_back(line.str());

            ++debug_exported_cells;

        }

        // ---------------------------------------------------------------------------
        // Step 5: Write CSV lines serially and close the file
        // ---------------------------------------------------------------------------

        for (const auto& line : csv_lines) {
            fpos << line;
        }

        fpos.close();

    }

  inline
  int RunFemSolver(
      bdm::Simulation& sim,
      const std::map<int, std::string>& cells,
      int time,
      bool* ran_any){
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
                --min_cell_radius <min_cell_radius <- legacy> \
                --max_cell_radius <max_cell_radius <- legacy> \
                --strut_radius <strut_radius> \
                --delta_F <delta_F> \
                --perturbance_dist <perturbance_dist> \
                --random_state <random_state> \
                --num_attachments <num_attachments> \
                --lattice_mesh_path "<lattice_mesh_path>" \
                --verbose \ <= IF TRUE 
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

            const double min_attachment_separation =
                this->params()->get<double>(
                    mech_base + "/min_attachment_separation"
                );

            if ((2.0 * min_attachment_separation) <
                this->params()->get<double>(CP_name + "/diameter/min")) {
                ABORT_(
                    "Model parameter '"
                    + mech_base
                    + "/min_attachment_separation' cannot be less than '"
                    + CP_name
                    + "/diameter/min'"
                );
                }

            const double max_attachment_separation =
                this->params()->get<double>(
                    mech_base + "/max_attachment_separation"
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

            const double min_contractile_force =
                this->params()->get<double>(
                    mech_base + "/min_contractile_force"
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
            + std::to_string(min_attachment_separation/2)
            + " ";

            cmd +=
            "--max_cell_radius "
            + std::to_string(max_attachment_separation/2)
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
            "--min_contractile_force "
            + std::to_string(min_contractile_force)
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

    inline
    void WriteScaffoldVersion2(
        const ObstacleScaffold& scaffold,
        const std::string& output_path) const{
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
  
    inline
    void PrepareActiveScaffoldSpatialIndex(
        const std::map<int, std::string>& cells,
        std::vector<ObstacleScaffold>& active_scaffolds) const {
        
        /*
        * Build the reusable spatial node and segment indices for the current
        * active scaffold.
        *
        * The index is rebuilt once after each scaffold replacement and is reused
        * by all cell candidate-radius searches during that timestep.
        */

        // -------------------------------------------------------------------------
        // Step 1: Exit when cell-matrix mechanics is disabled
        // -------------------------------------------------------------------------

        if (!this->AnyCellMatrixMechanicsEnabled(cells)) {
        return;
        }

        ASSERT_(
            active_scaffolds.size() == 1,
            "Scaffold spatial-index preparation requires exactly one active "
            "scaffold when cell-matrix mechanics is enabled"
        );

        // -------------------------------------------------------------------------
        // Step 2: Find the mechanics-enabled phenotype
        // -------------------------------------------------------------------------

        bool found_mechanics_phenotype = false;
        double max_attachment_separation = 0.0;

        for (const auto& cell_type : cells) {
        const int phenotype_id =
            cell_type.first;

        if (phenotype_id < 1) {
            continue;
        }

        const std::string& phenotype_name =
            cell_type.second;

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
            continue;
        }

        ASSERT_(
            !found_mechanics_phenotype,
            "Scaffold spatial-index preparation currently supports only one "
            "mechanics-enabled phenotype"
        );

        max_attachment_separation =
            this->params()->get<double>(
                mech_base + "/max_attachment_separation"
            );

        found_mechanics_phenotype = true;
        }

        ASSERT_(
            found_mechanics_phenotype,
            "Scaffold spatial-index preparation could not find a "
            "mechanics-enabled phenotype"
        );

        ASSERT_(
            std::isfinite(max_attachment_separation) &&
            max_attachment_separation > 0.0,
            "Scaffold spatial-index preparation requires a positive maximum "
            "attachment separation"
        );

        // -------------------------------------------------------------------------
        // Step 3: Build the index using the maximum pairwise search distance
        // -------------------------------------------------------------------------

        const double bucket_size =
            max_attachment_separation;

        active_scaffolds.front().BuildNodeSpatialIndex(
            bucket_size
        );

        active_scaffolds.front().BuildSegmentSpatialIndex(
            bucket_size
        );

    }    

    inline
    void ResetInitialAttachmentPathCache(
        const double bucket_size) {
        
            ASSERT_(
            std::isfinite(bucket_size) && bucket_size > 0.0,
            "Initial attachment path cache requires a positive bucket size"
        );

        // Reset initial attachment buckets
        initial_attachment_path_records_.clear();
        initial_attachment_path_spatial_index_.clear();

        initial_attachment_path_bucket_size_ = bucket_size;

        initial_attachment_path_max_length_ = 0.0;
    }

    

    inline
    void AddInitialAttachmentPathToCache(
        const unsigned int cell_uid,
        const bdm::Double3& path_start,
        const bdm::Double3& path_end) {

        /*
        * Add the initial attachment paths to cache for faster lookup/
        */

        ASSERT_(
            std::isfinite(initial_attachment_path_bucket_size_) &&
            initial_attachment_path_bucket_size_ > 0.0,
            "Initial attachment path cache has not been initialised"
        );

        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {

            ASSERT_(
                std::isfinite(path_start[coordinate]) && std::isfinite(path_end[coordinate]),
                "Initial attachment path cache received non-finite coordinates"
            );
        }

        const double path_length = L2norm(path_end - path_start);

        ASSERT_(
            std::isfinite(path_length) && path_length > 0.0,
            "Initial attachment path cache received a zero-length path"
        );

        initial_attachment_path_max_length_ =
            std::max(
                initial_attachment_path_max_length_,
                path_length
            );

        InitialAttachmentPathRecord record;

        record.cell_uid = cell_uid;
        record.start = path_start;
        record.end = path_end;
        record.midpoint = (path_start + path_end) * 0.5;
        record.length = path_length;

        const std::size_t record_index = initial_attachment_path_records_.size();

        initial_attachment_path_records_.push_back(record);

        const InitialAttachmentPathBucketKey bucket_key = {
            static_cast<int>(
                std::floor(
                    record.midpoint[0] /
                    initial_attachment_path_bucket_size_
                )
            ),
            static_cast<int>(
                std::floor(
                    record.midpoint[1] /
                    initial_attachment_path_bucket_size_
                )
            ),
            static_cast<int>(
                std::floor(
                    record.midpoint[2] /
                    initial_attachment_path_bucket_size_
                )
            )
        };

        initial_attachment_path_spatial_index_[bucket_key].push_back(record_index);
    }


    inline
    std::vector<std::size_t>
    GetNearbyInitialAttachmentPathIndices(
        const bdm::Double3& query_centre,
        const double search_radius) const {

        /*
        * Find neaby attachment paths to determine if cells are obstructing each other.
        */

        ASSERT_(
            std::isfinite(initial_attachment_path_bucket_size_) &&
            initial_attachment_path_bucket_size_ > 0.0,
            "Initial attachment path cache has not been initialised"
        );

        ASSERT_(
            std::isfinite(search_radius) && search_radius >= 0.0,
            "Initial attachment path query requires a non-negative radius"
        );

        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {

            ASSERT_(
                std::isfinite(query_centre[coordinate]),
                "Initial attachment path query received a non-finite centre"
            );
        }

        const int min_x = static_cast<int>(
            std::floor(
                (query_centre[0] - search_radius) /
                initial_attachment_path_bucket_size_
            )
        );

        const int max_x = static_cast<int>(
            std::floor(
                (query_centre[0] + search_radius) /
                initial_attachment_path_bucket_size_
            )
        );

        const int min_y = static_cast<int>(
            std::floor(
                (query_centre[1] - search_radius) /
                initial_attachment_path_bucket_size_
            )
        );

        const int max_y = static_cast<int>(
            std::floor(
                (query_centre[1] + search_radius) /
                initial_attachment_path_bucket_size_
            )
        );

        const int min_z = static_cast<int>(
            std::floor(
                (query_centre[2] - search_radius) /
                initial_attachment_path_bucket_size_
            )
        );

        const int max_z = static_cast<int>(
            std::floor(
                (query_centre[2] + search_radius) /
                initial_attachment_path_bucket_size_
            )
        );

        std::vector<std::size_t> nearby_path_indices;

        for (int x = min_x; x <= max_x; ++x) {
            for (int y = min_y; y <= max_y; ++y) {
                for (int z = min_z; z <= max_z; ++z) {

                    const InitialAttachmentPathBucketKey key = {
                        x, y, z
                    };

                    const auto bucket_it =
                        initial_attachment_path_spatial_index_.find(
                            key
                        );

                    if (bucket_it == initial_attachment_path_spatial_index_.end()) {
                        continue;
                    }

                    for (const std::size_t record_index : bucket_it->second) {

                        ASSERT_(
                            record_index < initial_attachment_path_records_.size(),
                            "Initial attachment path cache contains an invalid index"
                        );

                        const auto& record = initial_attachment_path_records_[
                            record_index
                        ];

                        const double midpoint_distance = L2norm(record.midpoint - query_centre);

                        if (midpoint_distance <= search_radius) {

                            nearby_path_indices.push_back(
                                record_index
                            );
                        }
                    }
                }
            }
        }

        return nearby_path_indices;
    }


    inline
    const InitialAttachmentPathRecord&
    GetInitialAttachmentPathRecord(
        const std::size_t record_index) const
    {
        ASSERT_(
            record_index <
                initial_attachment_path_records_.size(),
            "Requested initial attachment path index is out of range"
        );

        return initial_attachment_path_records_[
            record_index
        ];
    }

    inline
    bool InitialAttachmentPathObstructed(
        const bdm::BiologicalCell& candidate_cell,
        const bdm::Double3& candidate_start,
        const bdm::Double3& candidate_end,
        const double minimum_path_separation) const {

        /*
        * Scalable spatial query followed by segment-to-segment distance checks
        * only for nearby cell paths.
        */

        ASSERT_(
            std::isfinite(minimum_path_separation) && minimum_path_separation >= 0.0,
            "Initial attachment cell-path separation must be finite and non-negative"
        );

        if (initial_attachment_path_records_.empty()) {
            return false;
        }

        ASSERT_(
            std::isfinite(initial_attachment_path_max_length_) &&
            initial_attachment_path_max_length_ > 0.0,
            "Initial attachment path cache contains paths but has an invalid "
            "maximum path length"
        );

        const double candidate_length = L2norm(candidate_end - candidate_start);

        ASSERT_(
            std::isfinite(candidate_length) && candidate_length > 0.0,
            "Initial attachment cell-path obstruction requires a positive "
            "candidate path length"
        );

        const bdm::Double3 candidate_midpoint = (candidate_start + candidate_end) * 0.5;

        /*
        * If two finite segments can come within minimum_path_separation,
        * their midpoints cannot be farther apart than half the length of
        * each segment plus that separation.
        */
        const double broad_phase_radius = 
            0.5 * candidate_length + 0.5 * initial_attachment_path_max_length_ + minimum_path_separation;

        const std::vector<std::size_t> nearby_path_indices =
            this->GetNearbyInitialAttachmentPathIndices(
                candidate_midpoint,
                broad_phase_radius
            );

        for (const std::size_t path_index : nearby_path_indices) {

            const InitialAttachmentPathRecord& existing_path =
                this->GetInitialAttachmentPathRecord(
                    path_index
                );

            const double path_separation =
                candidate_cell.CalculateSegmentDistance(
                    candidate_start,
                    candidate_end,
                    existing_path.start,
                    existing_path.end
                );

            ASSERT_(
                std::isfinite(path_separation) &&
                path_separation >= 0.0,
                "Initial attachment cell-path obstruction calculated an "
                "invalid path separation"
            );

            const double tolerance = 1.0e-12 * std::max(1.0, minimum_path_separation);

            if (path_separation <= minimum_path_separation + tolerance) {
                return true;
            }
        }

        return false;
    }

    inline
    void FilterInitialAttachmentPairsByCellPathObstruction(
        const bdm::BiologicalCell& candidate_cell,
        const ObstacleScaffold& active_scaffold,
        const std::vector<std::pair<int, int>>& candidate_pairs,
        const std::vector<double>& candidate_weights,
        const double minimum_path_separation,
        std::vector<std::pair<int, int>>& valid_pairs,
        std::vector<double>& valid_weights) const {

        /*
        * Impossible pairs are removed together with their corresponding weights 
        * before selection probabilities are calculated.
        */

        ASSERT_(
            candidate_pairs.size() == candidate_weights.size(),
            "Initial attachment pair and weight counts do not match "
            "before cell-path obstruction filtering"
        );

        ASSERT_(
            std::isfinite(minimum_path_separation) && minimum_path_separation >= 0.0,
            "Initial attachment cell-path separation must be "
            "finite and non-negative"
        );

        valid_pairs.clear();
        valid_weights.clear();

        valid_pairs.reserve(candidate_pairs.size());
        valid_weights.reserve(candidate_weights.size());

        for (std::size_t i = 0; i < candidate_pairs.size(); ++i) {

            const auto& pair = candidate_pairs[i];

            ASSERT_(
                pair.first > 0 && pair.second > 0 && pair.first != pair.second,
                "Cell-path obstruction filtering received an invalid "
                "initial attachment pair"
            );

            ASSERT_(
                active_scaffold.HasNode(pair.first) &&
                active_scaffold.HasNode(pair.second),
                "Cell-path obstruction filtering could not find an "
                "attachment node"
            );

            ASSERT_(
                std::isfinite(candidate_weights[i]) && candidate_weights[i] > 0.0,
                "Cell-path obstruction filtering received an invalid weight"
            );

            const bdm::Double3& path_start =
                active_scaffold.GetNodePosition(
                    pair.first
                );

            const bdm::Double3& path_end =
                active_scaffold.GetNodePosition(
                    pair.second
                );

            const bool obstructed =
                this->InitialAttachmentPathObstructed(
                    candidate_cell,
                    path_start,
                    path_end,
                    minimum_path_separation
                );

            if (obstructed) {
                continue;
            }

            valid_pairs.push_back(pair);

            valid_weights.push_back(candidate_weights[i]);
        }

        ASSERT_(
            valid_pairs.size() == valid_weights.size(),
            "Cell-path obstruction filtering produced misaligned "
            "pairs and weights"
        );
    }

    inline
    InitialAttachmentGenerationResult
    GenerateInitialAttachmentForCell(
        bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const double minimum_path_separation,
        const unsigned int draw_index = 0) {

        /*
        * Initial cell attachment at step = 0, only for cells 
        * mechanically interacting with scaffold obstacle.
        */

        InitialAttachmentGenerationResult result;

        using LifecycleStatus = bdm::BiologicalCell::CellMatrixLifecycleStatus;

        ASSERT_(
            cell.GetCellMatrixLifecycleStatus() == LifecycleStatus::kNeedsInitialAttachment ||
            cell.GetCellMatrixLifecycleStatus() == LifecycleStatus::kInitialAttachmentRetry,
            "Initial attachment generation requires an "
            "uninitialised or retrying cell"
        );

        ASSERT_(
            !cell.HasAttachmentRecords(),
            "Initial attachment generation requires no "
            "existing attachment records"
        );

        ASSERT_(
            std::isfinite(minimum_path_separation) && minimum_path_separation >= 0.0,
            "Initial attachment minimum cell-path separation "
            "must be finite and non-negative"
        );

        ASSERT_(
            std::isfinite(initial_attachment_path_bucket_size_) &&
            initial_attachment_path_bucket_size_ > 0.0,
            "Initial attachment path cache must be initialised "
            "before generating attachments"
        );

        // Find nearby scaffold nodes
        const std::vector<int> candidate_node_ids =
            cell.GenerateInitialAttachmentCandidateNodeIds(
                active_scaffold
            );

        result.candidate_node_count = candidate_node_ids.size();

        if (candidate_node_ids.size() < 2) {
            return result;
        }

        // Generate separation-valid unordered pairs
        const std::vector<std::pair<int, int>>
            candidate_pairs =
                cell.GenerateInitialAttachmentCandidatePairs(
                    active_scaffold,
                    candidate_node_ids
                );

        result.candidate_pair_count = candidate_pairs.size();

        if (candidate_pairs.empty()) {
            return result;
        }

        // Calculate midpoint-distance weights
        const std::string& phenotype_name =
            this->params()->get<std::string>(
                "phenotype_ID/" +
                std::to_string(
                    cell.GetPhenotype()
                )
            );

        const std::string mech_base =
            phenotype_name + "/cell_matrix_mechanics";

        const double distance_sensitivity =
            this->params()->get<double>(
                mech_base + "/candidate_cell_distance_sensitivity"
            );

        const std::vector<double> candidate_weights =
            cell.CalculateInitialAttachmentPairWeights(
                active_scaffold,
                candidate_pairs,
                distance_sensitivity
            );

        ASSERT_(
            candidate_weights.size() ==
                candidate_pairs.size(),
            "Initial attachment pair and weight counts do not match"
        );

        // Remove scaffold-obstructed pairs
        std::vector<std::pair<int, int>>
            scaffold_valid_pairs;

        std::vector<double>
            scaffold_valid_weights;

        cell.FilterInitialAttachmentPairsByScaffoldObstruction(
            active_scaffold,
            candidate_pairs,
            candidate_weights,
            scaffold_valid_pairs,
            scaffold_valid_weights
        );

        result.scaffold_valid_pair_count =
            scaffold_valid_pairs.size();

        if (scaffold_valid_pairs.empty()) {
            return result;
        }

        // Remove pairs obstructed by already accepted cell paths
        std::vector<std::pair<int, int>> cell_valid_pairs;

        std::vector<double> cell_valid_weights;

        this->FilterInitialAttachmentPairsByCellPathObstruction(
            cell,
            active_scaffold,
            scaffold_valid_pairs,
            scaffold_valid_weights,
            minimum_path_separation,
            cell_valid_pairs,
            cell_valid_weights
        );

        result.cell_valid_pair_count =
            cell_valid_pairs.size();

        if (cell_valid_pairs.empty()) {
            return result;
        }

        // Convert surviving weights to probabilities
        const std::vector<double> selection_probabilities =
            cell.CalculateInitialAttachmentPairSelectionProbabilities(
                cell_valid_weights
            );

        ASSERT_(
            selection_probabilities.size() ==
                cell_valid_pairs.size(),
            "Initial attachment pair and probability counts do not match"
        );

        // Select one pair reproducibly
        const std::pair<int, int> selected_pair =
            cell.SelectInitialAttachmentPair(
                cell_valid_pairs,
                selection_probabilities,
                draw_index
            );

        // Commit the selected attachment pair
        cell.ApplySelectedInitialAttachmentPair(
            active_scaffold,
            selected_pair
        );

        // Add the accepted cell path to the spatial cache
          this->AddInitialAttachmentPathToCache(
            static_cast<unsigned int>(
                cell.GetUid().GetIndex()
            ),
            active_scaffold.GetNodePosition(
                selected_pair.first
            ),
            active_scaffold.GetNodePosition(
                selected_pair.second
            )
        );

        result.success =
            true;

        result.selected_pair =
            selected_pair;

        return result;
    }
    
    inline
    InitialAttachmentWorkflowResult
    GenerateInitialAttachmentWithRelocation(
        bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold,
        const double minimum_path_separation,
        const double relocation_search_radius,
        const std::size_t maximum_segment_attempts,
        const std::size_t positions_per_segment,
        const double scaffold_clearance,
        const double cell_clearance,
        const unsigned int draw_index = 0)
    {
        InitialAttachmentWorkflowResult result;

        ASSERT_(
            !cell.HasAttachmentRecords(),
            "Initial attachment workflow requires a cell "
            "without existing attachment records"
        );

        ASSERT_(
            std::isfinite(minimum_path_separation) &&
            minimum_path_separation >= 0.0,
            "Initial attachment workflow requires "
            "non-negative path separation"
        );

        ASSERT_(
            std::isfinite(relocation_search_radius) &&
            relocation_search_radius > 0.0,
            "Initial attachment workflow requires "
            "a positive relocation search radius"
        );

        ASSERT_(
            maximum_segment_attempts > 0,
            "Initial attachment workflow requires "
            "at least one relocation segment attempt"
        );

        ASSERT_(
            positions_per_segment > 0,
            "Initial attachment workflow requires "
            "at least one relocation position per segment"
        );

        const bdm::Double3 original_position =
            cell.GetPosition();

        const double cell_radius =
            0.5 *
            cell.GetDiameter();

        ASSERT_(
            std::isfinite(cell_radius) &&
            cell_radius > 0.0,
            "Initial attachment workflow requires "
            "a positive cell radius"
        );

        // ---------------------------------------------------------------------
        // Step 1: Check whether the original ABM position is already occupied
        // ---------------------------------------------------------------------

        const auto original_cell_overlap =
            this->CheckInitialAttachmentCellOverlap(
                static_cast<unsigned int>(
                    cell.GetUid().GetIndex()
                ),
                original_position,
                cell_radius,
                cell_clearance
            );

        result.original_position_blocked_by_cell =
            original_cell_overlap.overlaps;

        // ---------------------------------------------------------------------
        // Step 2: Try attachment at the original ABM position only when free
        // ---------------------------------------------------------------------

        if (!result.original_position_blocked_by_cell) {

            result.original_generation_result =
                this->GenerateInitialAttachmentForCell(
                    cell,
                    active_scaffold,
                    minimum_path_separation,
                    draw_index
                );

            if (result.original_generation_result.success) {

                result.success =
                    true;

                result.attached_at_original_position =
                    true;

                result.final_generation_result =
                    result.original_generation_result;

                this->AddInitialAttachmentCellToCache(
                    static_cast<unsigned int>(
                        cell.GetUid().GetIndex()
                    ),
                    cell.GetPosition(),
                    cell_radius
                );

                return result;
            }
        }

        ASSERT_(
            !cell.HasAttachmentRecords(),
            "Failed initial attachment attempt unexpectedly "
            "modified the cell attachment state"
        );

        // ---------------------------------------------------------------------
        // Step 2: Generate nearest-first relocation struts
        // ---------------------------------------------------------------------

        result.relocation_attempted =
            true;

        const std::vector<int> ordered_segment_ids =
            this->GenerateInitialAttachmentRelocationSegmentIds(
                cell,
                active_scaffold,
                relocation_search_radius,
                maximum_segment_attempts
            );

        // ---------------------------------------------------------------------
        // Step 3: Try each strut from nearest to furthest
        // ---------------------------------------------------------------------

        for (std::size_t segment_index = 0;
            segment_index < ordered_segment_ids.size();
            ++segment_index) {

            ++result.attempted_segment_count;

            const int segment_id =
                ordered_segment_ids[
                    segment_index
                ];

            const auto relocation =
                this->GenerateInitialAttachmentRelocationCandidateForSegment(
                    cell,
                    active_scaffold,
                    segment_id,
                    scaffold_clearance
                );

            ASSERT_(
                relocation.found,
                "Initial attachment workflow could not "
                "generate relocation geometry"
            );

            const auto candidate_positions =
                this->GenerateInitialAttachmentRelocationPositions(
                    active_scaffold,
                    relocation,
                    positions_per_segment
                );

            // -----------------------------------------------------------------
            // Step 4: Try every position around this strut
            // -----------------------------------------------------------------

            for (std::size_t position_index = 0;
                position_index < candidate_positions.size();
                ++position_index) {

                ++result.attempted_position_count;

                const bdm::Double3& candidate_position =
                    candidate_positions[
                        position_index
                    ];

                // -------------------------------------------------------------
                // Reject scaffold-overlapping positions
                // -------------------------------------------------------------

                const auto scaffold_overlap =
                    this->CheckInitialAttachmentScaffoldOverlap(
                        cell,
                        active_scaffold,
                        candidate_position,
                        scaffold_clearance
                    );

                if (scaffold_overlap.overlaps) {
                    continue;
                }

                // -------------------------------------------------------------
                // Reject positions occupied by already accepted cells
                // -------------------------------------------------------------

                const auto cell_overlap =
                    this->CheckInitialAttachmentCellOverlap(
                        static_cast<unsigned int>(
                            cell.GetUid().GetIndex()
                        ),
                        candidate_position,
                        cell_radius,
                        cell_clearance
                    );

                if (cell_overlap.overlaps) {
                    continue;
                }

                ++result.valid_relocation_position_count;

                // -------------------------------------------------------------
                // Step 5: Temporarily move the cell and retry attachment
                // -------------------------------------------------------------

                cell.SetPosition(
                    candidate_position
                );

                const unsigned int retry_draw_index =
                    draw_index +
                    static_cast<unsigned int>(
                        result.valid_relocation_position_count
                    );

                const auto attachment_result =
                    this->GenerateInitialAttachmentForCell(
                        cell,
                        active_scaffold,
                        minimum_path_separation,
                        retry_draw_index
                    );

                if (attachment_result.success) {

                    result.success =
                        true;

                    result.selected_segment_id =
                        segment_id;

                    result.selected_position_index =
                        position_index;

                    result.selected_relocation_position =
                        candidate_position;

                    result.final_generation_result =
                        attachment_result;

                    // GenerateInitialAttachmentForCell may update the final cell
                    // position after selecting the attachment pair.
                    this->AddInitialAttachmentCellToCache(
                        static_cast<unsigned int>(
                            cell.GetUid().GetIndex()
                        ),
                        cell.GetPosition(),
                        cell_radius
                    );

                    return result;
                }

                // -------------------------------------------------------------
                // Failed attachment attempt must leave no persistent state.
                // Restore the ABM position before trying another relocation.
                // -------------------------------------------------------------

                ASSERT_(
                    !cell.HasAttachmentRecords(),
                    "Failed relocated attachment attempt unexpectedly "
                    "modified persistent attachment state"
                );

                cell.SetPosition(
                    original_position
                );
            }
        }

        // ---------------------------------------------------------------------
        // Step 6: Complete failure
        //
        // No valid attachment could be formed at the original position or after
        // exhausting the permitted relocation attempts.
        // ---------------------------------------------------------------------

        cell.SetPosition(
            original_position
        );

        ASSERT_(
            !cell.HasAttachmentRecords(),
            "Initial attachment failure left persistent attachment records"
        );

        ASSERT_(
            !cell.RequiresMechanicsRecalculation(),
            "Initial attachment failure left pending mechanics recalculation"
        );

        // ---------------------------------------------------------------------
        // Step 7: Mark the cell for removal
        // ---------------------------------------------------------------------

        cell.MarkCellMatrixForDeath();
        cell.ValidateCellMatrixState();

        result.marked_for_death = true;

        result.removal_requested = true;

        // ---------------------------------------------------------------------
        // Step 8: Remove the failed cell from BioDynaMo
        //
        // Do not access 'cell' after this call.
        // ---------------------------------------------------------------------

        cell.RemoveFromSimulation();

        return result;
    }
    
    inline
    InitialAttachmentPopulationResult
    GenerateInitialAttachmentsForPopulation(
        bdm::Simulation& sim,
        const std::map<int, std::string>& cells,
        const std::vector<ObstacleScaffold>& active_scaffolds)
    {
        InitialAttachmentPopulationResult result;

        if (!this->AnyCellMatrixMechanicsEnabled(cells)) {
            return result;
        }

        ASSERT_(
            active_scaffolds.size() == 1,
            "Initial attachment population generation requires "
            "exactly one active scaffold"
        );

        const ObstacleScaffold& active_scaffold =
            active_scaffolds.front();

        // ---------------------------------------------------------------------
        // Step 1: Find the single mechanics-enabled phenotype
        // ---------------------------------------------------------------------

        int mechanics_phenotype_id = -1;

        std::string mechanics_base;

        for (const auto& cell_type : cells) {

            const int phenotype_id =
                cell_type.first;

            if (phenotype_id < 1) {
                continue;
            }

            const std::string candidate_mechanics_base =
                cell_type.second +
                "/cell_matrix_mechanics";

            const bool mechanics_enabled =
                this->params()->have_parameter<bool>(
                    candidate_mechanics_base + "/enabled"
                ) &&
                this->params()->get<bool>(
                    candidate_mechanics_base + "/enabled"
                );

            if (!mechanics_enabled) {
                continue;
            }

            ASSERT_(
                mechanics_phenotype_id < 0,
                "Initial attachment population generation currently "
                "supports only one mechanics-enabled phenotype"
            );

            mechanics_phenotype_id =
                phenotype_id;

            mechanics_base =
                candidate_mechanics_base;
        }

        ASSERT_(
            mechanics_phenotype_id > 0,
            "Initial attachment population generation could not find "
            "a mechanics-enabled phenotype"
        );

        const int maximum_segment_attempts_input =
            this->params()->get<int>(
                mechanics_base +
                "/initial_attachment_max_strut_candidates"
            );

        const int positions_per_segment_input =
            this->params()->get<int>(
                mechanics_base +
                "/initial_attachment_positions_per_strut"
            );

        ASSERT_(
            maximum_segment_attempts_input > 0,
            "Initial attachment population generation requires "
            "at least one permitted strut"
        );

        ASSERT_(
            positions_per_segment_input > 0,
            "Initial attachment population generation requires "
            "at least one relocation position per strut"
        );

        const std::size_t maximum_segment_attempts =
            static_cast<std::size_t>(
                maximum_segment_attempts_input
            );

        const std::size_t positions_per_segment =
            static_cast<std::size_t>(
                positions_per_segment_input
            );

        // ---------------------------------------------------------------------
        // Step 2: Read shared initialisation parameters
        // ---------------------------------------------------------------------

        const double max_attachment_separation =
            this->params()->get<double>(
                mechanics_base +
                "/max_attachment_separation"
            );

        const double scaffold_clearance =
            this->params()->get<double>(
                mechanics_base +
                "/candidate_scaffold_clearance"
            );

        const double cell_clearance =
            this->params()->get<double>(
                mechanics_base +
                "/candidate_cell_clearance"
            );

        /*
        * Current initial path-clearance definition.
        *
        * We can separate this into a dedicated input parameter later if required.
        */
        const double minimum_path_separation =
            cell_clearance;

        ASSERT_(
            std::isfinite(max_attachment_separation) &&
            max_attachment_separation > 0.0,
            "Initial attachment population generation requires "
            "positive maximum attachment separation"
        );

        // ---------------------------------------------------------------------
        // Step 3: Collect mechanics-enabled cells
        // ---------------------------------------------------------------------

        std::vector<bdm::BiologicalCell*>
            initial_cells;

        std::mutex initial_cells_mutex;

        auto* rm =
            sim.GetResourceManager();

        ASSERT_(
            rm != nullptr,
            "Initial attachment population generation requires "
            "a BioDynaMo resource manager"
        );

        rm->ForEachAgent(
            [&](bdm::Agent* agent) {

                auto* cell =
                    dynamic_cast<bdm::BiologicalCell*>(
                        agent
                    );

                if (cell == nullptr ||
                    cell->GetPhenotype() !=
                        mechanics_phenotype_id) {

                    return;
                }

                std::lock_guard<std::mutex> lock(
                    initial_cells_mutex
                );

                initial_cells.push_back(
                    cell
                );
            }
        );

        // ---------------------------------------------------------------------
        // Step 4: Enforce deterministic population order
        // ---------------------------------------------------------------------

        std::sort(
            initial_cells.begin(),
            initial_cells.end(),
            [](const bdm::BiologicalCell* a,
            const bdm::BiologicalCell* b) {

                return
                    a->GetUid().GetIndex() <
                    b->GetUid().GetIndex();
            }
        );

        // ---------------------------------------------------------------------
        // Step 5: Reset the shared initialisation caches ONCE
        // ---------------------------------------------------------------------

        this->ResetInitialAttachmentPathCache(
            max_attachment_separation
        );

        this->ResetInitialAttachmentCellCache(
            max_attachment_separation
        );

        using LifecycleStatus =
            bdm::BiologicalCell::
                CellMatrixLifecycleStatus;

        using LifecycleStatus =
            bdm::BiologicalCell::CellMatrixLifecycleStatus;

        bool requires_initial_generation =
            false;

        // -------------------------------------------------------------------------
        // Seed caches from cells that already have established attachments
        // -------------------------------------------------------------------------

        for (auto* cell : initial_cells) {

            const auto lifecycle_status =
                cell->GetCellMatrixLifecycleStatus();

            if (lifecycle_status ==
                LifecycleStatus::kEstablished) {

                ASSERT_(
                    cell->HasAttachmentRecords(),
                    "Established initialisation cell "
                    "has no attachment records"
                );

                this->SeedInitialAttachmentCellCacheFromEstablishedCell(
                    *cell
                );

                this->SeedInitialAttachmentPathCacheFromEstablishedCell(
                    *cell,
                    active_scaffold
                );

                ++result.seeded_established_cell_count;

                continue;
            }

            if (lifecycle_status ==
                    LifecycleStatus::kNeedsInitialAttachment ||
                lifecycle_status ==
                    LifecycleStatus::kInitialAttachmentRetry) {

                requires_initial_generation =
                    true;

                continue;
            }

            if (lifecycle_status ==
                LifecycleStatus::kMarkedForDeath) {

                continue;
            }

            ASSERT_(
                false,
                "Initial attachment population encountered "
                "an unexpected lifecycle state"
            );
        }

        // If every mechanics-enabled cell already has attachments,
        // the caches are now reconstructed and there is nothing to generate.
        if (!requires_initial_generation) {
            return result;
        }

        // ---------------------------------------------------------------------
        // Step 6: Process cells sequentially in UID order
        // ---------------------------------------------------------------------

        bool first_cell_recorded =
            false;

        for (auto* cell : initial_cells) {

            const auto lifecycle_status =
                cell->GetCellMatrixLifecycleStatus();

            if (lifecycle_status ==
                    LifecycleStatus::kEstablished ||
                lifecycle_status ==
                    LifecycleStatus::kMarkedForDeath) {

                continue;
            }

            ASSERT_(
                cell != nullptr,
                "Initial attachment population contains "
                "a null BiologicalCell pointer"
            );

            /*
            * This method is currently the generation route.
            * Imported/reused configurations will be handled separately.
            */
            ASSERT_(
                lifecycle_status ==
                    LifecycleStatus::kNeedsInitialAttachment ||
                lifecycle_status ==
                    LifecycleStatus::kInitialAttachmentRetry,
                "Initial attachment generation route encountered "
                "a mechanics-enabled cell that was already established"
            );

            cell->ValidateCellMatrixState();

            const unsigned int cell_uid =
                static_cast<unsigned int>(
                    cell->GetUid().GetIndex()
                );

            if (!first_cell_recorded) {

                result.first_processed_uid =
                    cell_uid;

                first_cell_recorded =
                    true;
            }

            result.last_processed_uid =
                cell_uid;

            ++result.processed_cell_count;

            // -------------------------------------------------------------
            // Run the complete validated per-cell workflow
            // -------------------------------------------------------------

            const InitialAttachmentWorkflowResult
                workflow_result =
                    this->GenerateInitialAttachmentWithRelocation(
                        *cell,
                        active_scaffold,
                        minimum_path_separation,
                        max_attachment_separation,
                        maximum_segment_attempts,
                        positions_per_segment,
                        scaffold_clearance,
                        cell_clearance,
                        0
                    );

            result.total_segment_attempt_count +=
                workflow_result.attempted_segment_count;

            result.total_position_attempt_count +=
                workflow_result.attempted_position_count;

            if (workflow_result
                    .original_position_blocked_by_cell) {

                ++result.original_position_blocked_count;
            }

            // -------------------------------------------------------------
            // Successful cell
            // -------------------------------------------------------------

            if (workflow_result.success) {

                const auto additional_attachment_result =
                    this->CompleteInitialAdditionalAttachmentsForCell(
                        *cell,
                        active_scaffold
                    );

                ++result.attached_cell_count;

                if (workflow_result
                        .attached_at_original_position) {

                    ++result.attached_at_original_count;

                } else {

                    ++result
                        .attached_after_relocation_count;
                }

                ASSERT_(
                    cell->GetCellMatrixLifecycleStatus() ==
                        LifecycleStatus::kEstablished,
                    "Successful initial attachment did not establish the cell"
                );

                const std::string& phenotype_name =
                    this->params()->get<std::string>(
                        "phenotype_ID/" +
                        std::to_string(
                            cell->GetPhenotype()
                        )
                    );

                const int maximum_attachment_count =
                    this->params()->get<int>(
                        phenotype_name +
                        "/cell_matrix_mechanics/num_attachments"
                    );

                ASSERT_(
                    maximum_attachment_count >= 2,
                    "Initial attachment population generation "
                    "requires num_attachments >= 2"
                );

                ASSERT_(
                    cell->GetNumberOfAttachmentRecords() >= 2,
                    "Successful initial attachment generation "
                    "must retain at least two attachments"
                );

                ASSERT_(
                    cell->GetNumberOfAttachmentRecords() <=
                        static_cast<std::size_t>(
                            maximum_attachment_count
                        ),
                    "Initial attachment generation exceeded "
                    "num_attachments"
                );

                ASSERT_(
                    additional_attachment_result.final_attachment_count ==
                        cell->GetNumberOfAttachmentRecords(),
                    "Additional initial attachment result "
                    "does not match the final cell state"
                );

                continue;
            }

            // -------------------------------------------------------------
            // Complete failure
            //
            // The cell has already been handed to BioDynaMo for removal.
            // Do not access 'cell' below this point.
            // -------------------------------------------------------------

            ASSERT_(
                workflow_result.marked_for_death &&
                workflow_result.removal_requested,
                "Failed initial attachment did not follow "
                "the expected cell-removal route"
            );

            ++result.removed_cell_count;
        }

        // ---------------------------------------------------------------------
        // Step 7: Validate population accounting
        // ---------------------------------------------------------------------

        ASSERT_(
            result.attached_cell_count +
                result.removed_cell_count ==
            result.processed_cell_count,
            "Initial attachment population accounting is inconsistent"
        );

        ASSERT_(
            result.attached_at_original_count +
                result.attached_after_relocation_count ==
            result.attached_cell_count,
            "Initial attachment success accounting is inconsistent"
        );

        return result;
    }
    
    inline
    std::size_t GetInitialAttachmentCellCount() const
    {
        return initial_attachment_cell_records_.size();
    }
    
    inline
    std::size_t GetInitialAttachmentPathCount() const
    {
        return initial_attachment_path_records_.size();
    }
    
    inline
    void SeedInitialAttachmentCellCacheFromEstablishedCell(
        const bdm::BiologicalCell& cell)
    {
        using LifecycleStatus =
            bdm::BiologicalCell::CellMatrixLifecycleStatus;

        ASSERT_(
            cell.GetCellMatrixLifecycleStatus() ==
                LifecycleStatus::kEstablished,
            "Initial cell-cache seeding requires "
            "an established cell"
        );

        ASSERT_(
            cell.HasAttachmentRecords(),
            "Established cell used for initial cell-cache "
            "seeding has no attachment records"
        );

        const double cell_radius =
            0.5 *
            cell.GetDiameter();

        ASSERT_(
            std::isfinite(cell_radius) &&
            cell_radius > 0.0,
            "Initial cell-cache seeding encountered "
            "an invalid cell radius"
        );

        this->AddInitialAttachmentCellToCache(
            static_cast<unsigned int>(
                cell.GetUid().GetIndex()
            ),
            cell.GetPosition(),
            cell_radius
        );
    }
    
    inline
    std::size_t SeedInitialAttachmentPathCacheFromEstablishedCell(
        const bdm::BiologicalCell& cell,
        const ObstacleScaffold& active_scaffold)
    {
        using LifecycleStatus =
            bdm::BiologicalCell::CellMatrixLifecycleStatus;

        ASSERT_(
            cell.GetCellMatrixLifecycleStatus() ==
                LifecycleStatus::kEstablished,
            "Initial attachment path-cache seeding "
            "requires an established cell"
        );

        ASSERT_(
            cell.HasAttachmentRecords(),
            "Established cell used for initial path-cache "
            "seeding has no attachment records"
        );

        const auto& attachments =
            cell.GetAttachmentRecords();

        std::size_t added_path_count =
            0;

        // Add every unordered attachment pair.
        for (std::size_t i = 0;
            i < attachments.size();
            ++i) {

            ASSERT_(
                active_scaffold.HasNode(
                    attachments[i].node_id
                ),
                "Initial path-cache seeding could not find "
                "an established attachment node"
            );

            for (std::size_t j = i + 1;
                j < attachments.size();
                ++j) {

                ASSERT_(
                    active_scaffold.HasNode(
                        attachments[j].node_id
                    ),
                    "Initial path-cache seeding could not find "
                    "an established attachment node"
                );

                this->AddInitialAttachmentPathToCache(
                    static_cast<unsigned int>(
                        cell.GetUid().GetIndex()
                    ),
                    active_scaffold.GetNodePosition(
                        attachments[i].node_id
                    ),
                    active_scaffold.GetNodePosition(
                        attachments[j].node_id
                    )
                );

                ++added_path_count;
            }
        }

        const std::size_t expected_path_count =
            (
                attachments.size() *
                (
                    attachments.size() - 1
                )
            ) /
            2;

        ASSERT_(
            added_path_count ==
                expected_path_count,
            "Established-cell path-cache seeding "
            "added an incorrect number of paths"
        );

        return added_path_count;
    }
    
    inline
    void SynchroniseAttachmentCoordinates(
        bdm::Simulation& sim,
        const std::map<int, std::string>& cells,
        const std::vector<ObstacleScaffold>& active_scaffolds) const
    {
        /*
        * Function goal
        * -------------
        * Update every retained attachment coordinate using its persistent scaffold
        * node ID and the current active scaffold geometry.
        *
        * Only the attachment position is refreshed. The node ID, k_ecm validity,
        * newly formed state and cell-level mechanics state remain unchanged.
        */

        // -------------------------------------------------------------------------
        // Step 1: Exit when cell-matrix mechanics is not enabled
        // -------------------------------------------------------------------------

        if (!this->AnyCellMatrixMechanicsEnabled(cells)) {
            return;
        }

        // -------------------------------------------------------------------------
        // Step 2: Validate the active coupled scaffold
        // -------------------------------------------------------------------------

        ASSERT_(
            this->params()->get<int>("simulation_obstacles") > 0,
            "Cell-matrix mechanics requires at least one scaffold obstacle"
        );

        ASSERT_(
            this->params()->get<bool>("simulation_obstacles/update"),
            "Cell-matrix mechanics requires simulation_obstacles/update to be true "
            "so attachment coordinates use the current scaffold geometry"
        );
        
        ASSERT_(
            active_scaffolds.size() == 1,
            "Attachment synchronisation requires exactly one active scaffold when "
            "cell-matrix mechanics is enabled"
        );

        const ObstacleScaffold& active_scaffold =
            active_scaffolds.front();

        ASSERT_(
            !active_scaffold.nodes_by_id.empty(),
            "Attachment synchronisation received an active scaffold with no nodes"
        );

        // -------------------------------------------------------------------------
        // Step 3: Collect BiologicalCell agents
        // -------------------------------------------------------------------------

        std::vector<bdm::BiologicalCell*> biological_cells;
        std::mutex biological_cells_mutex;

        auto* rm = sim.GetResourceManager();

        rm->ForEachAgent([&](bdm::Agent* agent) {
        auto* cell =
            dynamic_cast<bdm::BiologicalCell*>(agent);

        if (!cell) {
            return;
        }

        std::lock_guard<std::mutex> lock(
            biological_cells_mutex
        );

        biological_cells.push_back(cell);
        });

        // -------------------------------------------------------------------------
        // Step 4: Refresh each retained attachment coordinate
        // -------------------------------------------------------------------------

        for (auto* cell : biological_cells) {
            if (!cell->HasAttachmentRecords()) {
                continue;
            }

            std::ostringstream uid_stream;
            uid_stream << cell->GetUid();

            const std::string abm_cell_id =
                uid_stream.str();

            // Preserve the cell-level mechanics state so coordinate-only
            // synchronisation can be validated.
            const auto lifecycle_status_before =
                cell->GetCellMatrixLifecycleStatus();

            const bool requires_recalculation_before =
                cell->RequiresMechanicsRecalculation();

            const double k_ce_before =
                cell->GetKce();

            const double contractile_force_before =
                cell->GetContractileForce();

            // Work on a copy so only the attachment coordinates are changed.
            const auto original_attachment_records =
                cell->GetAttachmentRecords();

            const auto attachment_records =
                this->RefreshAttachmentRecordCoordinates(
                    original_attachment_records,
                    active_scaffold,
                    abm_cell_id
                );

            // SetAttachmentRecords validates persistent IDs and duplicate nodes.
            cell->SetAttachmentRecords(attachment_records);

            // -----------------------------------------------------------------------
            // Step 5: Confirm that only attachment coordinates changed
            // -----------------------------------------------------------------------

            ASSERT_(
                cell->GetNumberOfAttachmentRecords()
                    == original_attachment_records.size(),
                "Attachment synchronisation changed the attachment count for cell "
                + abm_cell_id
            );

            for (std::size_t i = 0;
                i < original_attachment_records.size();
                ++i) {

                const auto& before =
                    original_attachment_records[i];

                const auto& after =
                    cell->GetAttachmentRecord(i);

                ASSERT_(
                    after.node_id == before.node_id,
                    "Attachment synchronisation changed the persistent node ID for cell "
                    + abm_cell_id
                );

                ASSERT_(
                    after.newly_formed == before.newly_formed,
                    "Attachment synchronisation changed newly_formed for cell "
                    + abm_cell_id
                );

                ASSERT_(
                    after.attachment_age == before.attachment_age,
                    "Attachment synchronisation changed attachment age for cell "
                    + abm_cell_id
                );

                ASSERT_(
                    after.detachment_threshold == before.detachment_threshold,
                    "Attachment synchronisation changed detachment threshold for cell "
                    + abm_cell_id
                );


            }

            ASSERT_(
                cell->GetCellMatrixLifecycleStatus() == lifecycle_status_before,
                "Attachment synchronisation changed the lifecycle status for cell "
                + abm_cell_id
            );

            ASSERT_(
                cell->RequiresMechanicsRecalculation() == requires_recalculation_before,
                "Attachment synchronisation changed the mechanics-recalculation state "
                "for cell " + abm_cell_id
            );

            ASSERT_(
                cell->GetKce() == k_ce_before,
                "Attachment synchronisation changed k_ce for cell "
                + abm_cell_id
            );

            ASSERT_(
                cell->GetContractileForce() == contractile_force_before,
                "Attachment synchronisation changed the contractile force for cell "
                + abm_cell_id
            );

            // Coordinate refresh must not create an inconsistent lifecycle state.
            cell->ValidateCellMatrixState();

            // -----------------------------------------------------------------------
            // Step 6: Reposition established single-attachment cells
            // -----------------------------------------------------------------------

            using LifecycleStatus =
                bdm::BiologicalCell::CellMatrixLifecycleStatus;

            if (lifecycle_status_before == LifecycleStatus::kEstablished &&
                original_attachment_records.size() == 1) {

                const int attachment_node_id =
                    original_attachment_records.front().node_id;

                const bdm::Double3 previous_attachment_position =
                    original_attachment_records.front().position;

                const bdm::Double3 cell_position_before = cell->GetPosition();

                cell->FollowSingleAttachmentScaffold(
                    active_scaffold,
                    previous_attachment_position
                );

            }
        }
    }
    
    // ---------------------------------------------------------------------------
    // FEM-to-ABM import
    // ---------------------------------------------------------------------------

    inline
    void ImportFemCells(
        bdm::Simulation& sim,
        const std::map<int, std::string>& cells,
        int time) {

        /*
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
        // Build FEM result filename
        // ---------------------------------------------------------------------------

        char buf[2048];

        const int time_inc = time - 1;

        std::snprintf(buf, sizeof(buf),
                        "./results/FEM/step_%d/cell_mechanics_step_%d.dat",
                        time_inc, time_inc);

        const std::string fname(buf);

        // ---------------------------------------------------------------------------
        // Open FEM result file
        // ---------------------------------------------------------------------------

        std::ifstream fin(fname);

        ASSERT_(fin, "FEM import: could not open file " + fname);

        // ---------------------------------------------------------------------------
        // Read number of cells
        // ---------------------------------------------------------------------------

        int n_cells = -1;

        fin >> n_cells;

        ASSERT_(fin && n_cells >= 0,
                "FEM import: invalid n_cells in " + fname);

        // ---------------------------------------------------------------------------
        // Build ABM cell lookup using BioDynaMo UID string
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
        // Validate the expected FEM result count
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
        // Track only the cells expected from the FEM
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
            // Read attachment node IDs
            // -------------------------------------------------------------------------

            std::vector<int> attachment_node_ids(n_attach);

            for (int a = 0; a < n_attach; ++a) {
            fin >> attachment_node_ids[a];

            ASSERT_(fin,
                    "FEM import: not enough attachment node IDs for abm_cell_id = "
                    + abm_cell_id + " in " + fname);
            }

            // -------------------------------------------------------------------------
            // Read attachment coordinates
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
            // Read pairwise k_ce records
            // -------------------------------------------------------------------------

            int n_pairs = -1;

            fin >> n_pairs;

            ASSERT_(
                fin && n_pairs >= 0,
                "FEM import: invalid n_pairs for abm_cell_id = " + abm_cell_id + " in " + fname
            );

            const int expected_n_pairs = n_attach * (n_attach - 1) / 2;

            ASSERT_(
                n_pairs == expected_n_pairs,
                "FEM import: incorrect number of pairwise k_ce records for abm_cell_id = " + abm_cell_id
                + ". Expected " + std::to_string(expected_n_pairs)
                + ", received " + std::to_string(n_pairs)
            );

            std::vector<int> pair_node_a(n_pairs);
            std::vector<int> pair_node_b(n_pairs);
            std::vector<double> pair_kce(n_pairs);

            // Read the pairwise k_ce records
            for (int p = 0; p < n_pairs; ++p) {

                fin
                    >> pair_node_a[p]
                    >> pair_node_b[p]
                    >> pair_kce[p];

                ASSERT_(
                    fin,
                    "FEM import: failed reading pairwise k_ce record " + std::to_string(p)
                    + " for abm_cell_id = " + abm_cell_id
                );
            }

            // Validate the node pairs
            for (int p = 0; p < n_pairs; ++p) {

                const int node_a = pair_node_a[p];
                const int node_b = pair_node_b[p];

                const bool node_a_exists =
                    std::find(
                        attachment_node_ids.begin(),
                        attachment_node_ids.end(),
                        node_a
                    ) != attachment_node_ids.end();

                const bool node_b_exists =
                    std::find(
                        attachment_node_ids.begin(),
                        attachment_node_ids.end(),
                        node_b
                    ) != attachment_node_ids.end();

                ASSERT_(
                    node_a_exists && node_b_exists,
                    "FEM import: pairwise k_ce contains a node that is not part of "
                    "the attachment set for abm_cell_id = " + abm_cell_id
                );

                ASSERT_(
                    node_a != node_b,
                    "FEM import: pairwise k_ce contains a self-pair for abm_cell_id = " + abm_cell_id
                );

                ASSERT_(
                    std::isfinite(pair_kce[p]),
                    "FEM import: pairwise k_ce is non-finite for abm_cell_id = " + abm_cell_id
                );

                ASSERT_(
                    pair_kce[p] >= 0.0,
                    "FEM import: pairwise k_ce is negative for abm_cell_id = " + abm_cell_id
                );
            }

            // Check for duplicate pairs (order-independent)
            for (int p = 0; p < n_pairs; ++p) {

                const int a1 = std::min(pair_node_a[p], pair_node_b[p]);
                const int b1 = std::max(pair_node_a[p], pair_node_b[p]);

                for (int q = p + 1; q < n_pairs; ++q) {

                    const int a2 = std::min(pair_node_a[q], pair_node_b[q]);
                    const int b2 = std::max(pair_node_a[q], pair_node_b[q]);

                    ASSERT_(
                        !(a1 == a2 && b1 == b2),
                        "FEM import: duplicate pairwise k_ce record for pair ("
                        + std::to_string(a1) + ", " + std::to_string(b1)
                        + ") in abm_cell_id = " + abm_cell_id
                    );
                }
            }

            // Convert temporary FEM data into persistent pairwise k_ce records for the cell.
            std::vector<bdm::BiologicalCell::PairMechanicsRecord> pair_mechanics_records;
            pair_mechanics_records.reserve(n_pairs);

            for (int p = 0; p < n_pairs; ++p) {

                bdm::BiologicalCell::PairMechanicsRecord record;

                // Store node IDs canonically so A < B.
                record.node_id_a = std::min(pair_node_a[p], pair_node_b[p]);
                record.node_id_b = std::max(pair_node_a[p], pair_node_b[p]);
                record.k_ce = pair_kce[p];
                record.has_valid_k_ce = true;

                pair_mechanics_records.push_back(record);
            }

            // -------------------------------------------------------------------------
            // Build persistent attachment records
            // -------------------------------------------------------------------------

            std::vector<bdm::BiologicalCell::AttachmentRecord> attachment_records;

            const auto previous_attachment_records = cell->GetAttachmentRecords();

            attachment_records.reserve(n_attach);

            for (int a = 0; a < n_attach; ++a) {

                bdm::BiologicalCell::AttachmentRecord record;

                record.node_id = attachment_node_ids[a];
                record.position = attachment_points[a];

                // New FEM-returned attachments begin with zero age.
                // Their detachment threshold is assigned below if no previous
                // persistent attachment with the same node ID exists.
                record.attachment_age = 0.0;
                record.detachment_threshold = -1.0;

                bool retained_attachment = false;

                for (const auto& previous_record :
                    previous_attachment_records) {

                    if (previous_record.node_id ==
                        record.node_id) {

                        record.attachment_age =  previous_record.attachment_age;
                        record.detachment_threshold = previous_record.detachment_threshold;
                        retained_attachment = true;

                        break;
                    }
                }

                // A genuinely new FEM attachment receives its threshold exactly once.
                if (!retained_attachment) {

                    record.detachment_threshold =
                        cell->GenerateAttachmentDetachmentThreshold();
                }

                record.newly_formed = false;
                attachment_records.push_back(record);
            }

            // -------------------------------------------------------------------------
            // Update the cell and store the imported attachment state
            //
            // The legacy vectors remain populated for backward compatibility while the
            // persistent attachment records are introduced gradually.
            // -------------------------------------------------------------------------

            cell->SetPosition(bdm::Double3{x, y, z});

            // Store the persistent records and synchronise the compatibility vectors.
            cell->SetAttachmentRecords(attachment_records);

            // Store the pairwise mechanics belonging to this attachment set.
            cell->SetPairMechanicsRecords(pair_mechanics_records);

            // -------------------------------------------------------------------------
            // Validate the stored persistent attachment state
            // -------------------------------------------------------------------------

            ASSERT_(
                cell->GetNumberOfAttachmentRecords()
                    == attachment_records.size(),
                "FEM import: stored attachment-record count does not match the imported "
                "count for abm_cell_id = " + abm_cell_id
            );

            for (std::size_t a = 0;
                a < cell->GetNumberOfAttachmentRecords();
                ++a) {

                const auto& record =
                    cell->GetAttachmentRecord(a);

                ASSERT_(
                    record.node_id == attachment_node_ids[a],
                    "FEM import: stored attachment node ID does not match the imported "
                    "node ID for abm_cell_id = " + abm_cell_id
                );

                ASSERT_(
                    record.position[0] == attachment_points[a][0] &&
                    record.position[1] == attachment_points[a][1] &&
                    record.position[2] == attachment_points[a][2],
                    "FEM import: stored attachment position does not match the imported "
                    "position for abm_cell_id = " + abm_cell_id
                );

                ASSERT_(
                    std::isfinite(record.detachment_threshold) &&
                    record.detachment_threshold >= 0.0 &&
                    record.detachment_threshold <= 1.0,
                    "FEM import: stored attachment has an invalid detachment threshold "
                    "for abm_cell_id = " + abm_cell_id
                );

                ASSERT_(
                    !record.newly_formed,
                    "FEM import: FEM-returned attachment was incorrectly marked as newly "
                    "formed for abm_cell_id = " + abm_cell_id
                );
            }
            
            // Validate the stored persistent pairwise mechanics state
            ASSERT_(
                cell->GetNumberOfPairMechanicsRecords() == pair_mechanics_records.size(),
                "FEM import: stored pair-mechanics count does not match the imported "
                "count for abm_cell_id = " + abm_cell_id
            );

            for (std::size_t p = 0;
                p < cell->GetNumberOfPairMechanicsRecords();
                ++p) {

                const auto& stored_record = cell->GetPairMechanicsRecord(p);

                const auto& imported_record = pair_mechanics_records[p];

                ASSERT_(
                    stored_record.node_id_a == imported_record.node_id_a && stored_record.node_id_b == imported_record.node_id_b,
                    "FEM import: stored pair node IDs do not match the imported pair "
                    "for abm_cell_id = " + abm_cell_id
                );

                ASSERT_(
                    stored_record.k_ce == imported_record.k_ce,
                    "FEM import: stored pair k_ce does not match the imported value "
                    "for abm_cell_id = " + abm_cell_id
                );
            }

            cell->SetKce(k_ce);

            // -------------------------------------------------------------------------
            // Update adaptive cell-specific contractile force
            // -------------------------------------------------------------------------

            double adaptive_contractile_force = 0.0;

            const int phenotype_id = cell->GetPhenotype();

            auto phenotype_it = cells.find(phenotype_id);

            bool mechanics_verbose = false;

            if (phenotype_it != cells.end() && phenotype_id >= 1) {

                const std::string& CP_name = phenotype_it->second;

                const std::string mech_base = CP_name + "/cell_matrix_mechanics";

                mechanics_verbose =
                    this->params()->have_parameter<bool>(
                        mech_base + "/verbose"
                    ) &&
                    this->params()->get<bool>(
                        mech_base + "/verbose"
                    );
            }

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

                    // Calculate one contractile force for each attachment pair.
                    for (std::size_t p = 0;
                        p < cell->GetNumberOfPairMechanicsRecords();
                        ++p) {

                        const auto& pair_record = cell->GetPairMechanicsRecord(p);

                        const double pair_contractile_force =
                            calculate_contractile_force_from_kce(
                                pair_record.k_ce,
                                min_contractile_force,
                                max_contractile_force,
                                contractility_kce_sensitivity,
                                contractility_kce_midpoint
                            );

                        cell->SetPairContractileForce(
                            p,
                            pair_contractile_force
                        );
                    }
                        
                }
            }

            cell->SetContractileForce(adaptive_contractile_force);

            // -------------------------------------------------------------------------
            // Update the internal ABM mechanics state
            // -------------------------------------------------------------------------

            if (n_attach > 0) {
            // The cell now has an established attachment set managed by the ABM.
            cell->SetCellMatrixLifecycleStatus(
                bdm::BiologicalCell::CellMatrixLifecycleStatus::kEstablished
            );

            // The imported k_ce, contractile force and attachment mechanics correspond
            // to the attachment set returned by the FEM.
            cell->ClearMechanicsRecalculationRequirement();
            }

            // Confirm that the imported mechanics and lifecycle state agree.
            cell->ValidateCellMatrixState();

            // Do not update cell_state here.
            // cell_state_ represents the last ABM-to-FEM state exported by ABM.
            imported_flag_it->second = true;

            if (mechanics_verbose && i < debug_cells_to_print) {
                std::cout << "[FEM IMPORT] abm_cell_id=" << abm_cell_id
                    << " pos=(" << x << ", " << y << ", " << z << ")"
                    << " k_ce=" << k_ce
                    << " contractile_force=" << adaptive_contractile_force
                    << " n_attach=" << n_attach
                    << "\n";
                
                // Print the pairwise k_ce records for debugging
                for (const auto& pair_record : cell->GetPairMechanicsRecords()) {

                    std::cout
                        << "    [PAIR KCE] pair=(" << pair_record.node_id_a << ", " << pair_record.node_id_b << ")"
                        << " k_ce=" << pair_record.k_ce
                        << " k_ce_valid=" << pair_record.has_valid_k_ce
                        << " force=" << pair_record.contractile_force
                        << "\n";
                }
            }

        }

        // ---------------------------------------------------------------------------
        // Confirm every ABM cell was imported exactly once
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
    
    inline
    std::vector<bdm::BiologicalCell::AttachmentRecord>
    RefreshAttachmentRecordCoordinates(
        const std::vector<bdm::BiologicalCell::AttachmentRecord>& records,
        const ObstacleScaffold& active_scaffold,
        const std::string& abm_cell_id) const{
        
        /*
        * Function goal
        * -------------
        * Return a copy of the supplied attachment records with each coordinate
        * refreshed from the active scaffold using its persistent node ID.
        *
        * All non-coordinate attachment fields remain unchanged.
        */
        
        auto refreshed_records = records;

        for (auto& record : refreshed_records) {
        ASSERT_(
            record.node_id > 0,
            "Attachment synchronisation encountered a non-positive node ID "
            "for cell " + abm_cell_id
        );

        ASSERT_(
            active_scaffold.HasNode(record.node_id),
            "Attachment synchronisation could not find scaffold node ID "
            + std::to_string(record.node_id)
            + " for cell "
            + abm_cell_id
        );

        ASSERT_(
            !active_scaffold
                .GetConnectedNodeIds(record.node_id)
                .empty(),
            "Attachment synchronisation found an unconnected scaffold node ID "
            + std::to_string(record.node_id)
            + " for cell "
            + abm_cell_id
        );

        record.position =
            active_scaffold.GetNodePosition(record.node_id);
        }

        return refreshed_records;
    }

    double initial_attachment_path_bucket_size_ = 0.0;
    double initial_attachment_path_max_length_ = 0.0;

    std::vector<InitialAttachmentPathRecord> initial_attachment_path_records_;

    std::unordered_map<
        InitialAttachmentPathBucketKey,
        std::vector<std::size_t>,
        InitialAttachmentPathBucketKeyHash>
        initial_attachment_path_spatial_index_;

    double initial_attachment_cell_bucket_size_ = 0.0;
    double initial_attachment_cell_max_radius_ = 0.0;

    std::vector<InitialAttachmentCellRecord>
        initial_attachment_cell_records_;

    std::unordered_map<
        InitialAttachmentCellBucketKey,
        std::vector<std::size_t>,
        InitialAttachmentCellBucketKeyHash
    >
initial_attachment_cell_spatial_index_;
    
    // Model parameters are owned by ABM4bio, not by this class.
    Parameters* params_;

    std::unordered_set<std::string> exported_cell_ids_;
};

// =============================================================================
// Method definitions
// =============================================================================

#endif  // CELL_MATRIX_INTERACTION_H_