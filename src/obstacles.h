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
#ifndef _OBSTACLES_H_
#define _OBSTACLES_H_
// =============================================================================
#include "./global.h"
#include "core/container/math_array.h"
// =============================================================================
class Obstacle {
public:
  Obstacle() {}
  ~Obstacle() {}
//
public:
  unsigned int id;
  std::string type;
};
// -----------------------------------------------------------------------------
class ObstacleScaffold : public Obstacle {
 public:
  
  // ---------------------------------------------------------------------------
  // Spatial node-index bucket
  // ---------------------------------------------------------------------------

  struct SpatialBucketKey {
    int x;
    int y;
    int z;

    bool operator==(const SpatialBucketKey& other) const {
      return x == other.x &&
              y == other.y &&
              z == other.z;
    }
  };

  struct SpatialBucketKeyHash {
    std::size_t operator()(
        const SpatialBucketKey& key) const {
      /*
        * Combine the three integer bucket coordinates into one hash value.
        */

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

  // ---------------------------------------------------------------------------
  // Persistent scaffold node
  // ---------------------------------------------------------------------------

  struct ScaffoldNode {
    int id;
    bdm::Double3 position;
    double radius;
    std::vector<int> connected_node_ids;
  };

  // ---------------------------------------------------------------------------
  // Geometric segment representation
  //
  // This keeps the existing interface used by biological_cell-inline.h while
  // also retaining the persistent node and element identifiers.
  // ---------------------------------------------------------------------------

  struct Segment {
    int element_id;
    int node_id_1;
    int node_id_2;

    bdm::Double3 vertex_0;
    bdm::Double3 vertex_1;

    double length;
    double radius;
  };

 public:
  ObstacleScaffold() {}

  ~ObstacleScaffold() {}

  // ---------------------------------------------------------------------------
  // Read scaffold file
  //
  // Supported formats
  // -----------------
  //
  // Version 1.0:
  //
  //   1.0
  //   <number_of_nodes>
  //   <x> <y> <z> <radius>
  //   ...
  //   <number_of_elements>
  //   <zero_based_node_index_1> <zero_based_node_index_2>
  //
  // Version 2.0:
  //
  //   2.0
  //   <number_of_nodes>
  //   <node_id> <x> <y> <z> <radius>
  //   ...
  //   <number_of_elements>
  //   <element_id> <node_id_1> <node_id_2>
  // ---------------------------------------------------------------------------

  inline void init(const std::string& t, const std::string& fname) {
    this->type = t;

    // -------------------------------------------------------------------------
    // Step 1: Open scaffold file
    // -------------------------------------------------------------------------

    std::ifstream fin(fname);

    ASSERT_(fin, "could not open scaffold file " + fname);

    // -------------------------------------------------------------------------
    // Step 2: Clear previously stored scaffold data
    // -------------------------------------------------------------------------

    segment.clear();
    nodes_by_id.clear();
    segment_index_by_id.clear();
    node_spatial_index.clear();
    node_spatial_bucket_size = 0.0;

    // -------------------------------------------------------------------------
    // Step 3: Read and validate file version
    // -------------------------------------------------------------------------

    double version = 0.0;

    fin >> version;

    ASSERT_(fin,
            "could not read scaffold file version from " + fname);

    if (version != 1.0 && version != 2.0) {
      ABORT_(
        "unsupported scaffold file version "
        + std::to_string(version)
        + " in " + fname
      );
    }

    // -------------------------------------------------------------------------
    // Step 4: Read scaffold nodes
    // -------------------------------------------------------------------------

    unsigned int number_of_nodes = 0;

    fin >> number_of_nodes;

    ASSERT_(fin,
            "could not read number of scaffold nodes from " + fname);

    ASSERT_(number_of_nodes > 0,
            "scaffold file contains no nodes: " + fname);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
      int node_id = 0;

      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      double radius = 0.0;

      if (version == 1.0) {
        /*
         * Legacy version 1.0:
         *
         * Node IDs are not stored explicitly. Assign one-based persistent IDs
         * using the node row order.
         */

        node_id = static_cast<int>(i) + 1;

        fin >> x >> y >> z >> radius;

      } else {
        /*
         * Version 2.0:
         *
         * Read the persistent one-based node ID directly from the file.
         */

        fin >> node_id >> x >> y >> z >> radius;
      }

      ASSERT_(
        fin,
        "failed reading scaffold node row "
        + std::to_string(i + 1)
        + " from " + fname
      );

      ASSERT_(
        node_id > 0,
        "scaffold node IDs must be positive in " + fname
      );

      ASSERT_(
        radius >= 0.0,
        "scaffold node radius cannot be negative for node ID "
        + std::to_string(node_id)
        + " in " + fname
      );

      ASSERT_(
        nodes_by_id.find(node_id) == nodes_by_id.end(),
        "duplicate scaffold node ID "
        + std::to_string(node_id)
        + " in " + fname
      );

      ScaffoldNode node;

      node.id = node_id;
      node.position = bdm::Double3{x, y, z};
      node.radius = radius;
      node.connected_node_ids.clear();

      nodes_by_id.emplace(node_id, node);
    }

    // -------------------------------------------------------------------------
    // Step 5: Read scaffold elements
    // -------------------------------------------------------------------------

    unsigned int number_of_elements = 0;

    fin >> number_of_elements;

    ASSERT_(fin,
            "could not read number of scaffold elements from " + fname);

    segment.reserve(number_of_elements);

    for (unsigned int i = 0; i < number_of_elements; ++i) {
      int element_id = 0;
      int node_id_1 = 0;
      int node_id_2 = 0;

      if (version == 1.0) {
        /*
         * Legacy version 1.0:
         *
         * Connectivity is stored using zero-based node indices.
         * Convert these indices to the one-based persistent IDs assigned above.
         */

        int legacy_node_index_1 = -1;
        int legacy_node_index_2 = -1;

        fin >> legacy_node_index_1 >> legacy_node_index_2;

        element_id = static_cast<int>(i) + 1;

        node_id_1 = legacy_node_index_1 + 1;
        node_id_2 = legacy_node_index_2 + 1;

      } else {
        /*
         * Version 2.0:
         *
         * Read the persistent element ID and persistent node IDs directly.
         */

        fin >> element_id >> node_id_1 >> node_id_2;
      }

      ASSERT_(
        fin,
        "failed reading scaffold element row "
        + std::to_string(i + 1)
        + " from " + fname
      );

      ASSERT_(
        element_id > 0,
        "scaffold element IDs must be positive in " + fname
      );

      ASSERT_(
        segment_index_by_id.find(element_id) == segment_index_by_id.end(),
        "duplicate scaffold element ID "
        + std::to_string(element_id)
        + " in " + fname
      );

      ASSERT_(
        node_id_1 > 0 && node_id_2 > 0,
        "scaffold element "
        + std::to_string(element_id)
        + " contains an invalid node ID in " + fname
      );

      ASSERT_(
        node_id_1 != node_id_2,
        "scaffold element "
        + std::to_string(element_id)
        + " connects node "
        + std::to_string(node_id_1)
        + " to itself in " + fname
      );

      auto node_1_it = nodes_by_id.find(node_id_1);
      auto node_2_it = nodes_by_id.find(node_id_2);

      ASSERT_(
        node_1_it != nodes_by_id.end(),
        "scaffold element "
        + std::to_string(element_id)
        + " references missing node ID "
        + std::to_string(node_id_1)
        + " in " + fname
      );

      ASSERT_(
        node_2_it != nodes_by_id.end(),
        "scaffold element "
        + std::to_string(element_id)
        + " references missing node ID "
        + std::to_string(node_id_2)
        + " in " + fname
      );

      const bdm::Double3& position_1 = node_1_it->second.position;
      const bdm::Double3& position_2 = node_2_it->second.position;

      const double element_length =
        L2norm(position_2 - position_1);

      ASSERT_(
        element_length > 0.0,
        "scaffold element "
        + std::to_string(element_id)
        + " has zero length in " + fname
      );

      const double element_radius =
        0.5 * (
          node_1_it->second.radius
          + node_2_it->second.radius
        );

      // -----------------------------------------------------------------------
      // Step 5a: Store node connectivity
      // -----------------------------------------------------------------------

      node_1_it->second.connected_node_ids.push_back(node_id_2);
      node_2_it->second.connected_node_ids.push_back(node_id_1);

      // -----------------------------------------------------------------------
      // Step 5b: Preserve the existing segment representation
      // -----------------------------------------------------------------------

      Segment geometric_segment;

      geometric_segment.element_id = element_id;
      geometric_segment.node_id_1 = node_id_1;
      geometric_segment.node_id_2 = node_id_2;

      geometric_segment.vertex_0 = position_1;
      geometric_segment.vertex_1 = position_2;

      geometric_segment.length = element_length;
      geometric_segment.radius = element_radius;

      // Store the vector position before adding the segment.
      const std::size_t segment_index = segment.size();

      segment.push_back(geometric_segment);

      segment_index_by_id.emplace(element_id, segment_index);
    }

    fin.close();
  }

  // ---------------------------------------------------------------------------
  // Persistent-ID lookup functions
  // ---------------------------------------------------------------------------

  inline bool HasNode(const int node_id) const {
    return nodes_by_id.find(node_id) != nodes_by_id.end();
  }

  inline const ScaffoldNode& GetNode(const int node_id) const {
    auto node_it = nodes_by_id.find(node_id);

    ASSERT_(
      node_it != nodes_by_id.end(),
      "could not find scaffold node ID "
      + std::to_string(node_id)
    );

    return node_it->second;
  }

  inline const bdm::Double3& GetNodePosition(const int node_id) const {
    return GetNode(node_id).position;
  }

  inline double GetNodeRadius(const int node_id) const {
    return GetNode(node_id).radius;
  }

  inline const std::vector<int>& GetConnectedNodeIds(
      const int node_id) const {
    return GetNode(node_id).connected_node_ids;
  }

  inline bool HasSegment(const int element_id) const {
  return segment_index_by_id.find(element_id)
         != segment_index_by_id.end();
  }

  inline const Segment& GetSegment(const int element_id) const {
    auto index_it = segment_index_by_id.find(element_id);

    ASSERT_(
      index_it != segment_index_by_id.end(),
      "could not find scaffold element ID "
      + std::to_string(element_id)
    );

    ASSERT_(
      index_it->second < segment.size(),
      "stored segment index is out of range for element ID "
      + std::to_string(element_id)
    );

    return segment[index_it->second];
  }

  inline
  void BuildNodeSpatialIndex(
      const double bucket_size) {
    /*
     * Function goal
     * -------------
     * Group persistent scaffold node IDs into uniform spatial buckets so later
     * radius searches inspect only nearby scaffold regions.
     */

    ASSERT_(
        bucket_size > 0.0,
        "Scaffold spatial-index bucket size must be positive"
    );

    ASSERT_(
        !nodes_by_id.empty(),
        "Cannot build a spatial index for a scaffold with no nodes"
    );

    node_spatial_index.clear();
    node_spatial_bucket_size = bucket_size;

    for (const auto& node_entry : nodes_by_id) {
      const int node_id =
          node_entry.first;

      const bdm::Double3& position =
          node_entry.second.position;

      const SpatialBucketKey bucket_key = {
          static_cast<int>(
              std::floor(position[0] / bucket_size)
          ),
          static_cast<int>(
              std::floor(position[1] / bucket_size)
          ),
          static_cast<int>(
              std::floor(position[2] / bucket_size)
          )
      };

      node_spatial_index[bucket_key].push_back(
          node_id
      );
    }

    // Keep node processing deterministic within every bucket.
    for (auto& bucket_entry : node_spatial_index) {
      std::sort(
          bucket_entry.second.begin(),
          bucket_entry.second.end()
      );
    }
  }

  inline
  std::vector<int> GetNodeIdsWithinRadius(
      const bdm::Double3& centre,
      const double search_radius) const {
    /*
     * Function goal
     * -------------
     * Return the persistent scaffold node IDs lying within a specified radius
     * of a point using the prebuilt spatial node index.
     */

    // -------------------------------------------------------------------------
    // Step 1: Validate the spatial index and search radius
    // -------------------------------------------------------------------------

    ASSERT_(
        search_radius > 0.0,
        "Scaffold node-radius search requires a positive radius"
    );

    ASSERT_(
        node_spatial_bucket_size > 0.0,
        "Scaffold node-radius search requires a built spatial index"
    );

    ASSERT_(
        !node_spatial_index.empty(),
        "Scaffold node-radius search encountered an empty spatial index"
    );

    const double radius_squared =
        search_radius * search_radius;

    // -------------------------------------------------------------------------
    // Step 2: Identify the centre bucket and search range
    // -------------------------------------------------------------------------

    const SpatialBucketKey centre_bucket = {
        static_cast<int>(
            std::floor(
                centre[0] / node_spatial_bucket_size
            )
        ),
        static_cast<int>(
            std::floor(
                centre[1] / node_spatial_bucket_size
            )
        ),
        static_cast<int>(
            std::floor(
                centre[2] / node_spatial_bucket_size
            )
        )
    };

    const int bucket_search_radius =
        static_cast<int>(
            std::ceil(
                search_radius / node_spatial_bucket_size
            )
        );

    // -------------------------------------------------------------------------
    // Step 3: Inspect only nearby spatial buckets
    // -------------------------------------------------------------------------

    std::vector<int> nearby_node_ids;

    for (int dx = -bucket_search_radius;
         dx <= bucket_search_radius;
         ++dx) {

      for (int dy = -bucket_search_radius;
           dy <= bucket_search_radius;
           ++dy) {

        for (int dz = -bucket_search_radius;
             dz <= bucket_search_radius;
             ++dz) {

          const SpatialBucketKey bucket_key = {
              centre_bucket.x + dx,
              centre_bucket.y + dy,
              centre_bucket.z + dz
          };

          const auto bucket_it =
              node_spatial_index.find(bucket_key);

          if (bucket_it == node_spatial_index.end()) {
            continue;
          }

          // -----------------------------------------------------------------
          // Step 4: Apply the exact spherical-distance check
          // -----------------------------------------------------------------

          for (const int node_id : bucket_it->second) {
            const bdm::Double3& node_position =
                this->GetNodePosition(node_id);

            const double dx_node =
                node_position[0] - centre[0];

            const double dy_node =
                node_position[1] - centre[1];

            const double dz_node =
                node_position[2] - centre[2];

            const double distance_squared =
                dx_node * dx_node
                + dy_node * dy_node
                + dz_node * dz_node;

            if (distance_squared <= radius_squared) {
              nearby_node_ids.push_back(node_id);
            }
          }
        }
      }
    }

    // Preserve deterministic processing independent of bucket traversal order.
    std::sort(
        nearby_node_ids.begin(),
        nearby_node_ids.end()
    );

    return nearby_node_ids;
  }

 public:
  // Scaffold nodes stored by persistent one-based node ID.
  std::unordered_map<int, ScaffoldNode> nodes_by_id;

  // Persistent element ID mapped to its position in the segment vector.
  std::unordered_map<int, std::size_t> segment_index_by_id;

  // Existing segment representation used by obstacle calculations.
  std::vector<ObstacleScaffold::Segment> segment;

  // Persistent node IDs grouped into uniform spatial buckets.
  std::unordered_map<
      SpatialBucketKey,
      std::vector<int>,
      SpatialBucketKeyHash
  > node_spatial_index;

  // Width of one spatial-index bucket.
  double node_spatial_bucket_size = 0.0;

};
// -----------------------------------------------------------------------------
class ObstacleBox : public Obstacle {
public:
  ObstacleBox() {}
  ~ObstacleBox() {}
  //
  inline
  void init(const std::string& t, const std::vector<bdm::Double3>& v) {
    this->type = t;
    bdm::Double3 A = v[0], B = v[1], C = v[2], D = v[3],
                 E = v[4], F = v[5], G = v[6], H = v[7];
    // set the vertices
    vertex_0 = v[0]; vertex_1 = v[1]; vertex_2 = v[2]; vertex_3 = v[3];
    vertex_4 = v[4]; vertex_5 = v[5]; vertex_6 = v[6]; vertex_7 = v[7];
    // local axes of the box
    laxis_0 = (B-A); length_0 = L2norm(laxis_0);
    laxis_1 = (D-A); length_1 = L2norm(laxis_1);
    laxis_2 = (E-A); length_2 = L2norm(laxis_2);
    laxis_0.Normalize();
    laxis_1.Normalize();
    laxis_2.Normalize();
    center = (A+G)*0.5; // geometric center
    center2face_0 = ((B+C+F+G)*0.25); // X>center
    center2face_1 = ((A+D+E+H)*0.25); // X<center
    center2face_2 = ((C+D+G+H)*0.25); // Y>center
    center2face_3 = ((A+B+E+F)*0.25); // Y<center
    center2face_4 = ((E+F+G+H)*0.25); // Z>center
    center2face_5 = ((A+B+C+D)*0.25); // Z<center
    if ("box/inside"==this->type)
      {
        normal2face_0 = center2face_0 - center;
        normal2face_1 = center2face_1 - center;
        normal2face_2 = center2face_2 - center;
        normal2face_3 = center2face_3 - center;
        normal2face_4 = center2face_4 - center;
        normal2face_5 = center2face_5 - center;
      }
    else if ("box/outside"==this->type)
      {
        normal2face_0 = center - center2face_0;
        normal2face_1 = center - center2face_1;
        normal2face_2 = center - center2face_2;
        normal2face_3 = center - center2face_3;
        normal2face_4 = center - center2face_4;
        normal2face_5 = center - center2face_5;
      }
    else
      ABORT_("an exception is caught");
    //
    normal2face_0.Normalize();
    normal2face_1.Normalize();
    normal2face_2.Normalize();
    normal2face_3.Normalize();
    normal2face_4.Normalize();
    normal2face_5.Normalize();
  }
  //
  inline
  bool is_inside(const bdm::Double3& p) const {
    const bdm::Double3 V = p - center; // direction vector from point to center
    const double Va_0 = (2.0*fabs(V*laxis_0)),
                 Va_1 = (2.0*fabs(V*laxis_1)),
                 Va_2 = (2.0*fabs(V*laxis_2));
    return (Va_0<=length_0 && Va_1<=length_1 && Va_2<=length_2);
  }
//
public:
  bdm::Double3 vertex_0, vertex_1, vertex_2, vertex_3,
               vertex_4, vertex_5, vertex_6, vertex_7;
  bdm::Double3 laxis_0,  laxis_1,  laxis_2;
  double      length_0, length_1, length_2;
  bdm::Double3 center;
  bdm::Double3 center2face_0, center2face_1, center2face_2,
               center2face_3, center2face_4, center2face_5;
  bdm::Double3 normal2face_0, normal2face_1, normal2face_2,
               normal2face_3, normal2face_4, normal2face_5;
};
// -----------------------------------------------------------------------------
class ObstacleSphere : public Obstacle {
public:
  ObstacleSphere() {}
  ~ObstacleSphere() {}
  //
  inline
  void init(const std::string& t, const bdm::Double3& c, double r) {
    this->type = t;
    center = c; // geometric center
    radius = r; // ...and radius
  }
  //
  inline
  bool is_inside(const bdm::Double3& p) const {
    const bdm::Double3 V = p - center; // direction vector from point to center
    const double r = L2norm(V);
    return (r<=radius);
  }
//
public:
  bdm::Double3 center;
  double radius;
};
// -----------------------------------------------------------------------------
class ObstacleSTL : public Obstacle {
public:
  struct Triangle {
    bdm::Double3 vertex_0, vertex_1, vertex_2;
    bdm::Double3 center, normal, inside;
  };
//
public:
  ObstacleSTL() {}
  ~ObstacleSTL() {}
  //
  inline
  void init(const std::string& t, const std::string& fname) {
    this->type = t;
    // now open the STL file to process
    std::ifstream fin(fname);
    ASSERT_(fin.good(),"file \""+fname+"\" cannot be accessed");
    // enforce to clear memory
    triangle.clear();
    //
    std::string s;
    // read STL header and confirm it's valid
    {
      std::getline(fin, s);
      const std::vector<std::string> parsed_s = extract_words_vector(s);
      if ("solid"!=parsed_s[0])
        ABORT_("could not parse STL file here: "+s);
    }
    //
    while ( true )
      {
        bdm::Double3 v0, v1, v2, n;
        // read the normal or check if you have reached the end of the STL file
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          // check if you have almost finished reading the file
          if ("endsolid"==parsed_s[0]) break;
          // ...otherwise read this triangle, first the normal vector
          else if ("facet"!=parsed_s[0] || "normal"!=parsed_s[1] || 5!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
          n = { std::stod(parsed_s[2]), std::stod(parsed_s[3]), std::stod(parsed_s[4]) };
        }
        //
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          if ("outer"!=parsed_s[0] || "loop"!=parsed_s[1] || 2!=parsed_s.size())
            bdm::Log::Fatal(std::string(__FILE__),"error @line "+std::to_string(__LINE__)+" / "+s);
        }
        // read the 1st vertex
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          if ("vertex"!=parsed_s[0] || 4!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
          v0 = { std::stod(parsed_s[1]), std::stod(parsed_s[2]), std::stod(parsed_s[3]) };
        }
        // read the 2nd vertex
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          if ("vertex"!=parsed_s[0] || 4!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
          v1 = { std::stod(parsed_s[1]), std::stod(parsed_s[2]), std::stod(parsed_s[3]) };
        }
        // read the 3rd vertex
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          if ("vertex"!=parsed_s[0] || 4!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
          v2 = { std::stod(parsed_s[1]), std::stod(parsed_s[2]), std::stod(parsed_s[3]) };
        }
        //
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          if ("endloop"!=parsed_s[0] || 1!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
        }
        //
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          if ("endfacet"!=parsed_s[0] || 1!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
        }
        // completed reading this triangle, now process the data
        ObstacleSTL::Triangle tri3;
        tri3.vertex_0 = v0;
        tri3.vertex_1 = v1;
        tri3.vertex_2 = v2;
        // geometric center
        tri3.center = (v0+v1+v2)/3.0;
        // outward unit normal vector
        if (!normalize(n, tri3.normal))
          ABORT_("could not process normal of a facet in STL file");
        // internal point located right below the triangle
        {
          std::vector<double> len = { L2norm(tri3.vertex_0-tri3.center) ,
                                      L2norm(tri3.vertex_1-tri3.center) ,
                                      L2norm(tri3.vertex_2-tri3.center) };
          const double W = (*std::min_element(len.begin(), len.end()))
                         * 1.0e-3;
          tri3.inside = tri3.center - tri3.normal * W;
        }
        // load this triangle into member container
        triangle.push_back(tri3);
      }
    // now close the file stream
    fin.close();
  }
//
public:
  // list of all triangles
  std::vector<ObstacleSTL::Triangle> triangle;
};
// =============================================================================
class SimulationObstacles {
public:
  SimulationObstacles() {}
  ~SimulationObstacles() {}
  //
  inline
  void clear() { box.clear(); sphere.clear(); surface.clear(); scaffold.clear(); }
//
public:
  std::vector<ObstacleBox> box;
  std::vector<ObstacleSphere> sphere;
  std::vector<ObstacleSTL> surface;
  std::vector<ObstacleScaffold> scaffold;
};
// =============================================================================
#endif // _OBSTACLES_H_
// =============================================================================
