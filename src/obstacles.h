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
  struct Segment {
    bdm::Double3 vertex_0, vertex_1;
    double length, radius;
  };
//
public:
  ObstacleScaffold() {}
  ~ObstacleScaffold() {}
  //
  inline
  void init(const std::string& t, const std::string& fname) {
    this->type = t;
    // now open the file to process
    std::ifstream fin(fname);
    ASSERT_(fin,"could not open file "+fname);
    // enforce to clear memory
    segment.clear();
    //
    unsigned int n_pnt, n_segm;
    std::vector<bdm::Double3> vertices;
    std::vector<double> radii;
    //
    fin >> n_pnt;
    for (unsigned int i=0; i<n_pnt; i++)
      {
        double x, y, z, r;
        fin >> x >> y >> z >> r;
        vertices.push_back({x,y,z});
        radii.push_back(r);
      }
    //
    fin >> n_segm;
    segment.resize(n_segm);
    for (unsigned int i=0; i<n_segm; i++)
      {
        int n0, n1;
        fin >> n0 >> n1;
        segment[i].vertex_0 = vertices[n0];
        segment[i].vertex_1 = vertices[n1];
        segment[i].length = L2norm(vertices[n1]-vertices[n0]);
        double rm = 0.5*(radii[n0]+radii[n1]);
        segment[i].radius = rm;
      }
    // now close the file stream
    fin.close();
  }
  //
public:
  // list of all segments
  std::vector<ObstacleScaffold::Segment> segment;
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
