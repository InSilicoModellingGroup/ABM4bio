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
#ifndef _IOFLUX_H_
#define _IOFLUX_H_
// =============================================================================
#include "./global.h"
#include "core/container/math_array.h"
// =============================================================================
class SurfaceSTL {
public:
  struct Triangle {
    bdm::Double3 vertex_0, vertex_1, vertex_2;
    bdm::Double3 normal;
  };
//
public:
  SurfaceSTL() {}
  ~SurfaceSTL() {}
  //
  inline
  void read(const std::string& fname) {
    // open the STL file to process
    std::ifstream fin(fname);
    ASSERT_(fin,"could not open file "+fname);
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
          //
          // ...otherwise read this triangle, first the normal vector
          if ("facet"!=parsed_s[0] || "normal"!=parsed_s[1] || 5!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
          n = { std::stod(parsed_s[2]), std::stod(parsed_s[3]), std::stod(parsed_s[4]) };
        }
        //
        {
          std::getline(fin, s);
          const std::vector<std::string> parsed_s = extract_words_vector(s);
          if ("outer"!=parsed_s[0] || "loop"!=parsed_s[1] || 2!=parsed_s.size())
            ABORT_("could not parse STL file here: "+s);
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
        SurfaceSTL::Triangle tri3;
        tri3.vertex_0 = v0;
        tri3.vertex_1 = v1;
        tri3.vertex_2 = v2;
        // outward unit normal vector
        tri3.normal = normalize(n);
        // load this triangle into member container
        triangle.push_back(tri3);
      }
    // now close the file stream
    fin.close();
  }
//
public:
  // list of all triangles
  std::vector<SurfaceSTL::Triangle> triangle;
};
// =============================================================================
class SurfaceTraits {
public:
  SurfaceTraits() : time_step(1), phenotype(0) {}
  ~SurfaceTraits() {}
  //
public:
  // time interval between successive cell flux insertions
  int time_step;
  // list of time parameters defining the cell flux pattern
  std::vector<double> time_parameters;
  //
  int phenotype;
};
// =============================================================================
class Surface : public SurfaceSTL, public SurfaceTraits {
public:
  Surface() {}
  ~Surface() {}
  //
  inline
  void init(int DT, const std::vector<double>& F, int p, const std::string& fname) {
    this->time_step = DT;
    this->time_parameters = F;
    this->phenotype = p;
    // alse read the STL file
    this->read(fname);
  }
};
// =============================================================================
class SimulationIOFlux {
public:
  SimulationIOFlux() {}
  ~SimulationIOFlux() {}
  //
  inline
  void clear() { surface.clear(); }
//
public:
  std::vector<Surface> surface;
};
// =============================================================================
#endif // _IOFLUX_H_
// =============================================================================
