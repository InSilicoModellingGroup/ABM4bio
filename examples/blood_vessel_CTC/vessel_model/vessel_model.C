// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ include files that we need
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <stdio.h>
#include <math.h>

#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/newmark_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


// The main program
int main (int argc, char** argv)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init (argc, argv);

  // Check for proper usage.
  libmesh_error_msg_if(argc < 2, "Usage: " << argv[0] << " [meshfile]");

  // Tell the user what we are doing.
  libMesh::out << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  // Get the name of the mesh file
  // from the command line.
  std::string mesh_file = argv[1];
  libMesh::out << "Mesh file is: " << mesh_file << std::endl;

  // Create a ReplicatedMesh object, with dimension to be overridden
  // later, distributed across the default MPI communicator.
  Mesh mesh(init.comm());

  // Read the meshfile specified on the command line.
  GmshIO(mesh).read(mesh_file);

  // Print information about the mesh to the screen.
  mesh.prepare_for_use();
  mesh.print_info();

  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> dis__HUVEC(24.0, 35.0);
  std::uniform_real_distribution<> dis__RBC(7.5, 8.5);

  const double RAD = 40.0;
  const double HEI = 100.0;

  std::uniform_real_distribution<> r__RBC(0.05*RAD, 0.95*RAD);
  std::uniform_real_distribution<> z__RBC(-0.95*HEI, +0.95*HEI);
  std::uniform_real_distribution<> w__RBC(-M_PI, M_PI);

  std::ofstream fout;

  fout.open("initial_HUVEC.dat");
  fout << mesh.n_active_elem() << std::endl;
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      //libmesh_assert_equal_to(ElemType::TRI3, elem->type());

      const Point p0 = elem->point(0);
      const Point p1 = elem->point(1);
      const Point p2 = elem->point(2);
      const Point p10 = p1 - p0;
      const Point p20 = p2 - p0;

      Point e2 = p10.cross(p20).unit(),
            e0 = p10.unit(),
            e1 = e2.cross(e0);

      e0 *= 1.00 *2;
      e1 *= 1.00 *2;
      e2 *= 0.25 *2;

      double dia = dis__HUVEC(rng) /3;

      Point c = elem->vertex_average();

      fout << c(0) << ' ' << c(1) << ' ' << c(2)
           << ' ' << dia
           << ' ' << e0(0) << ' ' << e0(1) << ' ' << e0(2)
           << ' ' << e1(0) << ' ' << e1(1) << ' ' << e1(2)
           << ' ' << e2(0) << ' ' << e2(1) << ' ' << e2(2)
           << ' ' << 0
           << std::endl;
    }
  fout.close();

  fout.open("initial_RBC.dat");
  fout << 1290 << std::endl;
  for (size_t l=0; l<1290; l++)
    {
      double r = r__RBC(rng), w = w__RBC(rng),
             x = r*cos(w), y = r*sin(w), z = z__RBC(rng);

      double dia = dis__RBC(rng);

      Point c(x,y,z);

      fout << c(0) << ' ' << c(1) << ' ' << c(2)
           << ' ' << dia
           << ' ' << 1.0 << ' ' << 0.0 << ' ' << 0.0
           << ' ' << 0.0 << ' ' << 1.0 << ' ' << 0.0
           << ' ' << 0.0 << ' ' << 0.0 << ' ' << 1.0
           << ' ' << 0
           << std::endl;
    }
  fout.close();

  // All done.
  return 0;
}
