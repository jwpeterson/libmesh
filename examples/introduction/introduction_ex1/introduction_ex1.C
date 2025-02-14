// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>Introduction Example 1 - Creation of a Mesh Object</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This is the first example program.  It simply demonstrates
// how to create a mesh object.  A mesh is read from file,
// information is printed to the screen, and the mesh is then
// written.

// C++ include files that we need
#include <iostream>
// Functions to initialize the library.
#include "libmesh/libmesh.h"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

int main (int argc, char ** argv)
{
  // Initialize the library.  This is necessary because the library
  // may depend on a number of other libraries (i.e. MPI and PETSc)
  // that require initialization before use.  When the LibMeshInit
  // object goes out of scope, other libraries and resources are
  // finalized.
  LibMeshInit init (argc, argv);

  // Check for proper usage. The program is designed to be run
  // as follows:
  // ./ex1 -d DIM input_mesh_name [-o output_mesh_name]
  // where [output_mesh_name] is an optional parameter giving
  // a filename to write the mesh into.
  libmesh_error_msg_if(argc < 4, "Usage: " << argv[0] << " -d 2 in.mesh [-o out.mesh]");

  // Get the dimensionality of the mesh from the "-d" argument
  const unsigned int dim =
    libMesh::command_line_next("-d", libMesh::invalid_uint);
  libmesh_error_msg_if(dim > 3, "Usage: " << argv[0] << " -d 2 in.mesh [-o out.mesh]");

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  // We may need XDR support compiled in to read binary .xdr files
  const std::string input_filename = argv[3];
#ifndef LIBMESH_HAVE_XDR
  libmesh_example_requires(input_filename.rfind(".xdr") >=
                           input_filename.size(), "XDR support");
#endif

  // Read the input mesh.
  mesh.read (input_filename);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Write the output mesh if the user specified an
  // output file name.
  if (libMesh::on_command_line("-o"))
    {
      // We may need XDR support compiled in to read binary .xdr files
      const std::string output_filename =
        libMesh::command_line_next("-o",std::string());
#ifndef LIBMESH_HAVE_XDR
      libmesh_example_requires(output_filename.rfind(".xdr") >=
                               output_filename.size(), "XDR support");
#endif

      mesh.write (output_filename);
    }

  // All done.  libMesh objects are destroyed here.  Because the
  // LibMeshInit object was created first, its destruction occurs
  // last, and it's destructor finalizes any external libraries and
  // checks for leaked memory.
  return 0;
}
