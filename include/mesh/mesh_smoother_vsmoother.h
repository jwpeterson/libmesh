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



#ifndef LIBMESH_MESH_SMOOTHER_VSMOOTHER_H
#define LIBMESH_MESH_SMOOTHER_VSMOOTHER_H

#include "libmesh/libmesh_config.h"
#if defined(LIBMESH_ENABLE_VSMOOTHER)

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_smoother.h"
#include "libmesh/variational_smoother_system.h"
#include "petsc_diff_solver.h"

// C++ Includes
#include <cstddef>
#include <vector>
#include <map>
#include <fstream>

namespace libMesh
{

// Forward declarations
class UnstructuredMesh;

/**
 * This is an implementation of Larisa Branets' smoothing algorithms.
 * The initial implementation was done by her, the adaptation to
 * libmesh was completed by Derek Gaston.  The code was heavily
 * refactored into something more closely resembling C++ by John
 * Peterson in 2014.
 *
 * Here are the relevant publications:
 * 1) L. Branets, G. Carey, "Extension of a mesh quality metric for
 * elements with a curved boundary edge or surface",
 * Journal of Computing and Information Science in Engineering, vol. 5(4), pp.302-308, 2005.
 *
 * 2) L. Branets, G. Carey, "A local cell quality metric and variational grid
 * smoothing algorithm", Engineering with Computers, vol. 21, pp.19-28, 2005.
 *
 * 3) L. Branets, "A variational grid optimization algorithm based on a local
 * cell quality metric", Ph.D. thesis, The University of Texas at Austin, 2005.
 *
 * \author Derek R. Gaston
 * \date 2006
 */
class VariationalMeshSmoother : public MeshSmoother
{
public:

  /**
   * Simple constructor to use for smoothing purposes
   */
  VariationalMeshSmoother(UnstructuredMesh & mesh,
                          Real dilation_weight=0.5,
                          const bool preserve_subdomain_boundaries=true);

  /**
   * Destructor.
   */
  virtual ~VariationalMeshSmoother() = default;

  /**
   * Redefinition of the smooth function from the
   * base class.  All this does is call the smooth
   * function in this class which takes an int, using
   * a default value of 1.
   */
  virtual void smooth() override {this->smooth(1); }

  /**
   * The actual smoothing function, gets called whenever
   * the user specifies an actual number of smoothing
   * iterations.
   */
  void smooth(unsigned int n_iterations);


private:

  /**
   * Smoother control variables
   */
  const Real _dilation_weight;

  /// Whether subdomain boundaries are subject to change via smoothing
  const bool _preserve_subdomain_boundaries;
};

} // namespace libMesh

#endif // defined(LIBMESH_ENABLE_VSMOOTHER)

#endif // LIBMESH_MESH_SMOOTHER_VSMOOTHER_H
