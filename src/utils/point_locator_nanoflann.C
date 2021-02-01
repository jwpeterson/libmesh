// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local Includes
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/tree.h"

namespace libMesh
{

PointLocatorNanoflann::PointLocatorNanoflann (const MeshBase & mesh,
                                              const PointLocatorBase * master) :
  PointLocatorBase (mesh, master),
  _out_of_mesh_mode(false)
{
  this->init();
}



PointLocatorNanoflann::~PointLocatorNanoflann () = default;



void PointLocatorNanoflann::clear ()
{
}



void PointLocatorNanoflann::init ()
{
  // TODO

  this->_initialized = true;
}

const Elem *
PointLocatorNanoflann::operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator()", "PointLocatorNanoflann");

  // TODO
}


void
PointLocatorNanoflann::operator() (const Point & p,
                                   std::set<const Elem *> & candidate_elements,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator() returning set", "PointLocatorNanoflann");

  // TODO
}


void PointLocatorNanoflann::enable_out_of_mesh_mode ()
{
  // Out-of-mesh mode should now work properly even on meshes with
  // non-affine elements.
  _out_of_mesh_mode = true;
}


void PointLocatorNanoflann::disable_out_of_mesh_mode ()
{
  _out_of_mesh_mode = false;
}



} // namespace libMesh
