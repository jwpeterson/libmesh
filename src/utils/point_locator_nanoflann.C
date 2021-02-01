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


#include "libmesh/libmesh_config.h"

// This class is not defined unless libmesh is built with Nanoflann support
#ifdef LIBMESH_HAVE_NANOFLANN

// libmesh includes
#include "libmesh/point_locator_nanoflann.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"

// C++ includes
#include <array>

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



void
PointLocatorNanoflann::clear ()
{
  this->_initialized = false;
  this->_out_of_mesh_mode = false;
  _centroids.clear();
  _kd_tree.reset();
}



void
PointLocatorNanoflann::init ()
{
  // TODO: for the moment we ignore whether the "_master" flag is set or not.
  if (!_initialized)
    {
      // TODO: fill in the _centroids array


      // Construct the KD-Tree tree
      _kd_tree = libmesh_make_unique<kd_tree_t>
        (3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10));

      _kd_tree->buildIndex();

      // We are initialized now
      this->_initialized = true;
    }
}

const Elem *
PointLocatorNanoflann::operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator()", "PointLocatorNanoflann");

  // We are searching for Points close to the rotated centroid.
  std::array<Real, 3> query_pt = {p(0), p(1), p(2)};

  // The number of results we want to get back from Nanoflann.
  // Hopefully there is only 1 "obvious" nearest node, but it is
  // good to verify this.
  const std::size_t num_results = 3;

  // To catch values returned by the search
  std::array<std::size_t, num_results> ret_index;
  std::array<Real, num_results> out_dist_sqr;
  nanoflann::KNNResultSet<Real> result_set(num_results);

  // Initialize the result_set
  result_set.init(ret_index.data(), out_dist_sqr.data());

  // Do the search
  _kd_tree->findNeighbors(result_set, query_pt.data(), nanoflann::SearchParams(10));

  // Debugging: print the results
  for (unsigned r=0; r<result_set.size(); ++r)
    {
      // For indexing into original data structures.
      unsigned int i = ret_index[r];

      libMesh::out << "ret_index = " << i
                   << ", dist^2 = " << out_dist_sqr[r]
                   << std::endl;
    }

  // TODO: we found the closest Node, but that Node will be connected
  // to several elements, and we would need to know which connected
  // Elem the query point is actually inside, which would require a
  // bunch of contains_point() checks. So for this reason it is better
  // to make the KD-Tree from a list of Points which correspond to the
  // Elem centroids.
  return nullptr;
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


void
PointLocatorNanoflann::enable_out_of_mesh_mode ()
{
  // Out-of-mesh mode should now work properly even on meshes with
  // non-affine elements.
  _out_of_mesh_mode = true;
}


void
PointLocatorNanoflann::disable_out_of_mesh_mode ()
{
  _out_of_mesh_mode = false;
}

//
// Required Nanoflann APIs
//

std::size_t PointLocatorNanoflann::kdtree_get_point_count() const
{
  return _mesh.n_nodes();
}



PointLocatorNanoflann::coord_t
PointLocatorNanoflann::kdtree_distance(const coord_t * p1,
                                       const std::size_t idx_p2,
                                       std::size_t size) const
{
  // Construct a libmesh Point object from the input coord_t.  This
  // assumes LIBMESH_DIM==3.
  Point point1(p1[0],
               size > 1 ? p1[1] : 0.,
               size > 2 ? p1[2] : 0.);

  // Get the referred-to point from the Mesh
  const Point & point2 = _mesh.point(idx_p2);

  // Compute Euclidean distance, squared
  return (point1 - point2).norm_sq();
}



PointLocatorNanoflann::coord_t
PointLocatorNanoflann::kdtree_get_pt(const std::size_t idx, int dim) const
{
  libmesh_assert_less (idx, _mesh.n_nodes());
  libmesh_assert_less (dim, LIBMESH_DIM);

  return _mesh.point(idx)(dim);
}

} // namespace libMesh

#endif
