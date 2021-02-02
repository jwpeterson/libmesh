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
      // Fill in the _centroids array with the active element
      // centroids.
      //
      // TODO: Currently we assume replicated Mesh, but one way to
      // parallelize the Nanoflann PointLocator would be to create
      // separate KD-Trees for each partition, do the search locally,
      // and then gather the results. Since the KD-Tree seems to be
      // pretty fast and light-weight, we'll have to profile to see
      // whether this optimization is necessary...
      _centroids.clear();
      _centroids.reserve(_mesh.n_elem());
      for (const auto & elem : _mesh.active_element_ptr_range())
        _centroids.push_back(elem->centroid());

      // Construct the KD-Tree
      _kd_tree = libmesh_make_unique<kd_tree_t>
        (LIBMESH_DIM, *this, nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10));

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

  // We assume that the Point, if it is contained in _any_ element,
  // will be contained in the one whose centroid it is closest to.
  // If this assumption ever wrong?
  const Elem * elem = _mesh.elem_ptr(ret_index[0]);

  // Before we even check whether the Elem actually contains the
  // Point, we should check whether the Elem is from an allowed
  // subdomain.
  if (allowed_subdomains)
    {
      bool allowed = allowed_subdomains->count(elem->subdomain_id());

      // If the Elem with the nearest centroid is not from an allowed
      // subdomain, then we return either nullptr (when
      // out-of-mesh-mode is enabled) or throw an error.
      if (!allowed)
        {
          if (_out_of_mesh_mode)
            return nullptr;
          else
            libmesh_error_msg("The Elem closest to the searched-for Point is not in an "
                              "allowed subdomain, and _out_of_mesh_mode was not enabled.");
        }
    }

  // If we are using a tolerance pass it to the contains_point() call.
  bool inside = _use_contains_point_tol ?
    elem->contains_point(p, _contains_point_tol) :
    elem->contains_point(p);

  if (inside)
    return elem;

  else // outside
    {
      if (_out_of_mesh_mode)
        return nullptr;
      else
        libmesh_error_msg("Point was not contained within the Elem (to within the required tolerance) "
                          "whose centroid it was closest to, and _out_of_mesh_mode was not enabled.");
    }

  // We'll never get here, see above
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
  libmesh_not_implemented();
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
  return _centroids.size();
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

  // Compute Euclidean distance, squared
  return (point1 - _centroids[idx_p2]).norm_sq();
}



PointLocatorNanoflann::coord_t
PointLocatorNanoflann::kdtree_get_pt(const std::size_t idx, int dim) const
{
  libmesh_assert_less (idx, _mesh.n_nodes());
  libmesh_assert_less (dim, LIBMESH_DIM);

  return _centroids[idx](dim);
}

} // namespace libMesh

#endif
