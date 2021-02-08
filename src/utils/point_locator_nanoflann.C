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
#include "libmesh/int_range.h" // make_range

// C++ includes
#include <array>

namespace libMesh
{

PointLocatorNanoflann::PointLocatorNanoflann (const MeshBase & mesh,
                                              const PointLocatorBase * master) :
  PointLocatorBase (mesh, master),
  _out_of_mesh_mode(false),
  _num_results(3)
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

PointLocatorNanoflann::NanoflannResult
PointLocatorNanoflann::kd_tree_find_neighbors(const Point & p) const
{
  // We are searching for the Point(s) closest to Point p.
  //
  // TODO: The kd_tree's findNeighbors() routine needs a pointer to
  // Real of length LIBMESH_DIM. It might be convenient if libMesh
  // Points had a data() member that provided this, for now we just
  // copy the coordinates into a std::array of the appropriate size.
  std::array<Real, LIBMESH_DIM> query_pt;
  for (int i=0; i<LIBMESH_DIM; ++i)
    query_pt[i] = p(i);

  // To catch values returned by the search
  std::vector<std::size_t> ret_index(_num_results);
  std::vector<Real> out_dist_sqr(_num_results);
  nanoflann::KNNResultSet<Real> result_set(_num_results);

  // Initialize the result_set
  result_set.init(ret_index.data(), out_dist_sqr.data());

  // Do the search
  // We leave all the SearchParams ctor args on their defaults:
  // int ignored == 32
  // float eps == 0
  // bool sorted == true
  _kd_tree->findNeighbors(result_set, query_pt.data(), nanoflann::SearchParams());

  return std::make_tuple(ret_index, out_dist_sqr, result_set);
}

const Elem *
PointLocatorNanoflann::operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator()", "PointLocatorNanoflann");

  // Do the search
  auto t = this->kd_tree_find_neighbors(p);

  // References to the tuple contents.
  // TODO: In C++17 we can use structured bindings to replace this.
  const auto & ret_index = std::get<0>(t);
  const auto & out_dist_sqr = std::get<1>(t);
  const auto & result_set = std::get<2>(t);

  // Loop over the list of candidate centroids, returning the Elem associated with
  // the centroid to which the Point is both nearest, and contained by (to within the
  // _contains_point_tol).
  for (auto r : make_range(result_set.size()))
    {
      // For indexing into original data structures.
      unsigned int elem_id = ret_index[r];

      // Debugging: print the results
      // libMesh::out << "Centroid/Elem id = " << elem_id
      //              << ", dist^2 = " << out_dist_sqr[r]
      //              << std::endl;

      const Elem * candidate_elem = _mesh.elem_ptr(elem_id);

      // Before we even check whether the candidate Elem actually
      // contains the Point, we may need to check whether the
      // candidate Elem is from an allowed subdomain.  If the
      // candidate Elem is not from an allowed subdomain, we continue
      // to the next one.
      if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
        continue;

      // If we made it here, then the candidate Elem is from an
      // allowed subdomain, so let's next check whether it contains
      // the point. If the user set a custom tolerance, then we
      // actually check close_to_point() rather than contains_point(),
      // since this latter function warns about using non-default
      // tolerances, but otherwise does the same test.
      bool inside = _use_contains_point_tol ?
        candidate_elem->close_to_point(p, _contains_point_tol) :
        candidate_elem->contains_point(p);

      // If the point is inside an Elem from an allowed subdomain, we are done.
      if (inside)
        return candidate_elem;
    }

  // If we made it here, then either all the candidate elements were
  // from non-allowed subdomains, or the Point was not inside _any_
  // candidate Elem, so choose a return value based on the
  // _out_of_mesh_mode flag.
  if (_out_of_mesh_mode)
    return nullptr;
  else
    libmesh_error_msg("Point was not contained within the Elem (to within the required tolerance) "
                      "whose centroid it was closest to, and _out_of_mesh_mode was not enabled.");
}


void
PointLocatorNanoflann::operator() (const Point & p,
                                   std::set<const Elem *> & candidate_elements,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator() returning set", "PointLocatorNanoflann");

  // Do the search
  auto t = this->kd_tree_find_neighbors(p);

  // References to the tuple contents.
  // TODO: In C++17 we can use structured bindings to replace this.
  const auto & ret_index = std::get<0>(t);
  const auto & out_dist_sqr = std::get<1>(t);
  const auto & result_set = std::get<2>(t);

  // Loop over the list of candidate centroids, adding each one which
  // passes the test to the candidate_elements set.
  for (auto r : make_range(result_set.size()))
    {
      // For indexing into original data structures.
      unsigned int elem_id = ret_index[r];

      // Debugging: print the results
      // libMesh::out << "Centroid/Elem id = " << elem_id
      //              << ", dist^2 = " << out_dist_sqr[r]
      //              << std::endl;

      const Elem * candidate_elem = _mesh.elem_ptr(elem_id);

      // Before we even check whether the candidate Elem actually
      // contains the Point, we may need to check whether the
      // candidate Elem is from an allowed subdomain.  If the
      // candidate Elem is not from an allowed subdomain, we continue
      // to the next one.
      if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
        continue;

      // If we made it here, then the candidate Elem is from an
      // allowed subdomain, so let's next check whether it contains
      // the point. If the user set a custom tolerance, then we
      // actually check close_to_point() rather than contains_point(),
      // since this latter function warns about using non-default
      // tolerances, but otherwise does the same test.
      bool inside = _use_contains_point_tol ?
        candidate_elem->close_to_point(p, _contains_point_tol) :
        candidate_elem->contains_point(p);

      // If the point is inside an Elem from an allowed subdomain, add
      // it to the list.
      if (inside)
        candidate_elements.insert(candidate_elem);
    }
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
  // We only consider LIBMESH_DIM dimensional KD-Trees, so just make
  // sure we were called consistently.
  libmesh_assert(size == LIBMESH_DIM);

  // Construct a libmesh Point object from the LIBMESH_DIM-dimensional
  // input object, p1.
  Point point1;
  for (int i=0; i<LIBMESH_DIM; ++i)
    point1(i) = p1[i];

  // Compute Euclidean distance, squared
  return (point1 - _centroids[idx_p2]).norm_sq();
}



PointLocatorNanoflann::coord_t
PointLocatorNanoflann::kdtree_get_pt(const std::size_t idx, int dim) const
{
  libmesh_assert_less (idx, _centroids.size());
  libmesh_assert_less (dim, LIBMESH_DIM);

  return _centroids[idx](dim);
}

} // namespace libMesh

#endif
