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
#include "libmesh/utility.h" // libmesh_map_find
#include "libmesh/mesh_tools.h" // create_local_bounding_box

// C++ includes
#include <array>

namespace libMesh
{

PointLocatorNanoflann::PointLocatorNanoflann (const MeshBase & mesh,
                                              const PointLocatorBase * master) :
  PointLocatorBase (mesh, master),
  _out_of_mesh_mode(false),
  _initial_num_results(16),
  _max_num_results(512),
  _hmax(0.)
{
  this->init();
}



PointLocatorNanoflann::~PointLocatorNanoflann () = default;



void
PointLocatorNanoflann::clear ()
{
  this->_initialized = false;
  this->_out_of_mesh_mode = false;
  _ids.clear();
  _centroids.clear();
  _kd_tree.reset();
}



void
PointLocatorNanoflann::init ()
{
  LOG_SCOPE("init()", "PointLocatorNanoflann");

  // TODO: for the moment we ignore whether the "_master" flag is set or not.
  if (!_initialized)
    {
      // Fill in the _centroids data structure with active, local
      // element centroids.
      _ids.clear();
      _centroids.clear();

      // We can either reserve exactly the right amount of space or
      // let push_back() take care of it, not sure what would be
      // faster actually, since it takes some time to count the number
      // of active, local, etc. elements.
      //
      // Note: my current idea is that we will include all "active"
      // elements in the KD-Tree, instead of just local ones, since
      // that should result in fewer "failed" searches, which
      // currently are quite expensive for the Nanoflann PointLocator.
      auto n_active_elem = _mesh.n_active_elem();
      _ids.reserve(n_active_elem);
      _centroids.reserve(n_active_elem);

      for (const auto & elem : _mesh.active_element_ptr_range())
        {
          _ids.push_back(elem->id());
          _centroids.push_back(elem->centroid());

          // Only keep track of local elements' hmax.
          //
          // TODO: if we switch to building local, active KD-Trees, we
          // can drop this if-statement.
          if (elem->processor_id() == _mesh.comm().rank())
              _hmax = std::max(elem->hmax(), _hmax);
        }

      // Debugging:
      libMesh::out << "Local Elem hmax() = " << _hmax << std::endl;

      // Construct the KD-Tree
      _kd_tree = libmesh_make_unique<kd_tree_t>
        (LIBMESH_DIM, *this, nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10));

      _kd_tree->buildIndex();

      // A BoundingBox for the local elements
      _local_bbox = MeshTools::create_local_bounding_box (_mesh);

      // We are initialized now
      this->_initialized = true;
    }
}

nanoflann::KNNResultSet<Real>
PointLocatorNanoflann::kd_tree_find_neighbors(const Point & p,
                                              std::size_t num_results) const
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

  // Allocate storage for the indices and distances
  _ret_index.resize(num_results);
  _out_dist_sqr.resize(num_results);

  // nanoflann::KNNResultSet cannot be resized/reused easily, I think
  // it is just meant to be re-created for each search.
  nanoflann::KNNResultSet<Real> result_set(num_results);

  // Initialize the result_set
  result_set.init(_ret_index.data(), _out_dist_sqr.data());

  // Do the search
  // We leave all the SearchParams ctor args on their defaults:
  // int ignored == 32
  // float eps == 0
  // bool sorted == true
  _kd_tree->findNeighbors(result_set, query_pt.data(), nanoflann::SearchParams());

  return result_set;
}

void
PointLocatorNanoflann::kd_tree_radius_search(const Point & p, Real search_radius) const
{
  // Construct the search Point
  std::array<Real, LIBMESH_DIM> query_pt;
  for (int i=0; i<LIBMESH_DIM; ++i)
    query_pt[i] = p(i);

  // This type of search returns all KD-Tree Points within a given
  // search radius.  The output is stored in the "_ret_matches" data
  // structure. The search radius input to this function is not
  // squared, but Nanoflann is expecting squared distances, so we
  // square the input value.

  // std::size_t n_matches =
  _kd_tree->radiusSearch(query_pt.data(), search_radius * search_radius, _ret_matches, nanoflann::SearchParams());
}


const Elem *
PointLocatorNanoflann::operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator()", "PointLocatorNanoflann");

  std::size_t last_num_results = 0;
  std::size_t current_num_results = _initial_num_results;

  // Keep track of the number of elements checked in detail
  unsigned int n_elems_checked = 0;

  // Quick return if the Point is not in the BoundingBox of the local elements
  if (!_local_bbox.contains_point(p))
    goto done;

  // Do findNeighbors() search with expanding result set
//  while (current_num_results <= _max_num_results)
//    {
//      // Do the search
//      auto result_set = this->kd_tree_find_neighbors(p, current_num_results);
//
//      // Linear search over the list of candidate centroids starting from
//      // the index of the previous while loop, since we don't need to
//      // search those same centroids again.  Either returns from this
//      // function the Elem containing the searched for Point, or leaves
//      // this while loop.
//      for (std::size_t r = last_num_results; r < result_set.size(); ++r)
//        {
//          // Translate the Nanoflann index, which is from [0..n_centroids),
//          // into the corresponding Elem id from the mesh.
//          auto nanoflann_index = _ret_index[r];
//          auto elem_id = _ids[nanoflann_index];
//
//          // Debugging: print the results
//          // libMesh::out << "Centroid/Elem id = " << elem_id
//          //              << ", dist^2 = " << _out_dist_sqr[r]
//          //              << std::endl;
//
//          const Elem * candidate_elem = _mesh.elem_ptr(elem_id);
//
//          // Before we even check whether the candidate Elem actually
//          // contains the Point, we may need to check whether the
//          // candidate Elem is from an allowed subdomain.  If the
//          // candidate Elem is not from an allowed subdomain, we continue
//          // to the next one.
//          if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
//            {
//              // Debugging
//              // libMesh::out << "Elem " << elem_id << " was not from an allowed subdomain, continuing search." << std::endl;
//              continue;
//            }
//
//          // If we made it here, then the candidate Elem is from an
//          // allowed subdomain, so let's next check whether it contains
//          // the point. If the user set a custom tolerance, then we
//          // actually check close_to_point() rather than contains_point(),
//          // since this latter function warns about using non-default
//          // tolerances, but otherwise does the same test.
//          bool inside = _use_contains_point_tol ?
//            candidate_elem->close_to_point(p, _contains_point_tol) :
//            candidate_elem->contains_point(p);
//
//          // Increment the number of elements checked
//          n_elems_checked++;
//
//          // If the point is inside an Elem from an allowed subdomain, we are done.
//          if (inside)
//            {
//              // Debugging:
//              // libMesh::out << "Checked " << n_elems_checked << " nearby Elems before finding a containing Elem." << std::endl;
//
//              return candidate_elem;
//            }
//
//          // Debugging:
//          // libMesh::out << "Elem " << elem_id << " did not contain/was not close enough to Point " << p << std::endl;
//          // candidate_elem->print_info();
//        } // end for(r)
//
//      // Store the current number of results
//      last_num_results = current_num_results;
//
//      // Try the search again, requesting twice as many results as previously.
//      // current_num_results *= 2;
//
//      // Try the search again, requesting the max number of results.
//      current_num_results = _max_num_results;
//    } // end while


  // radiusSearch()
//  this->kd_tree_radius_search(p, _hmax/2);
//
//  // Debugging:
//  // std::cout << "Testing " << _ret_matches.size() << " elements within search radius." << std::endl;
//
//  for (const auto & pr : _ret_matches)
//    {
//      // Translate the Nanoflann index, which is from [0..n_centroids),
//      // into the corresponding Elem id from the mesh.
//      auto nanoflann_index = pr.first;
//      auto elem_id = _ids[nanoflann_index];
//
//      // Debugging: print the results
//      // libMesh::out << "Centroid/Elem id = " << elem_id
//      //              << ", dist^2 = " << pr.second
//      //              << std::endl;
//
//      const Elem * candidate_elem = _mesh.elem_ptr(elem_id);
//
//      // Before we even check whether the candidate Elem actually
//      // contains the Point, we may need to check whether the
//      // candidate Elem is from an allowed subdomain.  If the
//      // candidate Elem is not from an allowed subdomain, we continue
//      // to the next one.
//      if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
//        {
//          // Debugging
//          // libMesh::out << "Elem " << elem_id << " was not from an allowed subdomain, continuing search." << std::endl;
//          continue;
//        }
//
//      // If we made it here, then the candidate Elem is from an
//      // allowed subdomain, so let's next check whether it contains
//      // the point. If the user set a custom tolerance, then we
//      // actually check close_to_point() rather than contains_point(),
//      // since this latter function warns about using non-default
//      // tolerances, but otherwise does the same test.
//      bool inside = _use_contains_point_tol ?
//        candidate_elem->close_to_point(p, _contains_point_tol) :
//        candidate_elem->contains_point(p);
//
//      // Increment the number of elements checked
//      n_elems_checked++;
//
//      // If the point is inside an Elem from an allowed subdomain, we are done.
//      if (inside)
//        {
//          // Debugging:
//          // libMesh::out << "Checked " << n_elems_checked << " nearby Elems before finding a containing Elem." << std::endl;
//
//          return candidate_elem;
//        }
//    } // end for(r)


  // Hybrid approach: a "small" findNeighbors() search followed by an
  // exhaustive radiusSearch() based on the closest Elem's hmax() if
  // that fails.
  {
  auto result_set = this->kd_tree_find_neighbors(p, _initial_num_results);

  // We'll keep track of the largest "nearby" element's hmax in case
  // we need it later for a radiusSearch()
  Real largest_nearby_hmax = 0.;

  for (std::size_t r = last_num_results; r < result_set.size(); ++r)
    {
      // Translate the Nanoflann index, which is from [0..n_centroids),
      // into the corresponding Elem id from the mesh.
      auto nanoflann_index = _ret_index[r];
      auto elem_id = _ids[nanoflann_index];

      // Debugging: print the results
      // libMesh::out << "Centroid/Elem id = " << elem_id
      //              << ", dist^2 = " << _out_dist_sqr[r]
      //              << std::endl;

      const Elem * candidate_elem = _mesh.elem_ptr(elem_id);

      largest_nearby_hmax = std::max(candidate_elem->hmax(),
                                     largest_nearby_hmax);

      // Before we even check whether the candidate Elem actually
      // contains the Point, we may need to check whether the
      // candidate Elem is from an allowed subdomain.  If the
      // candidate Elem is not from an allowed subdomain, we continue
      // to the next one.
      if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
        {
          // Debugging
          // libMesh::out << "Elem " << elem_id << " was not from an allowed subdomain, continuing search." << std::endl;
          continue;
        }

      // If we made it here, then the candidate Elem is from an
      // allowed subdomain, so let's next check whether it contains
      // the point. If the user set a custom tolerance, then we
      // actually check close_to_point() rather than contains_point(),
      // since this latter function warns about using non-default
      // tolerances, but otherwise does the same test.
      bool inside = _use_contains_point_tol ?
        candidate_elem->close_to_point(p, _contains_point_tol) :
        candidate_elem->contains_point(p);

      // Increment the number of elements checked
      n_elems_checked++;

      // If the point is inside an Elem from an allowed subdomain, we are done.
      if (inside)
        {
          // Debugging:
          // libMesh::out << "Checked " << n_elems_checked << " nearby Elems before finding a containing Elem." << std::endl;

          return candidate_elem;
        }

      // Debugging:
      // libMesh::out << "Elem " << elem_id << " did not contain/was not close enough to Point " << p << std::endl;
      // candidate_elem->print_info();
    } // end for(r)

  // If we made it here without returning, try a more exhaustive
  // radiusSearch().

  // hmax of closest Elem
  // Real search_radius = _mesh.elem_ptr(_ids[_ret_index[0]])->hmax();

  // Some constant times the hmax of all "nearby" elements.
  // Using .5, 1, and 2 * local hmax = failed for my "challenging" test case
  // 10 passed but was too large

  // Let's try the _smaller_ of:
  // .) The global hmax() for the whole Mesh and
  // .) Some contant times the largest nearby hmax
  Real search_radius = std::min(_hmax, 6*largest_nearby_hmax);

  // Debugging:
  // libMesh::out << "Executing radiusSearch() with radius (not squared) = "
  //              << search_radius
  //              << std::endl;

  this->kd_tree_radius_search(p, search_radius);

  for (const auto & pr : _ret_matches)
    {
      // Translate the Nanoflann index, which is from [0..n_centroids),
      // into the corresponding Elem id from the mesh.
      auto nanoflann_index = pr.first;
      auto elem_id = _ids[nanoflann_index];

      // Debugging: print the results
      // libMesh::out << "Centroid/Elem id = " << elem_id
      //              << ", dist^2 = " << pr.second
      //              << std::endl;

      const Elem * candidate_elem = _mesh.elem_ptr(elem_id);

      // Before we even check whether the candidate Elem actually
      // contains the Point, we may need to check whether the
      // candidate Elem is from an allowed subdomain.  If the
      // candidate Elem is not from an allowed subdomain, we continue
      // to the next one.
      if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
        {
          // Debugging
          // libMesh::out << "Elem " << elem_id << " was not from an allowed subdomain, continuing search." << std::endl;
          continue;
        }

      // If we made it here, then the candidate Elem is from an
      // allowed subdomain, so let's next check whether it contains
      // the point. If the user set a custom tolerance, then we
      // actually check close_to_point() rather than contains_point(),
      // since this latter function warns about using non-default
      // tolerances, but otherwise does the same test.
      bool inside = _use_contains_point_tol ?
        candidate_elem->close_to_point(p, _contains_point_tol) :
        candidate_elem->contains_point(p);

      // Increment the number of elements checked
      n_elems_checked++;

      // If the point is inside an Elem from an allowed subdomain, we are done.
      if (inside)
        {
          // Debugging:
          // libMesh::out << "Checked " << n_elems_checked << " nearby Elems before finding a containing Elem." << std::endl;

          return candidate_elem;
        }
    } // end for(r)

  } // end scope for hybrid search approach

 done:
  // If we made it here, then at least one of the following happened:
  // .) The search Point was not in the BoundingBox of local Elems.
  // .) All the _max_num_results candidate elements were from non-allowed subdomains.
  // .) The Point was not inside _any_ of the _max_num_results candidate Elems.
  // Thus, if we are not in _out_of_mesh_mode, throw an error,
  // otherwise return nullptr to indicate that no suitable element was
  // found.
  if (!_out_of_mesh_mode)
    libmesh_error_msg("Point was not contained within the closest " << n_elems_checked <<
                      " elems (by centroid distance), and _out_of_mesh_mode was not enabled.");

  return nullptr;
}


void
PointLocatorNanoflann::operator() (const Point & p,
                                   std::set<const Elem *> & candidate_elements,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator() returning set", "PointLocatorNanoflann");

  std::size_t last_num_results = 0;
  std::size_t current_num_results = _initial_num_results;

  // Keep track of the number of elements checked in detail
  unsigned int n_elems_checked = 0;

  while (current_num_results < _max_num_results)
    {
      // Do the KD-Tree search
      auto result_set = this->kd_tree_find_neighbors(p, current_num_results);

      // Linear search over the list of candidate centroids starting from
      // the index of the previous while loop, since we don't need to
      // search those same centroids again.
      for (std::size_t r = last_num_results; r < result_set.size(); ++r)
        {
          // Translate the Nanoflann index, which is from [0..n_centroids),
          // into the corresponding Elem id from the mesh.
          auto nanoflann_index = _ret_index[r];
          auto elem_id = _ids[nanoflann_index];

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

          // Increment the number of elements checked
          n_elems_checked++;

          // If the point is contained in/close to an Elem from an
          // allowed subdomain, add it to the list.
          if (inside)
            candidate_elements.insert(candidate_elem);
        } // end for(r)

      // Possibly repeat the search, requesting twice as many results as previously.
      last_num_results = current_num_results;
      current_num_results *= 2;

      // If the candidate_elements set is non-empty, stop searching now.
      if (!candidate_elements.empty())
        break; // out of while loop
    } // end while

  // Debugging: for performance reasons, it may be useful to print the
  // number of Elems actually checked during the search.
  // libMesh::out << "Checked " << n_elems_checked << " nearby Elems before finding a containing Elem." << std::endl;
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



std::size_t
PointLocatorNanoflann::get_initial_num_results() const
{
  return _initial_num_results;
}



void
PointLocatorNanoflann::set_initial_num_results(std::size_t val)
{
  // Must request at least 1 result
  _initial_num_results = std::max(static_cast<std::size_t>(1), val);
}

std::size_t
PointLocatorNanoflann::get_max_num_results() const
{
  return _max_num_results;
}



void
PointLocatorNanoflann::set_max_num_results(std::size_t val)
{
  // Must be >= _initial_num_results
  _max_num_results = std::max(_initial_num_results, val);
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
