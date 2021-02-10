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



#ifndef LIBMESH_POINT_LOCATOR_NANOFLANN_H
#define LIBMESH_POINT_LOCATOR_NANOFLANN_H

#include "libmesh/libmesh_config.h"

// This class is not defined unless libmesh is built with Nanoflann support
#ifdef LIBMESH_HAVE_NANOFLANN

// libmesh includes
#include "libmesh/point_locator_base.h"
#include "libmesh/point.h"

// contrib includes
#include "libmesh/nanoflann.hpp"

// C++ includes
#include <vector>
#include <memory>
#include <unordered_map>

namespace libMesh
{

// Forward Declarations
class MeshBase;
class Point;
class Elem;

/**
 * This is a PointLocator that uses Nanoflann for its implementation.
 * Nanoflann is distributed with libmesh (see: contrib/nanoflann) and
 * libmesh must be built with nanoflann enabled for this class to
 * work.
 *
 * \author John W. Peterson
 * \date 2021
 */
class PointLocatorNanoflann : public PointLocatorBase
{
public:
  /**
   * Constructor. Needs the \p mesh in which the points should be
   * located. Optionally takes a pointer to a "master" PointLocator
   * object. If non-nullptr, this object simply forwards its calls
   * onto the master, so we can have multiple pointers that use the
   * same Nanoflann KD-Tree data structure.
   */
  PointLocatorNanoflann (const MeshBase & mesh,
                         const PointLocatorBase * master = nullptr);

  /**
   * Destructor.
   */
  virtual ~PointLocatorNanoflann ();

  /**
   * Restore to PointLocator to a just-constructed state.
   */
  virtual void clear() override;

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used. This function allocates dynamic memory with "new".
   */
  virtual void init() override;

  /**
   * Locates the element in which the point with global coordinates \p
   * p is located, optionally restricted to a set of allowed
   * subdomains.
   */
  virtual const Elem * operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains = nullptr) const override;

  /**
   * Locates a set of elements in proximity to the point with global
   * coordinates \p p. Optionally allows the user to restrict the
   * subdomains searched.  The idea here is that if a Point lies on
   * the boundary between two elements, so that it is "in" both
   * elements (to within some tolerance) then the candidate_elements
   * set will contain both of these elements instead of just picking
   * one or the other.
   */
  virtual void operator() (const Point & p,
                           std::set<const Elem *> & candidate_elements,
                           const std::set<subdomain_id_type> * allowed_subdomains = nullptr) const override;

  /**
   * Enables out-of-mesh mode. In this mode, if a searched-for Point
   * is not contained in any element of the Mesh, return nullptr
   * instead of throwing an error. By default, this mode is off.
   */
  virtual void enable_out_of_mesh_mode () override;

  /**
   * Disables out-of-mesh mode (default). See above.
   */
  virtual void disable_out_of_mesh_mode () override;

  /**
   * Set/get the "initial" number of results returned from each
   * Nanoflann search.  If a containing Elem (to within
   * _contains_point_tol) is not found after linear searching the
   * first _initial_num_results Elems returned by Nanoflann, then we
   * begin increasing this number until either a containing Elem is
   * found, or _max_num_results are searched.
   *
   * Remarks:
   * 1.) Although we do a linear search through the Nanoflann results,
   * the Nanoflann results are always sorted by distance to the
   * searched-for Point, so it's likely that the containing Elem will
   * be found at the beginning of the linear search rather than the
   * end. For nicely shaped elements, the containing Elem is often the
   * first or second one found in the Nanoflann results.
   *
   * 2.) I'm not sure about the relative cost of requesting more
   * results from the Nanoflann search than one actually needs, but
   * presumably reducing _initial_num_results will result in better
   * performance for your particular application up to a point. If, on
   * the other hand, _initial_num_results is too small, you will waste
   * time repeating Nanoflann searches with the same Point (again,
   * there might be some caching within Nanoflann itself that makes
   * repeated searches faster, I'm not sure).
   */
  std::size_t get_initial_num_results() const;
  void set_initial_num_results(std::size_t val);

  /**
   * Get/set the max number of results returned by the Nanoflann
   * search before giving up. Currently, the number of results
   * requested is doubled after each failed search until
   * _max_num_results is exceeded. Therefore, it should be both
   * reasonably robust and efficient to use a large value for
   * _max_num_results, with a hard maximum of the number of elements
   * in the Mesh.
   */
  std::size_t get_max_num_results() const;
  void set_max_num_results(std::size_t val);

  //
  // Required Nanoflann typedefs and APIs
  //

  /**
   * Floating point type used for storing coordinates
   */
  typedef Real coord_t;

  /**
   * Must return the number of data points
   */
  std::size_t kdtree_get_point_count() const;

  /**
   * Returns the squared distance between the vector (p1[0], p1[1], p1[2])
   * and the centroid of Elem "idx_p2" stored in _mesh
   */
  coord_t kdtree_distance(const coord_t * p1,
                          const std::size_t idx_p2,
                          std::size_t size) const;

  /**
   * Returns the dim'th component of the idx'th centroid.
   */
  coord_t kdtree_get_pt(const std::size_t idx, int dim) const;

  /**
   * Optional bounding-box computation: return false to default to a
   * standard bbox computation loop.  Return true if the BBOX can
   * be computed more efficiently (and returned in "bb") than the
   * standard bounding box computation. The BBOX template parameter
   * must at least support bb.size() to find out the expected
   * dimensionality of the box (e.g. 2 or 3 for point clouds).
   */
  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /*bb*/) const { return false; }

protected:

  /**
   * \p true if out-of-mesh mode is enabled.  See \p
   * enable_out_of_mesh_mode() for details.
   */
  bool _out_of_mesh_mode;

  /**
   * Lists of Mesh Elem ids and centroids which is created when init()
   * is called. These are the points used in the KD-Tree. We keep two
   * separate vectors since the local, active Elems are not numbered
   * contiguously in general.
   */
  std::vector<dof_id_type> _ids;
  std::vector<Point> _centroids;

  /**
   * Defaults to 30. This is the number of results initially returned by Nanoflann
   * when operator() is called. If a containing Elem is not found with the first
   * _initial_num_results, the KD-Tree search is repeated with twice as many
   * results requested, until _max_num_results is reached.
   */
  std::size_t _initial_num_results;

  /**
   * Defaults to 1000. We will continue to request more and more
   * results from KD-Tree searches until we either reach this number,
   * or find a containg Elem.
   */
  std::size_t _max_num_results;

  // kd_tree will be initialized during init() and then automatically
  // cleaned up by the destructor. We always create a LIBMESH_DIM
  // dimensional Nanoflann object.
  typedef nanoflann::L2_Simple_Adaptor<Real, PointLocatorNanoflann> adapter_t;
  typedef nanoflann::KDTreeSingleIndexAdaptor<adapter_t, PointLocatorNanoflann, LIBMESH_DIM> kd_tree_t;
  std::unique_ptr<kd_tree_t> _kd_tree;

  /**
   * Helper function that wraps the call to the KDTree's
   * findNeighbors() routine.  Must be passed the Point to search for
   * and the number of results to return. Stores the results of the
   * search in the _ret_index and _out_dist_sqr class members and
   * returns a KNNResultSet, which has pointers to the index and
   * distance data.
   */
  nanoflann::KNNResultSet<Real>
  kd_tree_find_neighbors(const Point & p,
                         std::size_t num_results) const;

  /**
   * The operator() functions on PointLocator-derived classes are
   * const, so to avoid re-allocating these result data structures every
   * time operator() is called, they have to be mutable.
   */
  mutable std::vector<std::size_t> _ret_index;
  mutable std::vector<Real> _out_dist_sqr;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_NANOFLANN
#endif // LIBMESH_POINT_LOCATOR_NANOFLANN_H
