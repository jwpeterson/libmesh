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
 * \date 2020
 */
class PointLocatorNanoflann : public PointLocatorBase
{
public:
  /**
   * Constructor. Needs the \p mesh in which the points should be
   * located. Optionally takes a pointer to a "primary" PointLocator
   * object. If non-nullptr, this object simply forwards its calls
   * onto the primary, so we can have multiple pointers that use the
   * same Nanoflann KD-Tree data structure.
   */
  PointLocatorNanoflann (const MeshBase & mesh,
                         const PointLocatorBase * primary = nullptr);

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

protected:

  /**
   * \p true if out-of-mesh mode is enabled.  See \p
   * enable_out_of_mesh_mode() for details.
   */
  bool _out_of_mesh_mode;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_NANOFLANN
#endif // LIBMESH_POINT_LOCATOR_NANOFLANN_H
