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



#ifndef LIBMESH_ENUM_PARTITIONER_TYPE_H
#define LIBMESH_ENUM_PARTITIONER_TYPE_H

namespace libMesh {

/**
 * Defines an \p enum for mesh partitioner types
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum EigenSolverType : int;
 * reducing header file dependencies.
 */
enum PartitionerType : int {
                      CENTROID_PARTITIONER=0,
                      LINEAR_PARTITIONER,
                      SFC_PARTITIONER,
                      HILBERT_SFC_PARTITIONER,
                      MORTON_SFC_PARTITIONER,
                      METIS_PARTITIONER,
                      PARMETIS_PARTITIONER,
                      SUBDOMAIN_PARTITIONER,
                      MAPPED_SUBDOMAIN_PARTITIONER,
                      // Invalid
                      INVALID_PARTITIONER};

} // namespace libMesh

#endif
