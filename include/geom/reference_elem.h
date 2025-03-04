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



#ifndef LIBMESH_REFERENCE_ELEM_H
#define LIBMESH_REFERENCE_ELEM_H

// Local includes
#include "libmesh/libmesh_common.h"

namespace libMesh
{

// forward declarations
class Elem;
enum ElemType : int;

/**
 * This namespace implements singleton reference elements for each
 * fundamental element type supported by \p libMesh.
 *
 * \author Benjamin S. Kirk
 * \date 2013
 * \brief Namespace providing access to reference geometric element types.
 */
namespace ReferenceElem
{
/**
 * \returns A constant reference to the reference element of
 * the user-requested type.
 */
const Elem & get (const ElemType type_in);

} // namespace ReferenceElem


} // namespace libMesh


#endif // LIBMESH_REFERENCE_ELEM_H
