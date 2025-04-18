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



#ifndef LIBMESH_ENUM_XDR_MODE_H
#define LIBMESH_ENUM_XDR_MODE_H

namespace libMesh {

/**
 * Defines an \p enum for read/write mode in Xdr format.
 * \p READ, \p WRITE perform reading and writing in ASCII format,
 * and \p DECODE, \p ENCODE do the same in binary format.
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum XdrMODE : int;
 * reducing header file dependencies.
 */
enum XdrMODE : int
  {
    UNKNOWN = -1,
    ENCODE=0,
    DECODE,
    WRITE,
    READ
  };
}

#endif
