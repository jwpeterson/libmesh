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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA


// Local includes
#include "libmesh/number_lookups.h"

namespace libMesh
{

// These numbers need to go up to at least maximum_totalorder - 2

// triangular_number_*: indices for triangle interiors
const unsigned char triangular_number_row[] = {
  0,
  1, 1,
  2, 2, 2,
  3, 3, 3, 3,
  4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9
};

const unsigned char triangular_number_column[] = {
  0,
  0, 1,
  0, 1, 2,
  0, 1, 2, 3,
  0, 1, 2, 3, 4,
  0, 1, 2, 3, 4, 5,
  0, 1, 2, 3, 4, 5, 6,
  0, 1, 2, 3, 4, 5, 6, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9
};


// square_number_*: indices for square interiors, cube faces
const unsigned char square_number_column[] = {
  0,
  0, 1, 1,
  0, 1, 2, 2, 2,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9
};

const unsigned char square_number_row[] = {
  0,
  1, 1, 0,
  2, 2, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
};


// cube_number_*: indices for cube interiors
const unsigned char cube_number_column[] = {
  0,

  0, 1, 1,
  0, 1, 1,
  0,

  0, 1, 2, 2, 2,
  0, 1, 2, 2, 2,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0,

  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0,

  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0,

  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0,

  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0,

  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0,

  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0,

  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7,
  0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6,
  0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5,
  0, 1, 2, 3, 4, 4, 4, 4, 4,
  0, 1, 2, 3, 3, 3, 3,
  0, 1, 2, 2, 2,
  0, 1, 1,
  0
};

const unsigned char cube_number_row[] = {
  0,

  1, 1, 0,
  1, 1, 0,
  0,

  2, 2, 2, 1, 0,
  2, 2, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0,

  3, 3, 3, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0,

  4, 4, 4, 4, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0,

  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0,

  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0,

  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0,

  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0,

  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  7, 7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0,
  6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0,
  5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0,
  4, 4, 4, 4, 4, 3, 2, 1, 0,
  3, 3, 3, 3, 2, 1, 0,
  2, 2, 2, 1, 0,
  1, 1, 0,
  0
};

const unsigned char cube_number_page[] = {
  0,

  0, 0, 0,
  1, 1, 1,
  1,

  0, 0, 0, 0, 0,
  1, 1, 1, 1, 1,
  2, 2, 2, 2, 2,
  2, 2, 2,
  2,

  0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3,
  3, 3, 3,
  3,

  0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4,
  4, 4, 4,
  4,

  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5,
  5, 5, 5,
  5,

  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6,
  6, 6, 6,
  6,

  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7,
  7, 7, 7,
  7,

  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8,
  8, 8, 8,
  8,

  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9,
  9, 9, 9,
  9
};

// tetrahedral_number_*: indices for tetrahedron interiors
const unsigned char tetrahedral_number_column[] = {
  0,

  0, 1,
  0,

  0, 1, 2,
  0, 1,
  0,

  0, 1, 2, 3,
  0, 1, 2,
  0, 1,
  0,

  0, 1, 2, 3, 4,
  0, 1, 2, 3,
  0, 1, 2,
  0, 1,
  0,

  0, 1, 2, 3, 4, 5,
  0, 1, 2, 3, 4,
  0, 1, 2, 3,
  0, 1, 2,
  0, 1,
  0,

  0, 1, 2, 3, 4, 5, 6,
  0, 1, 2, 3, 4, 5,
  0, 1, 2, 3, 4,
  0, 1, 2, 3,
  0, 1, 2,
  0, 1,
  0,

  0, 1, 2, 3, 4, 5, 6, 7,
  0, 1, 2, 3, 4, 5, 6,
  0, 1, 2, 3, 4, 5,
  0, 1, 2, 3, 4,
  0, 1, 2, 3,
  0, 1, 2,
  0, 1,
  0,

  0, 1, 2, 3, 4, 5, 6, 7, 8,
  0, 1, 2, 3, 4, 5, 6, 7,
  0, 1, 2, 3, 4, 5, 6,
  0, 1, 2, 3, 4, 5,
  0, 1, 2, 3, 4,
  0, 1, 2, 3,
  0, 1, 2,
  0, 1,
  0,

  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
  0, 1, 2, 3, 4, 5, 6, 7, 8,
  0, 1, 2, 3, 4, 5, 6, 7,
  0, 1, 2, 3, 4, 5, 6,
  0, 1, 2, 3, 4, 5,
  0, 1, 2, 3, 4,
  0, 1, 2, 3,
  0, 1, 2,
  0, 1,
  0
};

const unsigned char tetrahedral_number_row[] = {
  0,

  1, 1,
  1,

  2, 2, 2,
  2, 2,
  2,

  3, 3, 3, 3,
  3, 3, 3,
  3, 3,
  3,

  4, 4, 4, 4, 4,
  4, 4, 4, 4,
  4, 4, 4,
  4, 4,
  4,

  5, 5, 5, 5, 5,
  5, 5, 5, 5,
  5, 5, 5,
  5, 5,
  5,

  6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6,
  6, 6, 6, 6,
  6, 6, 6,
  6, 6,
  6,

  7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7,
  7, 7, 7, 7,
  7, 7, 7,
  7, 7,
  7,

  8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8,
  8, 8, 8, 8,
  8, 8, 8,
  8, 8,
  8,

  9, 9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9,
  9, 9, 9, 9,
  9, 9, 9,
  9, 9,
  9
};

const unsigned char tetrahedral_number_page[] = {
  0,

  0, 0,
  1,

  0, 0, 0,
  1, 1,
  2,

  0, 0, 0, 0,
  1, 1, 1,
  2, 2,
  3,

  0, 0, 0, 0, 0,
  1, 1, 1, 1,
  2, 2, 2,
  3, 3,
  4,

  0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1,
  2, 2, 2, 2,
  3, 3, 3,
  4, 4,
  5,

  0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2,
  3, 3, 3, 3,
  4, 4, 4,
  5, 5,
  6,

  0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3,
  4, 4, 4, 4,
  5, 5, 5,
  6, 6,
  7,

  0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4,
  5, 5, 5, 5,
  6, 6, 6,
  7, 7,
  8,

  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5,
  6, 6, 6, 6,
  7, 7, 7,
  8, 8,
  9
};

// prism_number_*: indices for prism interiors
//
// These are split into a linear and a triangular index, rather than
// into three linear indices
const unsigned char prism_number_triangle[] = {

  0,
  0,

  1, 2,
  1, 2,
  0,
  1, 2,

  3, 4, 5,
  3, 4, 5,
  3, 4, 5,
  0,
  1, 2,
  3, 4, 5,

  6, 7, 8, 9,
  6, 7, 8, 9,
  6, 7, 8, 9,
  6, 7, 8, 9,
  0,
  1, 2,
  3, 4, 5,
  6, 7, 8, 9,

  10, 11, 12, 13, 14,
  10, 11, 12, 13, 14,
  10, 11, 12, 13, 14,
  10, 11, 12, 13, 14,
  10, 11, 12, 13, 14,
  0,
  1, 2,
  3, 4, 5,
  6, 7, 8, 9,
  10, 11, 12, 13, 14,

  15, 16, 17, 18, 19, 20,
  15, 16, 17, 18, 19, 20,
  15, 16, 17, 18, 19, 20,
  15, 16, 17, 18, 19, 20,
  15, 16, 17, 18, 19, 20,
  15, 16, 17, 18, 19, 20,
  15, 16, 17, 18, 19, 20,
  0,
  1, 2,
  3, 4, 5,
  6, 7, 8, 9,
  10, 11, 12, 13, 14,
  15, 16, 17, 18, 19, 20,

  21, 22, 23, 24, 25, 26, 27,
  21, 22, 23, 24, 25, 26, 27,
  21, 22, 23, 24, 25, 26, 27,
  21, 22, 23, 24, 25, 26, 27,
  21, 22, 23, 24, 25, 26, 27,
  21, 22, 23, 24, 25, 26, 27,
  21, 22, 23, 24, 25, 26, 27,
  21, 22, 23, 24, 25, 26, 27,
  0,
  1, 2,
  3, 4, 5,
  6, 7, 8, 9,
  10, 11, 12, 13, 14,
  15, 16, 17, 18, 19, 20,
  21, 22, 23, 24, 25, 26, 27,

  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  28, 29, 30, 31, 32, 33, 34, 35,
  0,
  1, 2,
  3, 4, 5,
  6, 7, 8, 9,
  10, 11, 12, 13, 14,
  15, 16, 17, 18, 19, 20,
  21, 22, 23, 24, 25, 26, 27,
  28, 29, 30, 31, 32, 33, 34, 35,

  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  0,
  1, 2,
  3, 4, 5,
  6, 7, 8, 9,
  10, 11, 12, 13, 14,
  15, 16, 17, 18, 19, 20,
  21, 22, 23, 24, 25, 26, 27,
  28, 29, 30, 31, 32, 33, 34, 35,
  36, 37, 38, 39, 40, 41, 42, 43, 44
};

const unsigned char prism_number_page[] = {
  0,
  1,

  0, 0,
  1, 1,
  2,
  2, 2,

  0, 0, 0,
  1, 1, 1,
  2, 2, 2,
  3,
  3, 3,
  3, 3, 3,

  0, 0, 0, 0,
  1, 1, 1, 1,
  2, 2, 2, 2,
  3, 3, 3, 3,
  4,
  4, 4,
  4, 4, 4,
  4, 4, 4, 4,

  0, 0, 0, 0, 0,
  1, 1, 1, 1, 1,
  2, 2, 2, 2, 2,
  3, 3, 3, 3, 3,
  4, 4, 4, 4, 4,
  5,
  5, 5,
  5, 5, 5,
  5, 5, 5, 5,
  5, 5, 5, 5, 5,

  0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5,
  6,
  6, 6,
  6, 6, 6,
  6, 6, 6, 6,
  6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6,

  0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6,
  7,
  7, 7,
  7, 7, 7,
  7, 7, 7, 7,
  7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7,

  0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7,
  8,
  8, 8,
  8, 8, 8,
  8, 8, 8, 8,
  8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8,

  0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8,
  9,
  9, 9,
  9, 9, 9,
  9, 9, 9, 9,
  9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9,
  9, 9, 9, 9, 9, 9, 9, 9, 9
};

} // namespace libMesh
