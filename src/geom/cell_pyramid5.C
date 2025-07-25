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


// Local includes
#include "libmesh/cell_pyramid5.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/cell_hex8.h"

namespace libMesh
{




// ------------------------------------------------------------
// Pyramid5 class static member initializations
const int Pyramid5::num_nodes;
const int Pyramid5::nodes_per_side;
const int Pyramid5::nodes_per_edge;

const unsigned int Pyramid5::side_nodes_map[Pyramid5::num_sides][Pyramid5::nodes_per_side] =
  {
    {0, 1, 4, 99}, // Side 0
    {1, 2, 4, 99}, // Side 1
    {2, 3, 4, 99}, // Side 2
    {3, 0, 4, 99}, // Side 3
    {0, 3, 2,  1}  // Side 4
  };

const unsigned int Pyramid5::edge_nodes_map[Pyramid5::num_edges][Pyramid5::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {0, 4}, // Edge 4
    {1, 4}, // Edge 5
    {2, 4}, // Edge 6
    {3, 4}  // Edge 7
  };

// ------------------------------------------------------------
// Pyramid5 class member functions

bool Pyramid5::is_vertex(const unsigned int libmesh_dbg_var(n)) const
{
  libmesh_assert_not_equal_to (n, invalid_uint);
  return true;
}

bool Pyramid5::is_edge(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_face(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Pyramid5::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s == 4) ? 0 : 1;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

std::vector<unsigned>
Pyramid5::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Pyramid5::is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



bool Pyramid5::has_affine_map() const
{
  //  Point v = this->point(3) - this->point(0);
  //  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
  return false;
}



Order Pyramid5::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> Pyramid5::build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  switch (i)
    {
    case 0: // triangular face 1
    case 1: // triangular face 2
    case 2: // triangular face 3
    case 3: // triangular face 4
      {
        face = std::make_unique<Tri3>();
        break;
      }
    case 4: // the quad face at z=0
      {
        face = std::make_unique<Quad4>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(Pyramid5::side_nodes_map[i][n]));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
}



void Pyramid5::build_side_ptr (std::unique_ptr<Elem> & side,
                               const unsigned int i)
{
  this->side_ptr(side, i);
  side->set_interior_parent(this);
  side->inherit_data_from(*this);
}



std::unique_ptr<Elem> Pyramid5::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge2,Pyramid5>(i);
}



void Pyramid5::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Pyramid5>(edge, i, EDGE2);
}



void Pyramid5::connectivity(const unsigned int libmesh_dbg_var(sc),
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(8);
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(3)+1;
        conn[4] = this->node_id(4)+1;
        conn[5] = this->node_id(4)+1;
        conn[6] = this->node_id(4)+1;
        conn[7] = this->node_id(4)+1;
        return;
      }

    case VTK:
      {
        conn.resize(5);
        conn[0] = this->node_id(3);
        conn[1] = this->node_id(2);
        conn[2] = this->node_id(1);
        conn[3] = this->node_id(0);
        conn[4] = this->node_id(4);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}


Point Pyramid5::true_centroid () const
{
  // Call Hex8 static helper function, passing 4 copies of the final
  // vertex point, effectively treating the Pyramid as a degenerate
  // hexahedron.  In my testing, this still seems to give correct
  // results.
  return Hex8::centroid_from_points(
    point(0), point(1), point(2), point(3),
    point(4), point(4), point(4), point(4));
}



Real Pyramid5::volume () const
{
  // The pyramid with a bilinear base has volume given by the
  // formula in: "Calculation of the Volume of a General Hexahedron
  // for Flow Predictions", AIAA Journal v.23, no.6, 1984, p.954-
  Point
    x0 = point(0), x1 = point(1), x2 = point(2),
    x3 = point(3), x4 = point(4);

  // Construct various edge and diagonal vectors.
  Point v40 = x0 - x4;
  Point v13 = x3 - x1;
  Point v02 = x2 - x0;
  Point v03 = x3 - x0;
  Point v01 = x1 - x0;

  // Finally, ready to return the volume!
  return
    triple_product(v40, v13, v02) / 6. +
    triple_product(v02, v01, v03) / 12.;
}

BoundingBox
Pyramid5::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}

void Pyramid5::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4neighbors(0,1,2,3);
    }
}


void Pyramid5::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2neighbors(1,3);
  swap2boundarysides(1,3,boundary_info);
  swap2boundaryedges(1,3,boundary_info);
  swap2boundaryedges(4,5,boundary_info);
  swap2boundaryedges(6,7,boundary_info);
}


ElemType Pyramid5::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s < 4)
    return TRI3;
  return QUAD4;
}


} // namespace libMesh
