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
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes cont'd
#include "libmesh/cell_inf_hex8.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfHex8 class static member initializations
const int InfHex8::num_nodes;
const int InfHex8::nodes_per_side;
const int InfHex8::nodes_per_edge;

const unsigned int InfHex8::side_nodes_map[InfHex8::num_sides][InfHex8::nodes_per_side] =
  {
    { 0, 1, 2, 3}, // Side 0
    { 0, 1, 4, 5}, // Side 1
    { 1, 2, 5, 6}, // Side 2
    { 2, 3, 6, 7}, // Side 3
    { 3, 0, 7, 4}  // Side 4
  };

const unsigned int InfHex8::edge_nodes_map[InfHex8::num_edges][InfHex8::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {0, 4}, // Edge 4
    {1, 5}, // Edge 5
    {2, 6}, // Edge 6
    {3, 7}  // Edge 7
  };

// ------------------------------------------------------------
// InfHex8 class member functions

bool InfHex8::is_node_on_side(const unsigned int n,
                              const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
InfHex8::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
InfHex8::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool InfHex8::is_node_on_edge(const unsigned int n,
                              const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



Order InfHex8::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> InfHex8::build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (i)
    {
    case 0: // the base face
      {
        face = std::make_unique<Quad4>();
        break;
      }

      // connecting to another infinite element
    case 1:
    case 2:
    case 3:
    case 4:
      {
        face = std::make_unique<InfQuad4>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(InfHex8::side_nodes_map[i][n]));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
}


void InfHex8::build_side_ptr (std::unique_ptr<Elem> & side,
                              const unsigned int i)
{
  this->side_ptr(side, i);
  side->set_interior_parent(this);
  side->inherit_data_from(*this);
}



std::unique_ptr<Elem> InfHex8::build_edge_ptr (const unsigned int i)
{
  if (i < 4) // base edges
    return this->simple_build_edge_ptr<Edge2,InfHex8>(i);

  // infinite edges
  return this->simple_build_edge_ptr<InfEdge2,InfHex8>(i);
}



void InfHex8::build_edge_ptr (std::unique_ptr<Elem> & edge,
                              const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  switch (i)
    {
      // the base edges
    case 0:
    case 1:
    case 2:
    case 3:
      {
        if (!edge.get() || edge->type() != EDGE2)
          {
            edge = this->build_edge_ptr(i);
            return;
          }
        break;
      }

      // the infinite edges
    case 4:
    case 5:
    case 6:
    case 7:
      {
        if (!edge.get() || edge->type() != INFEDGE2)
          {
            edge = this->build_edge_ptr(i);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid edge i = " << i);
    }

  edge->inherit_data_from(*this);

  // Set the nodes
  for (auto n : edge->node_index_range())
    edge->set_node(n, this->node_ptr(InfHex8::edge_nodes_map[i][n]));
}



void InfHex8::connectivity(const unsigned int libmesh_dbg_var(sc),
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
        conn[5] = this->node_id(5)+1;
        conn[6] = this->node_id(6)+1;
        conn[7] = this->node_id(7)+1;
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



#ifdef LIBMESH_ENABLE_AMR

const Real InfHex8::_embedding_matrix[InfHex8::num_children][InfHex8::num_nodes][InfHex8::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //     0      1      2      3      4      5      6      7 th parent N.(ode)
      {    1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
      {    0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
      {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 2
      {    0.5,   0.0,   0.0,   0.5,   0.0,   0.0,   0.0,   0.0}, // 3
      {    0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0}, // 4
      {    0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0,   0.0}, // 5
      {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 6
      {    0.0,   0.0,   0.0,   0.0,   0.5,   0.0,   0.0,   0.5}  // 7
    },

    // embedding matrix for child 1
    {
      //     0      1      2      3      4      5      6      7 th parent N.(ode)
      {    0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
      {    0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
      {    0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0}, // 2
      {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 3
      {    0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0,   0.0}, // 4
      {    0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0}, // 5
      {    0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0}, // 6
      {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}  // 7
    },

    // embedding matrix for child 2
    {
      //     0      1      2      3      4      5      6      7 th parent N.(ode)
      {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
      {    0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
      {    0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 2
      {    0.0,   0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0}, // 3
      {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 4
      {    0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0}, // 5
      {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0}, // 6
      {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5}  // 7
    },

    // embedding matrix for child 3
    {
      //     0      1      2      3      4      5      6      7 th parent N.(ode)
      {    0.5,   0.0,   0.0,   0.5,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
      {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 1
      {    0.0,   0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0}, // 2
      {    0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0}, // 3
      {    0.0,   0.0,   0.0,   0.0,   0.5,   0.0,   0.0,   0.5}, // 4
      {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 5
      {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5}, // 6
      {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0}  // 7
    }
  };


#endif


void
InfHex8::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4nodes(4,5,6,7);
      swap4neighbors(1,2,3,4);
    }
}


void
InfHex8::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2nodes(4,5);
  swap2nodes(6,7);
  swap2neighbors(0,4);
  swap2boundarysides(0,4,boundary_info);
  swap2boundaryedges(1,3,boundary_info);
  swap2boundaryedges(4,5,boundary_info);
  swap2boundaryedges(6,7,boundary_info);
}


ElemType
InfHex8::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s == 0)
    return QUAD4;
  return INFQUAD4;
}


} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
