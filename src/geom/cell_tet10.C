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
#include "libmesh/cell_tet10.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tet10 class static member initializations
const int Tet10::num_nodes;
const int Tet10::nodes_per_side;
const int Tet10::nodes_per_edge;

const unsigned int Tet10::side_nodes_map[Tet10::num_sides][Tet10::nodes_per_side] =
  {
    {0, 2, 1, 6, 5, 4}, // Side 0
    {0, 1, 3, 4, 8, 7}, // Side 1
    {1, 2, 3, 5, 9, 8}, // Side 2
    {2, 0, 3, 6, 7, 9}  // Side 3
  };

const unsigned int Tet10::edge_nodes_map[Tet10::num_edges][Tet10::nodes_per_edge] =
  {
    {0, 1, 4}, // Edge 0
    {1, 2, 5}, // Edge 1
    {0, 2, 6}, // Edge 2
    {0, 3, 7}, // Edge 3
    {1, 3, 8}, // Edge 4
    {2, 3, 9}  // Edge 5
  };

// ------------------------------------------------------------
// Tet10 class member functions

bool Tet10::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

bool Tet10::is_edge(const unsigned int i) const
{
  if (i < 4)
    return false;
  return true;
}

bool Tet10::is_face(const unsigned int) const
{
  return false;
}

bool Tet10::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Tet10::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Tet10::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Tet10::is_node_on_edge(const unsigned int n,
                            const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}


#ifdef LIBMESH_ENABLE_AMR

// This function only works if LIBMESH_ENABLE_AMR...
bool Tet10::is_child_on_side(const unsigned int c,
                             const unsigned int s) const
{
  // Table of local IDs for the midege nodes on the side opposite a given node.
  // See the ASCII art in the header file for this class to confirm this.
  const unsigned int midedge_nodes_opposite[4][3] =
    {
      {5,8,9}, // midedge nodes opposite node 0
      {6,7,9}, // midedge nodes opposite node 1
      {4,7,8}, // midedge nodes opposite node 2
      {4,5,6}  // midedge nodes opposite node 3
    };

  // Call the base class helper function
  return Tet::is_child_on_side_helper(c, s, midedge_nodes_opposite);
}

#else

bool Tet10::is_child_on_side(const unsigned int /*c*/,
                             const unsigned int /*s*/) const
{
  libmesh_not_implemented();
  return false;
}

#endif //LIBMESH_ENABLE_AMR



bool Tet10::has_affine_map() const
{
  // Make sure edges are straight
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals
      ((this->point(4) - this->point(0))*2, affine_tol))
    return false;
  v = this->point(2) - this->point(1);
  if (!v.relative_fuzzy_equals
      ((this->point(5) - this->point(1))*2, affine_tol))
    return false;
  v = this->point(2) - this->point(0);
  if (!v.relative_fuzzy_equals
      ((this->point(6) - this->point(0))*2, affine_tol))
    return false;
  v = this->point(3) - this->point(0);
  if (!v.relative_fuzzy_equals
      ((this->point(7) - this->point(0))*2, affine_tol))
    return false;
  v = this->point(3) - this->point(1);
  if (!v.relative_fuzzy_equals
      ((this->point(8) - this->point(1))*2, affine_tol))
    return false;
  v = this->point(3) - this->point(2);
  if (!v.relative_fuzzy_equals
      ((this->point(9) - this->point(2))*2, affine_tol))
    return false;
  return true;
}



Order Tet10::default_order() const
{
  return SECOND;
}



unsigned int Tet10::local_side_node(unsigned int side,
                                    unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Tet10::nodes_per_side);

  return Tet10::side_nodes_map[side][side_node];
}



unsigned int Tet10::local_edge_node(unsigned int edge,
                                    unsigned int edge_node) const
{
  libmesh_assert_less (edge, this->n_edges());
  libmesh_assert_less (edge_node, Tet10::nodes_per_edge);

  return Tet10::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> Tet10::build_side_ptr (const unsigned int i)
{
  return this->simple_build_side_ptr<Tri6, Tet10>(i);
}



void Tet10::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  this->simple_build_side_ptr<Tet10>(side, i, TRI6);
}



std::unique_ptr<Elem> Tet10::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge3,Tet10>(i);
}



void Tet10::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Tet10>(edge, i, EDGE3);
}



void Tet10::connectivity(const unsigned int sc,
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
        switch (sc)
          {


            // Linear sub-tet 0
          case 0:

            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(4)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(7)+1;
            conn[5] = this->node_id(7)+1;
            conn[6] = this->node_id(7)+1;
            conn[7] = this->node_id(7)+1;

            return;

            // Linear sub-tet 1
          case 1:

            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(5)+1;
            conn[3] = this->node_id(5)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;

            // Linear sub-tet 2
          case 2:

            conn[0] = this->node_id(5)+1;
            conn[1] = this->node_id(2)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(9)+1;
            conn[5] = this->node_id(9)+1;
            conn[6] = this->node_id(9)+1;
            conn[7] = this->node_id(9)+1;

            return;

            // Linear sub-tet 3
          case 3:

            conn[0] = this->node_id(7)+1;
            conn[1] = this->node_id(8)+1;
            conn[2] = this->node_id(9)+1;
            conn[3] = this->node_id(9)+1;
            conn[4] = this->node_id(3)+1;
            conn[5] = this->node_id(3)+1;
            conn[6] = this->node_id(3)+1;
            conn[7] = this->node_id(3)+1;

            return;

            // Linear sub-tet 4
          case 4:

            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(8)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(7)+1;
            conn[5] = this->node_id(7)+1;
            conn[6] = this->node_id(7)+1;
            conn[7] = this->node_id(7)+1;

            return;

            // Linear sub-tet 5
          case 5:

            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(5)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;

            // Linear sub-tet 6
          case 6:

            conn[0] = this->node_id(5)+1;
            conn[1] = this->node_id(9)+1;
            conn[2] = this->node_id(6)+1;
            conn[3] = this->node_id(6)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;

            // Linear sub-tet 7
          case 7:

            conn[0] = this->node_id(7)+1;
            conn[1] = this->node_id(6)+1;
            conn[2] = this->node_id(9)+1;
            conn[3] = this->node_id(9)+1;
            conn[4] = this->node_id(8)+1;
            conn[5] = this->node_id(8)+1;
            conn[6] = this->node_id(8)+1;
            conn[7] = this->node_id(8)+1;

            return;


          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    case VTK:
      {
        // VTK connectivity for VTK_QUADRATIC_TETRA matches libMesh's own.
        conn.resize(Tet10::num_nodes);
        for (auto i : index_range(conn))
          conn[i] = this->node_id(i);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



const unsigned short int Tet10::_second_order_vertex_child_number[10] =
  {
    99,99,99,99, // Vertices
    0,1,0,0,1,2  // Edges
  };



const unsigned short int Tet10::_second_order_vertex_child_index[10] =
  {
    99,99,99,99, // Vertices
    1,2,2,3,3,3  // Edges
  };



std::pair<unsigned short int, unsigned short int>
Tet10::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



unsigned short int Tet10::second_order_adjacent_vertex (const unsigned int n,
                                                        const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



const unsigned short int Tet10::_second_order_adjacent_vertices[6][2] =
  {
    {0, 1}, // vertices adjacent to node 4
    {1, 2}, // vertices adjacent to node 5
    {0, 2}, // vertices adjacent to node 6
    {0, 3}, // vertices adjacent to node 7
    {1, 3}, // vertices adjacent to node 8
    {2, 3}  // vertices adjacent to node 9
  };





#ifdef LIBMESH_ENABLE_AMR

const Real Tet10::_embedding_matrix[Tet10::num_children][Tet10::num_nodes][Tet10::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
      { 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 4
      {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 5
      { 0.375,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.,    0.}, // 6
      { 0.375,    0.,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.}, // 7
      {    0.,-0.125,    0.,-0.125,   0.5,    0.,    0.,   0.5,  0.25,    0.}, // 8
      {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}  // 9
    },

    // embedding matrix for child 1
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
      {-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 4
      {    0., 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 5
      {-0.125,    0.,-0.125,    0.,   0.5,   0.5,  0.25,    0.,    0.,    0.}, // 6
      {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 7
      {    0., 0.375,    0.,-0.125,    0.,    0.,    0.,    0.,  0.75,    0.}, // 8
      {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}  // 9
    },

    // embedding matrix for child 2
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 3
      {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 4
      {    0.,-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 5
      {-0.125,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.,    0.}, // 6
      {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 7
      {    0.,-0.125,    0.,-0.125,    0.,   0.5,    0.,    0.,  0.25,   0.5}, // 8
      {    0.,    0., 0.375,-0.125,    0.,    0.,    0.,    0.,    0.,  0.75}  // 9
    },

    // embedding matrix for child 3
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
      {    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.}, // 3
      {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 4
      {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}, // 5
      {-0.125,    0.,-0.125,    0.,    0.,    0.,  0.25,   0.5,    0.,   0.5}, // 6
      {-0.125,    0.,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.}, // 7
      {    0.,-0.125,    0., 0.375,    0.,    0.,    0.,    0.,  0.75,    0.}, // 8
      {    0.,    0.,-0.125, 0.375,    0.,    0.,    0.,    0.,    0.,  0.75}  // 9
    },

    // embedding matrix for child 4
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
      {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 4
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 5
      {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 6
      {    0.,-0.125,    0.,-0.125,   0.5,    0.,    0.,   0.5,  0.25,    0.}, // 7
      {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 8
      {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}  // 9
    },

    // embedding matrix for child 5
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
      {-0.125,    0.,-0.125,    0.,   0.5,   0.5,  0.25,    0.,    0.,    0.}, // 4
      {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 5
      {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 6
      {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 7
      {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}, // 8
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}  // 9
    },

    // embedding matrix for child 6
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
      {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 4
      {    0.,-0.125,    0.,-0.125,    0.,   0.5,    0.,    0.,  0.25,   0.5}, // 5
      {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 6
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 7
      {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}, // 8
      {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}  // 9
    },

    // embedding matrix for child 7
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 4
      {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}, // 5
      {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 6
      {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}, // 7
      {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 8
      {-0.125,    0.,-0.125,    0.,    0.,    0.,  0.25,   0.5,    0.,   0.5}  // 9
    }
  };



Real Tet10::embedding_matrix (const unsigned int i,
                              const unsigned int j,
                              const unsigned int k) const
{
  // Choose an optimal diagonal, if one has not already been selected
  this->choose_diagonal();

  // Permuted j and k indices
  unsigned int
    jp=j,
    kp=k;

  if ((i>3) && (this->_diagonal_selection!=DIAG_02_13))
    {
      // Just the enum value cast to an unsigned int...
      const unsigned ds = static_cast<unsigned int>(this->_diagonal_selection); // == 1 or 2

      // Instead of doing a lot of arithmetic, use these
      // straightforward arrays for the permutations.  Note that 3 ->
      // 3, and the first array consists of "forward" permutations of
      // the sets {0,1,2}, {4,5,6}, and {7,8,9} while the second array
      // consists of "reverse" permutations of the same sets.
      const unsigned int perms[2][10] =
        {
          {1, 2, 0, 3, 5, 6, 4, 8, 9, 7},
          {2, 0, 1, 3, 6, 4, 5, 9, 7, 8}
        };

      // Permute j
      jp = perms[ds-1][j];
      //      if (jp<3)
      //        jp = (jp+ds)%3;
      //      else if (jp>3)
      //        jp = (jp-1+ds)%3 + 1 + 3*((jp-1)/3);

      // Permute k
      kp = perms[ds-1][k];
      //      if (kp<3)
      //        kp = (kp+ds)%3;
      //      else if (kp>3)
      //        kp = (kp-1+ds)%3 + 1 + 3*((kp-1)/3);
    }

  // Debugging:
  // libMesh::err << "Selected diagonal " << _diagonal_selection << std::endl;
  // libMesh::err << "j=" << j << std::endl;
  // libMesh::err << "k=" << k << std::endl;
  // libMesh::err << "jp=" << jp << std::endl;
  // libMesh::err << "kp=" << kp << std::endl;

  // Call embedding matrix with permuted indices
  return this->_embedding_matrix[i][jp][kp];
}

#endif // #ifdef LIBMESH_ENABLE_AMR



Real Tet10::volume () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0), x1 = point(1), x2 = point(2), x3 = point(3), x4 = point(4),
    x5 = point(5), x6 = point(6), x7 = point(7), x8 = point(8), x9 = point(9);

  // The constant components of the dx/dxi vector, linear in xi, eta, zeta.
  // These were copied directly from the output of a Python script.
  Point dx_dxi[4] =
    {
      -3*x0 - x1 + 4*x4,         // constant
      4*x0 - 4*x4 - 4*x7 + 4*x8, // zeta
      4*x0 - 4*x4 + 4*x5 - 4*x6, // eta
      4*x0 + 4*x1 - 8*x4         // xi
    };

  // The constant components of the dx/deta vector, linear in xi, eta, zeta.
  // These were copied directly from the output of a Python script.
  Point dx_deta[4] =
    {
      -3*x0 - x2 + 4*x6,         // constant
      4*x0 - 4*x6 - 4*x7 + 4*x9, // zeta
      4*x0 + 4*x2 - 8*x6,        // eta
      4*x0 - 4*x4 + 4*x5 - 4*x6  // xi
    };

  // The constant components of the dx/dzeta vector, linear in xi, eta, zeta.
  // These were copied directly from the output of a Python script.
  Point dx_dzeta[4] =
    {
      -3*x0 - x3 + 4*x7,         // constant
      4*x0 + 4*x3 - 8*x7,        // zeta
      4*x0 - 4*x6 - 4*x7 + 4*x9, // eta
      4*x0 - 4*x4 - 4*x7 + 4*x8  // xi
    };

  // 2x2x2 conical quadrature rule.  Note: there is also a five point
  // rule for tets with a negative weight which would be cheaper, but
  // we'll use this one to preclude any possible issues with
  // cancellation error.
  const int N = 8;
  static const Real w[N] =
    {
      3.6979856358852914509238091810505e-02_R,
      1.6027040598476613723156741868689e-02_R,
      2.1157006454524061178256145400082e-02_R,
      9.1694299214797439226823542540576e-03_R,
      3.6979856358852914509238091810505e-02_R,
      1.6027040598476613723156741868689e-02_R,
      2.1157006454524061178256145400082e-02_R,
      9.1694299214797439226823542540576e-03_R
    };

  static const Real xi[N] =
    {
      1.2251482265544137786674043037115e-01_R,
      5.4415184401122528879992623629551e-01_R,
      1.2251482265544137786674043037115e-01_R,
      5.4415184401122528879992623629551e-01_R,
      1.2251482265544137786674043037115e-01_R,
      5.4415184401122528879992623629551e-01_R,
      1.2251482265544137786674043037115e-01_R,
      5.4415184401122528879992623629551e-01_R
    };

  static const Real eta[N] =
    {
      1.3605497680284601717109468420738e-01_R,
      7.0679724159396903069267439165167e-02_R,
      5.6593316507280088053551297149570e-01_R,
      2.9399880063162286589079157179842e-01_R,
      1.3605497680284601717109468420738e-01_R,
      7.0679724159396903069267439165167e-02_R,
      5.6593316507280088053551297149570e-01_R,
      2.9399880063162286589079157179842e-01_R
    };

  static const Real zeta[N] =
    {
      1.5668263733681830907933725249176e-01_R,
      8.1395667014670255076709592007207e-02_R,
      6.5838687060044409936029672711329e-02_R,
      3.4202793236766414300604458388142e-02_R,
      5.8474756320489429588282763292971e-01_R,
      3.0377276481470755305409673253211e-01_R,
      2.4571332521171333166171692542182e-01_R,
      1.2764656212038543100867773351792e-01_R
    };

  Real vol = 0.;
  for (int q=0; q<N; ++q)
    {
      // Compute dx_dxi, dx_deta, dx_dzeta at the current quadrature point.
      Point
        dx_dxi_q   = dx_dxi[0]   + zeta[q]*dx_dxi[1]   + eta[q]*dx_dxi[2]   + xi[q]*dx_dxi[3],
        dx_deta_q  = dx_deta[0]  + zeta[q]*dx_deta[1]  + eta[q]*dx_deta[2]  + xi[q]*dx_deta[3],
        dx_dzeta_q = dx_dzeta[0] + zeta[q]*dx_dzeta[1] + eta[q]*dx_dzeta[2] + xi[q]*dx_dzeta[3];

      // Compute scalar triple product, multiply by weight, and accumulate volume.
      vol += w[q] * triple_product(dx_dxi_q, dx_deta_q, dx_dzeta_q);
    }

  return vol;
}


void Tet10::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 12);

  const unsigned int side = perm_num % 4;
  const unsigned int rotate = perm_num / 4;

  for (unsigned int i = 0; i != rotate; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(4,5,6);
      swap3nodes(7,8,9);
      swap3neighbors(1,2,3);
    }

  switch (side) {
  case 0:
    break;
  case 1:
    swap3nodes(0,2,3);
    swap3nodes(4,5,8);
    swap3nodes(6,9,7);
    swap3neighbors(0,2,1);
    break;
  case 2:
    swap3nodes(2,0,3);
    swap3nodes(5,4,8);
    swap3nodes(6,7,9);
    swap3neighbors(0,1,2);
    break;
  case 3:
    swap3nodes(2,1,3);
    swap3nodes(5,8,9);
    swap3nodes(6,4,7);
    swap3neighbors(0,1,3);
    break;
  default:
    libmesh_error();
  }
}


void Tet10::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,2);
  swap2nodes(4,5);
  swap2nodes(7,9);
  swap2neighbors(1,2);
  swap2boundarysides(1,2,boundary_info);
  swap2boundaryedges(0,1,boundary_info);
  swap2boundaryedges(3,5,boundary_info);
}


ElemType Tet10::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, 4);
  return TRI6;
}


} // namespace libMesh
