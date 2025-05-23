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
#include "libmesh/cell_tet.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_elem_quality.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tet class static member initializations
const int Tet::num_sides;
const int Tet::num_edges;
const int Tet::num_children;

const Real Tet::_master_points[14][3] =
  {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {0.5, 0, 0},
    {0.5, 0.5, 0},
    {0, 0.5, 0},
    {0, 0, 0.5},
    {0.5, 0, 0.5},
    {0, 0.5, 0.5},
    {1/Real(3), 1/Real(3), 0},
    {1/Real(3), 0, 1/Real(3)},
    {1/Real(3), 1/Real(3), 1/Real(3)},
    {0, 1/Real(3), 1/Real(3)}
  };

const unsigned int Tet::edge_sides_map[6][2] =
  {
    {0, 1}, // Edge 0
    {0, 2}, // Edge 1
    {0, 3}, // Edge 2
    {1, 3}, // Edge 3
    {1, 2}, // Edge 4
    {2, 3}  // Edge 5
  };

const unsigned int Tet::adjacent_edges_map[/*num_vertices*/4][/*n_adjacent_edges*/3] =
  {
    {0, 2, 3},  // Edges adjacent to node 0
    {0, 1, 4},  // Edges adjacent to node 1
    {1, 2, 5},  // Edges adjacent to node 2
    {3, 4, 5},  // Edges adjacent to node 3
  };

// ------------------------------------------------------------
// Tet class member functions
dof_id_type Tet::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Tet4::side_nodes_map[s][0]),
                           this->node_id(Tet4::side_nodes_map[s][1]),
                           this->node_id(Tet4::side_nodes_map[s][2]));
}



dof_id_type Tet::low_order_key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Tet4::side_nodes_map[s][0]),
                           this->node_id(Tet4::side_nodes_map[s][1]),
                           this->node_id(Tet4::side_nodes_map[s][2]));
}



unsigned int Tet::local_side_node(unsigned int side,
                                  unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Tet4::nodes_per_side);

  return Tet4::side_nodes_map[side][side_node];
}



unsigned int Tet::local_edge_node(unsigned int edge,
                                  unsigned int edge_node) const
{
  libmesh_assert_less (edge, this->n_edges());
  libmesh_assert_less (edge_node, Tet4::nodes_per_edge);

  return Tet4::edge_nodes_map[edge][edge_node];
}


std::unique_ptr<Elem> Tet::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face = std::make_unique<Tri3>();

  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(Tet4::side_nodes_map[i][n]));

  return face;
}



void Tet::side_ptr (std::unique_ptr<Elem> & side,
                    const unsigned int i)
{
  this->simple_side_ptr<Tet,Tet4>(side, i, TRI3);
}



void Tet::select_diagonal (const Diagonal diag) const
{
  libmesh_assert_equal_to (_diagonal_selection, INVALID_DIAG);
  _diagonal_selection = diag;
}





#ifdef LIBMESH_ENABLE_AMR


bool Tet::is_child_on_side_helper(const unsigned int c,
                                  const unsigned int s,
                                  const unsigned int checked_nodes[][3]) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  // For the 4 vertices, child c touches vertex c, so we can return
  // true if that vertex is on side s
  for (unsigned int i = 0; i != 3; ++i)
    if (Tet4::side_nodes_map[s][i] == c)
      return true;

  // If we are a "vertex child" and we didn't already return true,
  // we must not be on the side in question
  if (c < 4)
    return false;

  // For the 4 non-vertex children, the child ordering depends on the
  // diagonal selection.  We'll let the embedding matrix figure that
  // out: if this child has three nodes that don't depend on the
  // position of the node_facing_side[s], then we're on side s.  Which
  // three nodes those are depends on the subclass, so their responsibility
  // is to call this function with the proper check_nodes array
  const unsigned int node_facing_side[4] = {3, 2, 0, 1};
  const unsigned int n = node_facing_side[s];

  // Add up the absolute values of the entries of the embedding matrix for the
  // nodes opposite node n.  If it is equal to zero, then the child in question is
  // on side s, so return true.
  Real embedding_sum = 0.;
  for (unsigned i=0; i<3; ++i)
    embedding_sum += std::abs(this->embedding_matrix(c, checked_nodes[n][i], n));

  return ( std::abs(embedding_sum) < 1.e-3 );
}

#else

bool Tet::is_child_on_side_helper(const unsigned int /*c*/,
                                  const unsigned int /*s*/,
                                  const unsigned int /*checked_nodes*/[][3]) const
{
  libmesh_not_implemented();
  return false;
}

#endif //LIBMESH_ENABLE_AMR




void Tet::choose_diagonal() const
{
  // If uninitialized diagonal selection, select the shortest octahedron diagonal
  if (this->_diagonal_selection==INVALID_DIAG)
    {
      Real diag_01_23 = (this->point(0) + this->point(1) - this->point(2) - this->point(3)).norm_sq();
      Real diag_02_13 = (this->point(0) - this->point(1) + this->point(2) - this->point(3)).norm_sq();
      Real diag_03_12 = (this->point(0) - this->point(1) - this->point(2) + this->point(3)).norm_sq();

      std::array<Real, 3> D = {diag_02_13, diag_03_12, diag_01_23};

      this->_diagonal_selection =
        Diagonal(std::distance(D.begin(), std::min_element(D.begin(), D.end())));
    }
}



bool Tet::is_edge_on_side(const unsigned int e,
                          const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (edge_sides_map[e][0] == s || edge_sides_map[e][1] == s);
}



std::vector<unsigned int> Tet::sides_on_edge(const unsigned int e) const
{
  libmesh_assert_less (e, this->n_edges());
  return {edge_sides_map[e][0], edge_sides_map[e][1]};
}



bool
Tet::is_flipped() const
{
  return (triple_product(this->point(1)-this->point(0),
                        this->point(2)-this->point(0),
                        this->point(3)-this->point(0)) < 0);
}


std::vector<unsigned int>
Tet::edges_adjacent_to_node(const unsigned int n) const
{
  libmesh_assert_less(n, this->n_nodes());

  // For vertices, we use the Tet::adjacent_sides_map, otherwise each
  // of the mid-edge nodes is adjacent only to the edge it is on, and the
  // mid-face nodes are not adjacent to any edges.
  if (this->is_vertex(n))
    return {std::begin(adjacent_edges_map[n]), std::end(adjacent_edges_map[n])};
  else if (this->is_edge(n))
    return {n - this->n_vertices()};

  // Current Tets have only vertex, edge, and face nodes.
  libmesh_assert(this->is_face(n));
  return {};
}


Real Tet::quality(const ElemQuality q) const
{
  return Elem::quality(q); // Not implemented
}




std::pair<Real, Real> Tet::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {

    case ASPECT_RATIO_BETA:
    case ASPECT_RATIO_GAMMA:
      bounds.first  = 1.;
      bounds.second = 3.;
      break;

    case SIZE:
    case SHAPE:
      bounds.first  = 0.2;
      bounds.second = 1.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 3.;
      break;

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;

    case JACOBIAN:
    case SCALED_JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.414;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}



bool Tet::on_reference_element(const Point & p,
                               const Real eps) const
{
  const Real & xi = p(0);
  const Real & eta = p(1);
  const Real & zeta = p(2);
        // The reference tetrahedral is isosceles
        // and is bound by xi=0, eta=0, zeta=0,
        // and xi+eta+zeta=1.
  return ((xi   >= 0.-eps) &&
          (eta  >= 0.-eps) &&
          (zeta >= 0.-eps) &&
          ((xi + eta + zeta) <= 1.+eps));
}


} // namespace libMesh
