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
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tri3 class static member initializations
const int Tri3::num_nodes;
const int Tri3::nodes_per_side;

const unsigned int Tri3::side_nodes_map[Tri3::num_sides][Tri3::nodes_per_side] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 0}  // Side 2
  };

#ifdef LIBMESH_ENABLE_AMR

const Real Tri3::_embedding_matrix[Tri3::num_children][Tri3::num_nodes][Tri3::num_nodes] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2
      {1.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0}, // 1
      {0.5, 0.0, 0.5}  // 2
    },

    // embedding matrix for child 1
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 1.0, 0.0}, // 1
      {0.0, 0.5, 0.5}  // 2
    },

    // embedding matrix for child 2
    {
      // 0    1    2
      {0.5, 0.0, 0.5}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.0, 0.0, 1.0}  // 2
    },

    // embedding matrix for child 3
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.5, 0.0, 0.5}  // 2
    }
  };

#endif



// ------------------------------------------------------------
// Tri3 class member functions

bool Tri3::is_vertex(const unsigned int libmesh_dbg_var(n)) const
{
  libmesh_assert_not_equal_to (n, invalid_uint);
  return true;
}

bool Tri3::is_edge(const unsigned int) const
{
  return false;
}

bool Tri3::is_face(const unsigned int) const
{
  return false;
}

bool Tri3::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Tri3::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Tri3::nodes_on_edge(const unsigned int e) const
{
  return nodes_on_side(e);
}

Order Tri3::default_order() const
{
  return FIRST;
}

bool Tri3::has_invertible_map(Real tol) const
{
  return this->volume() > tol;
}

std::unique_ptr<Elem> Tri3::build_side_ptr (const unsigned int i)
{
  return this->simple_build_side_ptr<Edge2, Tri3>(i);
}



void Tri3::build_side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  this->simple_build_side_ptr<Tri3>(side, i, EDGE2);
}



void Tri3::connectivity(const unsigned int libmesh_dbg_var(sf),
                        const IOPackage iop,
                        std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(4);
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(2)+1;
        return;
      }

    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}






Point Tri3::true_centroid () const
{
  return Elem::vertex_average();
}

Real Tri3::volume () const
{
  // 3-node triangles have the following formula for computing the area
  return 0.5 * cross_norm(point(1) - point(0),
                          point(2) - point(0));
}



std::pair<Real, Real> Tri3::min_and_max_angle() const
{
  Point v10 ( this->point(1) - this->point(0) );
  Point v20 ( this->point(2) - this->point(0) );
  Point v21 ( this->point(2) - this->point(1) );

  const Real
    len_10=v10.norm(),
    len_20=v20.norm(),
    len_21=v21.norm();

  const Real
    theta0=std::acos(( v10*v20)/len_10/len_20),
    theta1=std::acos((-v10*v21)/len_10/len_21),
    theta2=libMesh::pi - theta0 - theta1
    ;

  libmesh_assert_greater (theta0, 0.);
  libmesh_assert_greater (theta1, 0.);
  libmesh_assert_greater (theta2, 0.);

  return std::make_pair(std::min(theta0, std::min(theta1,theta2)),
                        std::max(theta0, std::max(theta1,theta2)));
}

bool Tri3::contains_point (const Point & p, Real tol) const
{
  // See "Barycentric Technique" section at
  // http://www.blackpawn.com/texts/pointinpoly for details.

  // Compute vectors
  Point v0 = this->point(1) - this->point(0);
  Point v1 = this->point(2) - this->point(0);
  Point v2 = p - this->point(0);

  // Compute dot products
  Real dot00 = v0 * v0;
  Real dot01 = v0 * v1;
  Real dot02 = v0 * v2;
  Real dot11 = v1 * v1;
  Real dot12 = v1 * v2;

  // Out of plane check
  if (std::abs(triple_product(v2, v0, v1)) / std::max(dot00, dot11) > tol)
    return false;

  // Compute barycentric coordinates
  Real invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  Real u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  Real v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // Check if point is in triangle
  return (u > -tol) && (v > -tol) && (u + v < 1 + tol);
}

BoundingBox
Tri3::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}


void Tri3::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 3);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap3nodes(0,1,2);
      swap3neighbors(0,1,2);
    }
}


void Tri3::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2neighbors(1,2);
  swap2boundarysides(1,2,boundary_info);
  swap2boundaryedges(1,2,boundary_info);
}


ElemType
Tri3::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, 3);
  return EDGE2;
}

} // namespace libMesh
