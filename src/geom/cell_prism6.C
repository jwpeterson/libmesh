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
#include "libmesh/cell_prism6.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/fe_lagrange_shape_1D.h"

// C++ includes
#include <array>

// Helper functions for computing volume() and true_centroid()
namespace {

using namespace libMesh;

// Templated numerical quadrature implementation function used by both
// Prism6::volume() and Prism6::true_centroid() routines.
template <typename T>
void cell_integral(const Elem * elem, T & t);

/**
 * Object passed to cell_integral() routine for computing the cell volume.
 */
struct VolumeIntegral
{
  // API
  void accumulate_qp(Real /*xi*/, Real /*eta*/, Real /*zeta*/, Real JxW)
  {
    vol += JxW;
  };

  Real vol = 0.;
};

/**
 * Object passed to cell_integral() routine for computing the nodal
 * volumes (one value per node)
 */
struct NodalVolumeIntegral
{
  // API
  void accumulate_qp(Real xi, Real eta, Real zeta, Real JxW)
  {
    // Copied from fe_lagrange_shape_3D.C
    static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
    static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};

    // Copied from fe_lagrange_shape_2D.C
    auto tri3_phi = [](int i, Real x, Real y)
    {
      switch(i)
      {
        case 0:
          return 1 - x - y;

        case 1:
          return x;

        case 2:
          return y;

        default:
          libmesh_error_msg("Invalid shape function index i = " << i);
      }
    };

    for (int i=0; i<6; ++i)
      {
        Real phi =
          tri3_phi(i1[i], xi, eta) * fe_lagrange_1D_linear_shape(i0[i], zeta);

        V[i] += phi * JxW;
      }
  };

  // nodal volumes
  std::array<Real, 6> V {};
};

}

namespace libMesh
{



// ------------------------------------------------------------
// Prism6 class static member initializations
const int Prism6::num_nodes;
const int Prism6::nodes_per_side;
const int Prism6::nodes_per_edge;

const unsigned int Prism6::side_nodes_map[Prism6::num_sides][Prism6::nodes_per_side] =
  {
    {0, 2, 1, 99}, // Side 0
    {0, 1, 4,  3}, // Side 1
    {1, 2, 5,  4}, // Side 2
    {2, 0, 3,  5}, // Side 3
    {3, 4, 5, 99}  // Side 4
  };

const unsigned int Prism6::side_elems_map[Prism6::num_sides][Prism6::nodes_per_side] =
  {
    {0, 1, 2, 3}, // Side 0
    {0, 1, 4, 5}, // Side 1
    {1, 2, 5, 6}, // Side 2
    {0, 2, 4, 6}, // Side 3
    {4, 5, 6, 7}  // Side 4
  };

const unsigned int Prism6::edge_nodes_map[Prism6::num_edges][Prism6::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {0, 2}, // Edge 2
    {0, 3}, // Edge 3
    {1, 4}, // Edge 4
    {2, 5}, // Edge 5
    {3, 4}, // Edge 6
    {4, 5}, // Edge 7
    {3, 5}  // Edge 8
  };

// ------------------------------------------------------------
// Prism6 class member functions

bool Prism6::is_vertex(const unsigned int libmesh_dbg_var(n)) const
{
  libmesh_assert_not_equal_to (n, invalid_uint);
  return true;
}

bool Prism6::is_edge(const unsigned int) const
{
  return false;
}


bool Prism6::is_face(const unsigned int) const
{
  return false;
}

bool Prism6::is_node_on_side(const unsigned int n,
                             const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Prism6::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s > 0 && s < 4) ? 0 : 1;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

std::vector<unsigned>
Prism6::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool Prism6::is_node_on_edge(const unsigned int n,
                             const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



bool Prism6::has_affine_map() const
{
  // Make sure z edges are affine
  Point v = this->point(3) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(1), affine_tol) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(2), affine_tol))
    return false;
  return true;
}



Order Prism6::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> Prism6::build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  switch (i)
    {
    case 0: // the triangular face at z=-1
    case 4: // the triangular face at z=1
      {
        face = std::make_unique<Tri3>();
        break;
      }
    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        face = std::make_unique<Quad4>();
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(Prism6::side_nodes_map[i][n]));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
}



void Prism6::build_side_ptr (std::unique_ptr<Elem> & side,
                             const unsigned int i)
{
  this->side_ptr(side, i);
  side->set_interior_parent(this);
  side->inherit_data_from(*this);
}



std::unique_ptr<Elem> Prism6::build_edge_ptr (const unsigned int i)
{
  return this->simple_build_edge_ptr<Edge2,Prism6>(i);
}



void Prism6::build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i)
{
  this->simple_build_edge_ptr<Prism6>(edge, i, EDGE2);
}



void Prism6::connectivity(const unsigned int libmesh_dbg_var(sc),
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
        conn[3] = this->node_id(2)+1;
        conn[4] = this->node_id(3)+1;
        conn[5] = this->node_id(4)+1;
        conn[6] = this->node_id(5)+1;
        conn[7] = this->node_id(5)+1;
        return;
      }

    case VTK:
      {
        conn.resize(6);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(2);
        conn[2] = this->node_id(1);
        conn[3] = this->node_id(3);
        conn[4] = this->node_id(5);
        conn[5] = this->node_id(4);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



#ifdef LIBMESH_ENABLE_AMR

const Real Prism6::_embedding_matrix[Prism6::num_children][Prism6::num_nodes][Prism6::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //  0     1     2     3     4     5
      { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 3
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 4
      { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
    },

    // embedding matrix for child 1
    {
      //  0     1     2     3     4     5
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
      { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 4
      { 0.0,  .25,  .25,  0.0,  .25,  .25}  // 5
    },

    // embedding matrix for child 2
    {
      //  0     1     2     3     4     5
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0}, // 2
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 3
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}  // 5
    },

    // embedding matrix for child 3
    {
      //  0     1     2     3     4     5
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
      { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
    },

    // embedding matrix for child 4
    {
      //  0     1     2     3     4     5
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 0
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 1
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
      { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
    },

    // embedding matrix for child 5
    {
      //  0     1     2     3     4     5
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
      { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 1
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}  // 5
    },

    // embedding matrix for child 6
    {
      //  0     1     2     3     4     5
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 0
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // 5
    },

    // embedding matrix for child 7
    {
      //  0     1     2     3     4     5
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
    }
  };

#endif



Point Prism6::true_centroid () const
{
  NodalVolumeIntegral nvi;
  cell_integral(this, nvi);

  // Compute centroid
  Point cp;
  Real vol = 0.;

  for (int i=0; i<6; ++i)
    {
      cp += this->point(i) * nvi.V[i];
      vol += nvi.V[i];
    }

  return cp / vol;
}

Real Prism6::volume () const
{
  VolumeIntegral vi;
  cell_integral(this, vi);
  return vi.vol;
}

BoundingBox
Prism6::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}


void
Prism6::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 6);
  const unsigned int side = perm_num % 2;
  const unsigned int rotate = perm_num / 2;

  for (unsigned int i = 0; i != rotate; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(3,4,5);
      swap3neighbors(1,2,3);
    }

  switch (side) {
  case 0:
    break;
  case 1:
    swap2nodes(1,3);
    swap2nodes(0,4);
    swap2nodes(2,5);
    swap2neighbors(0,4);
    swap2neighbors(2,3);
    break;
  default:
    libmesh_error();
  }

}


void
Prism6::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(3,4);
  swap2neighbors(2,3);
  swap2boundarysides(2,3,boundary_info);
  swap2boundaryedges(0,1,boundary_info);
  swap2boundaryedges(3,4,boundary_info);
  swap2boundaryedges(7,8,boundary_info);
}


ElemType
Prism6::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 5);
  if (s == 0 || s == 4)
    return TRI3;
  return QUAD4;
}

} // namespace libMesh

namespace // anonymous
{

template <typename T>
void cell_integral(const Elem * elem, T & t)
{
  // Conveient references to our points.
  const Point
    &x0 = elem->point(0), &x1 = elem->point(1), &x2 = elem->point(2),
    &x3 = elem->point(3), &x4 = elem->point(4), &x5 = elem->point(5);

  // constant and zeta terms only.  These are copied directly from a
  // Python script.
  Point dx_dxi[2] =
    {
      -x0/2 + x1/2 - x3/2 + x4/2, // constant
      x0/2 - x1/2 - x3/2 + x4/2,  // zeta
    };

  // constant and zeta terms only.  These are copied directly from a
  // Python script.
  Point dx_deta[2] =
    {
      -x0/2 + x2/2 - x3/2 + x5/2, // constant
      x0/2 - x2/2 - x3/2 + x5/2,  // zeta
    };

  // Constant, xi, and eta terms
  Point dx_dzeta[3] =
    {
      -x0/2 + x3/2,              // constant
      x0/2 - x2/2 - x3/2 + x5/2, // eta
      x0/2 - x1/2 - x3/2 + x4/2  // xi
    };

  // The quadrature rule for the Prism6 is a tensor product between a
  // four-point TRI3 rule (in xi, eta) and a two-point EDGE2 rule (in
  // zeta) which is capable of integrating cubics exactly.

  // Number of points in the 2D quadrature rule.
  const int N2D = 4;

  static const Real w2D[N2D] =
    {
      1.5902069087198858469718450103758e-01_R,
      9.0979309128011415302815498962418e-02_R,
      1.5902069087198858469718450103758e-01_R,
      9.0979309128011415302815498962418e-02_R
    };

  static const Real xi[N2D] =
    {
      1.5505102572168219018027159252941e-01_R,
      6.4494897427831780981972840747059e-01_R,
      1.5505102572168219018027159252941e-01_R,
      6.4494897427831780981972840747059e-01_R
    };

  static const Real eta[N2D] =
    {
      1.7855872826361642311703513337422e-01_R,
      7.5031110222608118177475598324603e-02_R,
      6.6639024601470138670269327409637e-01_R,
      2.8001991549907407200279599420481e-01_R
    };

  // Number of points in the 1D quadrature rule.  The weights of the
  // 1D quadrature rule are equal to 1.
  const int N1D = 2;

  // Points of the 1D quadrature rule
  static const Real zeta[N1D] =
    {
      -std::sqrt(Real(3))/3,
      std::sqrt(Real(3))/3.
    };

  for (int i=0; i<N2D; ++i)
    {
      // dx_dzeta depends only on the 2D quadrature rule points.
      Point dx_dzeta_q = dx_dzeta[0] + eta[i]*dx_dzeta[1] + xi[i]*dx_dzeta[2];

      for (int j=0; j<N1D; ++j)
        {
          // dx_dxi and dx_deta only depend on the 1D quadrature rule points.
          Point
            dx_dxi_q  = dx_dxi[0]  + zeta[j]*dx_dxi[1],
            dx_deta_q = dx_deta[0] + zeta[j]*dx_deta[1];

          // Compute t-specific contributions at current qp
          t.accumulate_qp(xi[i], eta[i], zeta[j], w2D[i] * triple_product(dx_dxi_q, dx_deta_q, dx_dzeta_q));
        }
    }
}

}
