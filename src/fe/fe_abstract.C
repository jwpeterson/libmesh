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

// libmesh includes
#include "libmesh/fe.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/remote_elem.h"
#include "libmesh/tensor_value.h"
#include "libmesh/threads.h"
#include "libmesh/enum_to_string.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "libmesh/inf_fe.h"
#include "libmesh/fe_interface_macros.h"
#endif

namespace libMesh
{

FEAbstract::FEAbstract(const unsigned int d,
                       const FEType & fet) :
  _fe_map( FEMap::build(fet) ),
  dim(d),
  calculations_started(false),
  calculate_dual(false),
  calculate_default_dual_coeff(true),
  calculate_nothing(false),
  calculate_map(false),
  calculate_phi(false),
  calculate_dphi(false),
  calculate_d2phi(false),
  calculate_curl_phi(false),
  calculate_div_phi(false),
  calculate_dphiref(false),
  fe_type(fet),
  _elem_type(INVALID_ELEM),
  _elem(nullptr),
  _elem_p_level(0),
  _p_level(0),
  qrule(nullptr),
  shapes_on_quadrature(false),
  _n_total_qp(0),
  _add_p_level_in_reinit(true)
{
}


FEAbstract::~FEAbstract() = default;


std::unique_ptr<FEAbstract> FEAbstract::build(const unsigned int dim,
                                              const FEType & fet)
{
  switch (dim)
    {
      // 0D
    case 0:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return std::make_unique<FE<0,CLOUGH>>(fet);

          case HERMITE:
            return std::make_unique<FE<0,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<0,LAGRANGE>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FE<0,LAGRANGE_VEC>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<0,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FE<0,L2_LAGRANGE_VEC>>(fet);

          case HIERARCHIC_VEC:
            return std::make_unique<FE<0,HIERARCHIC_VEC>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<0,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<0,L2_HIERARCHIC>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FE<0,L2_HIERARCHIC_VEC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<0,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<0,MONOMIAL>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FE<0,MONOMIAL_VEC>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<0,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<0,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<0,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<0>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<0>>(fet);

          case NEDELEC_ONE:
            return std::make_unique<FENedelecOne<0>>(fet);

          case RAVIART_THOMAS:
            return std::make_unique<FERaviartThomas<0>>(fet);

          case L2_RAVIART_THOMAS:
            return std::make_unique<FEL2RaviartThomas<0>>(fet);

          case SUBDIVISION:
            return std::make_unique<FESubdivision>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << Utility::enum_to_string(fet.family));
          }
      }
      // 1D
    case 1:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return std::make_unique<FE<1,CLOUGH>>(fet);

          case HERMITE:
            return std::make_unique<FE<1,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<1,LAGRANGE>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FE<1,LAGRANGE_VEC>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<1,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FE<1,L2_LAGRANGE_VEC>>(fet);

          case HIERARCHIC_VEC:
            return std::make_unique<FE<1,HIERARCHIC_VEC>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<1,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<1,L2_HIERARCHIC>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FE<1,L2_HIERARCHIC_VEC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<1,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<1,MONOMIAL>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FE<1,MONOMIAL_VEC>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<1,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<1,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<1,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<1>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<1>>(fet);

          case NEDELEC_ONE:
            return std::make_unique<FENedelecOne<1>>(fet);

          case RAVIART_THOMAS:
            return std::make_unique<FERaviartThomas<1>>(fet);

          case L2_RAVIART_THOMAS:
            return std::make_unique<FEL2RaviartThomas<1>>(fet);

          case SUBDIVISION:
            return std::make_unique<FESubdivision>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << Utility::enum_to_string(fet.family));
          }
      }


      // 2D
    case 2:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return std::make_unique<FE<2,CLOUGH>>(fet);

          case HERMITE:
            return std::make_unique<FE<2,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<2,LAGRANGE>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FE<2,LAGRANGE_VEC>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<2,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FE<2,L2_LAGRANGE_VEC>>(fet);

          case HIERARCHIC_VEC:
            return std::make_unique<FE<2,HIERARCHIC_VEC>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<2,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<2,L2_HIERARCHIC>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FE<2,L2_HIERARCHIC_VEC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<2,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<2,MONOMIAL>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FE<2,MONOMIAL_VEC>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<2,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<2,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<2,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<2>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<2>>(fet);

          case NEDELEC_ONE:
            return std::make_unique<FENedelecOne<2>>(fet);

          case RAVIART_THOMAS:
            return std::make_unique<FERaviartThomas<2>>(fet);

          case L2_RAVIART_THOMAS:
            return std::make_unique<FEL2RaviartThomas<2>>(fet);

          case SUBDIVISION:
            return std::make_unique<FESubdivision>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << Utility::enum_to_string(fet.family));
          }
      }


      // 3D
    case 3:
      {
        switch (fet.family)
          {
          case CLOUGH:
            libmesh_error_msg("ERROR: Clough-Tocher elements currently only support 1D and 2D");

          case HERMITE:
            return std::make_unique<FE<3,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<3,LAGRANGE>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FE<3,LAGRANGE_VEC>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<3,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FE<3,L2_LAGRANGE_VEC>>(fet);

          case HIERARCHIC_VEC:
            return std::make_unique<FE<3,HIERARCHIC_VEC>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<3,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<3,L2_HIERARCHIC>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FE<3,L2_HIERARCHIC_VEC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<3,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<3,MONOMIAL>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FE<3,MONOMIAL_VEC>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<3,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<3,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<3,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<3>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<3>>(fet);

          case NEDELEC_ONE:
            return std::make_unique<FENedelecOne<3>>(fet);

          case RAVIART_THOMAS:
            return std::make_unique<FERaviartThomas<3>>(fet);

          case L2_RAVIART_THOMAS:
            return std::make_unique<FEL2RaviartThomas<3>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << Utility::enum_to_string(fet.family));
          }
      }

    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    }
}



void FEAbstract::get_refspace_nodes(const ElemType itemType, std::vector<Point> & nodes)
{
  const unsigned int n_nodes = Elem::type_to_n_nodes_map[itemType];
  if (n_nodes == invalid_uint)
    libmesh_error_msg("Number of nodes is not well-defined for " <<
                      Utility::enum_to_string(itemType));

  nodes.resize(n_nodes);
  switch(itemType)
    {
    case NODEELEM:
      {
        nodes[0] = Point (0.,0.,0.);
        return;
      }
    case EDGE3:
      {
        nodes[2] = Point (0.,0.,0.);
        libmesh_fallthrough();
      }
    case EDGE2:
      {
        nodes[0] = Point (-1.,0.,0.);
        nodes[1] = Point (1.,0.,0.);
        return;
      }
    case EDGE4: // not nested with EDGE3
      {
        nodes[0] = Point (-1.,0.,0.);
        nodes[1] = Point (1.,0.,0.);
        nodes[2] = Point (-1./3.,0.,0.);
        nodes[3] - Point (1./3.,0.,0.);
        return;
      }
    case TRI7:
      {
        nodes[6] = Point (1./3.,1./3.,0.);
        libmesh_fallthrough();
      }
    case TRI6:
      {
        nodes[3] = Point (.5,0.,0.);
        nodes[4] = Point (.5,.5,0.);
        nodes[5] = Point (0.,.5,0.);
        libmesh_fallthrough();
      }
    case TRI3:
    case TRISHELL3:
      {
        nodes[0] = Point (0.,0.,0.);
        nodes[1] = Point (1.,0.,0.);
        nodes[2] = Point (0.,1.,0.);
        return;
      }
    case QUAD9:
    case QUADSHELL9:
      {
        nodes[8] = Point (0.,0.,0.);
        libmesh_fallthrough();
      }
    case QUAD8:
    case QUADSHELL8:
      {
        nodes[4] = Point (0.,-1.,0.);
        nodes[5] = Point (1.,0.,0.);
        nodes[6] = Point (0.,1.,0.);
        nodes[7] = Point (-1.,0.,0.);
        libmesh_fallthrough();
      }
    case QUAD4:
    case QUADSHELL4:
      {
        nodes[0] = Point (-1.,-1.,0.);
        nodes[1] = Point (1.,-1.,0.);
        nodes[2] = Point (1.,1.,0.);
        nodes[3] = Point (-1.,1.,0.);
        return;
      }
    case TET14:
      {
        nodes[10] = Point (1/Real(3),1/Real(3),0.);
        nodes[11] = Point (1/Real(3),0.,1/Real(3));
        nodes[12] = Point (1/Real(3),1/Real(3),1/Real(3));
        nodes[13] = Point (0.,1/Real(3),1/Real(3));
        libmesh_fallthrough();
      }
    case TET10:
      {
        nodes[4] = Point (.5,0.,0.);
        nodes[5] = Point (.5,.5,0.);
        nodes[6] = Point (0.,.5,0.);
        nodes[7] = Point (0.,0.,.5);
        nodes[8] = Point (.5,0.,.5);
        nodes[9] = Point (0.,.5,.5);
        libmesh_fallthrough();
      }
    case TET4:
      {
        nodes[0] = Point (0.,0.,0.);
        nodes[1] = Point (1.,0.,0.);
        nodes[2] = Point (0.,1.,0.);
        nodes[3] = Point (0.,0.,1.);
        return;
      }
    case HEX27:
      {
        nodes[20] = Point (0.,0.,-1.);
        nodes[21] = Point (0.,-1.,0.);
        nodes[22] = Point (1.,0.,0.);
        nodes[23] = Point (0.,1.,0.);
        nodes[24] = Point (-1.,0.,0.);
        nodes[25] = Point (0.,0.,1.);
        nodes[26] = Point (0.,0.,0.);
        libmesh_fallthrough();
      }
    case HEX20:
      {
        nodes[8] = Point (0.,-1.,-1.);
        nodes[9] = Point (1.,0.,-1.);
        nodes[10] = Point (0.,1.,-1.);
        nodes[11] = Point (-1.,0.,-1.);
        nodes[12] = Point (-1.,-1.,0.);
        nodes[13] = Point (1.,-1.,0.);
        nodes[14] = Point (1.,1.,0.);
        nodes[15] = Point (-1.,1.,0.);
        nodes[16] = Point (0.,-1.,1.);
        nodes[17] = Point (1.,0.,1.);
        nodes[18] = Point (0.,1.,1.);
        nodes[19] = Point (-1.,0.,1.);
        libmesh_fallthrough();
      }
    case HEX8:
      {
        nodes[0] = Point (-1.,-1.,-1.);
        nodes[1] = Point (1.,-1.,-1.);
        nodes[2] = Point (1.,1.,-1.);
        nodes[3] = Point (-1.,1.,-1.);
        nodes[4] = Point (-1.,-1.,1.);
        nodes[5] = Point (1.,-1.,1.);
        nodes[6] = Point (1.,1.,1.);
        nodes[7] = Point (-1.,1.,1.);
        return;
      }
    case PRISM21:
      {
        nodes[20] = Point (1/Real(3),1/Real(3),0);
        libmesh_fallthrough();
      }
    case PRISM20:
      {
        nodes[18] = Point (1/Real(3),1/Real(3),-1);
        nodes[19] = Point (1/Real(3),1/Real(3),1);
        libmesh_fallthrough();
      }
    case PRISM18:
      {
        nodes[15] = Point (.5,0.,0.);
        nodes[16] = Point (.5,.5,0.);
        nodes[17] = Point (0.,.5,0.);
        libmesh_fallthrough();
      }
    case PRISM15:
      {
        nodes[6] = Point (.5,0.,-1.);
        nodes[7] = Point (.5,.5,-1.);
        nodes[8] = Point (0.,.5,-1.);
        nodes[9] = Point (0.,0.,0.);
        nodes[10] = Point (1.,0.,0.);
        nodes[11] = Point (0.,1.,0.);
        nodes[12] = Point (.5,0.,1.);
        nodes[13] = Point (.5,.5,1.);
        nodes[14] = Point (0.,.5,1.);
        libmesh_fallthrough();
      }
    case PRISM6:
      {
        nodes[0] = Point (0.,0.,-1.);
        nodes[1] = Point (1.,0.,-1.);
        nodes[2] = Point (0.,1.,-1.);
        nodes[3] = Point (0.,0.,1.);
        nodes[4] = Point (1.,0.,1.);
        nodes[5] = Point (0.,1.,1.);
        return;
      }
    case PYRAMID18:
      {
        // triangle centers
        nodes[14] = Point (-2/Real(3),0.,1/Real(3));
        nodes[15] = Point (0.,2/Real(3),1/Real(3));
        nodes[16] = Point (2/Real(3),0.,1/Real(3));
        nodes[17] = Point (0.,-2/Real(3),1/Real(3));

        libmesh_fallthrough();
      }
    case PYRAMID14:
      {
        // base center
        nodes[13] = Point (0.,0.,0.);

        libmesh_fallthrough();
      }
    case PYRAMID13:
      {
        // base midedge
        nodes[5] = Point (0.,-1.,0.);
        nodes[6] = Point (1.,0.,0.);
        nodes[7] = Point (0.,1.,0.);
        nodes[8] = Point (-1,0.,0.);

        // lateral midedge
        nodes[9] = Point (-.5,-.5,.5);
        nodes[10] = Point (.5,-.5,.5);
        nodes[11] = Point (.5,.5,.5);
        nodes[12] = Point (-.5,.5,.5);

        libmesh_fallthrough();
      }
    case PYRAMID5:
      {
        // base corners
        nodes[0] = Point (-1.,-1.,0.);
        nodes[1] = Point (1.,-1.,0.);
        nodes[2] = Point (1.,1.,0.);
        nodes[3] = Point (-1.,1.,0.);
        // apex
        nodes[4] = Point (0.,0.,1.);
        return;
      }

    default:
      libmesh_error_msg("ERROR: Unknown element type " << Utility::enum_to_string(itemType));
    }
}



#ifdef LIBMESH_ENABLE_DEPRECATED
bool FEAbstract::on_reference_element(const Point & p, const ElemType t, const Real eps)
{
  // Use Elem::on_reference_element() instead
  libmesh_deprecated();

  libmesh_assert_greater_equal (eps, 0.);

  const Real xi   = p(0);
#if LIBMESH_DIM > 1
  const Real eta  = p(1);
#else
  const Real eta  = 0.;
#endif
#if LIBMESH_DIM > 2
  const Real zeta = p(2);
#else
  const Real zeta  = 0.;
#endif

  switch (t)
    {
    case NODEELEM:
      {
        return (!xi && !eta && !zeta);
      }
    case EDGE2:
    case EDGE3:
    case EDGE4:
      {
        // The reference 1D element is [-1,1].
        if ((xi >= -1.-eps) &&
            (xi <=  1.+eps))
          return true;

        return false;
      }


    case TRI3:
    case TRISHELL3:
    case TRI6:
    case TRI7:
      {
        // The reference triangle is isosceles
        // and is bound by xi=0, eta=0, and xi+eta=1.
        if ((xi  >= 0.-eps) &&
            (eta >= 0.-eps) &&
            ((xi + eta) <= 1.+eps))
          return true;

        return false;
      }


    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      {
        // The reference quadrilateral element is [-1,1]^2.
        if ((xi  >= -1.-eps) &&
            (xi  <=  1.+eps) &&
            (eta >= -1.-eps) &&
            (eta <=  1.+eps))
          return true;

        return false;
      }


    case TET4:
    case TET10:
    case TET14:
      {
        // The reference tetrahedral is isosceles
        // and is bound by xi=0, eta=0, zeta=0,
        // and xi+eta+zeta=1.
        if ((xi   >= 0.-eps) &&
            (eta  >= 0.-eps) &&
            (zeta >= 0.-eps) &&
            ((xi + eta + zeta) <= 1.+eps))
          return true;

        return false;
      }


    case HEX8:
    case HEX20:
    case HEX27:
      {
        /*
          if ((xi   >= -1.) &&
          (xi   <=  1.) &&
          (eta  >= -1.) &&
          (eta  <=  1.) &&
          (zeta >= -1.) &&
          (zeta <=  1.))
          return true;
        */

        // The reference hexahedral element is [-1,1]^3.
        if ((xi   >= -1.-eps) &&
            (xi   <=  1.+eps) &&
            (eta  >= -1.-eps) &&
            (eta  <=  1.+eps) &&
            (zeta >= -1.-eps) &&
            (zeta <=  1.+eps))
          {
            //    libMesh::out << "Strange Point:\n";
            //    p.print();
            return true;
          }

        return false;
      }

    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      {
        // Figure this one out...
        // inside the reference triangle with zeta in [-1,1]
        if ((xi   >=  0.-eps) &&
            (eta  >=  0.-eps) &&
            (zeta >= -1.-eps) &&
            (zeta <=  1.+eps) &&
            ((xi + eta) <= 1.+eps))
          return true;

        return false;
      }


    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      {
        // Check that the point is on the same side of all the faces
        // by testing whether:
        //
        // n_i.(x - x_i) <= 0
        //
        // for each i, where:
        //   n_i is the outward normal of face i,
        //   x_i is a point on face i.
        if ((-eta - 1. + zeta <= 0.+eps) &&
            (  xi - 1. + zeta <= 0.+eps) &&
            ( eta - 1. + zeta <= 0.+eps) &&
            ( -xi - 1. + zeta <= 0.+eps) &&
            (            zeta >= 0.-eps))
          return true;

        return false;
      }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
      {
        // The reference infhex8 is a [-1,1]^3.
        if ((xi   >= -1.-eps) &&
            (xi   <=  1.+eps) &&
            (eta  >= -1.-eps) &&
            (eta  <=  1.+eps) &&
            (zeta >= -1.-eps) &&
            (zeta <=  1.+eps))
          {
            return true;
          }
        return false;
      }

    case INFPRISM6:
    case INFPRISM12:
      {
        // inside the reference triangle with zeta in [-1,1]
        if ((xi   >=  0.-eps) &&
            (eta  >=  0.-eps) &&
            (zeta >= -1.-eps) &&
            (zeta <=  1.+eps) &&
            ((xi + eta) <= 1.+eps))
          {
            return true;
          }

        return false;
      }
#endif

    default:
      libmesh_error_msg("ERROR: Unknown element type " << Utility::enum_to_string(t));
    }

  // If we get here then the point is _not_ in the
  // reference element.   Better return false.

  return false;
}
#endif // LIBMESH_ENABLE_DEPRECATED



void FEAbstract::print_JxW(std::ostream & os) const
{
  this->_fe_map->print_JxW(os);
}



void FEAbstract::print_xyz(std::ostream & os) const
{
  this->_fe_map->print_xyz(os);
}


void FEAbstract::print_info(std::ostream & os) const
{
  os << "phi[i][j]: Shape function i at quadrature pt. j" << std::endl;
  this->print_phi(os);

  os << "dphi[i][j]: Shape function i's gradient at quadrature pt. j" << std::endl;
  this->print_dphi(os);

  os << "XYZ locations of the quadrature pts." << std::endl;
  this->print_xyz(os);

  os << "Values of JxW at the quadrature pts." << std::endl;
  this->print_JxW(os);
}


std::ostream & operator << (std::ostream & os, const FEAbstract & fe)
{
  fe.print_info(os);
  return os;
}



#ifdef LIBMESH_ENABLE_AMR

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
void FEAbstract::compute_node_constraints (NodeConstraints & constraints,
                                           const Elem * elem)
{
  libmesh_assert(elem);

  const unsigned int Dim = elem->dim();

  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  // Only constrain active and ancestor elements
  if (elem->subactive())
    return;


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (elem->infinite())
    {
      const FEType fe_t(elem->default_order(), FEMap::map_fe_type(*elem));

      // expand the infinite_compute_constraint in its template-arguments.
      switch(Dim)
      {
         case 2:
            {
            inf_fe_family_mapping_switch(2, inf_compute_node_constraints (constraints, elem) , ,; break;);
            break;
          }
         case 3:
            {
            inf_fe_family_mapping_switch(3, inf_compute_node_constraints (constraints, elem) , ,; break;);
            break;
            }
         default:
           libmesh_error_msg("Invalid dim = " << Dim);
      }
      return;
    }

#endif
  const FEFamily mapping_family = FEMap::map_fe_type(*elem);
  const FEType fe_type(elem->default_side_order(), mapping_family);

  // Pull objects out of the loop to reduce heap operations
  std::vector<const Node *> my_nodes, parent_nodes;
  std::unique_ptr<const Elem> my_side, parent_side;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (auto s : elem->side_index_range())
    if (elem->neighbor_ptr(s) != nullptr &&
        elem->neighbor_ptr(s) != remote_elem)
      if (elem->neighbor_ptr(s)->level() < elem->level()) // constrain dofs shared between
        {                                                 // this element and ones coarser
          // than this element.
          // Get pointers to the elements of interest and its parent.
          const Elem * parent = elem->parent();

          // This can't happen...  Only level-0 elements have nullptr
          // parents, and no level-0 elements can be at a higher
          // level than their neighbors!
          libmesh_assert(parent);

          elem->build_side_ptr(my_side, s);
          parent->build_side_ptr(parent_side, s);

          const unsigned int n_side_nodes = my_side->n_nodes();

          my_nodes.clear();
          my_nodes.reserve (n_side_nodes);
          parent_nodes.clear();
          parent_nodes.reserve (n_side_nodes);

          for (unsigned int n=0; n != n_side_nodes; ++n)
            my_nodes.push_back(my_side->node_ptr(n));

          for (unsigned int n=0; n != n_side_nodes; ++n)
            parent_nodes.push_back(parent_side->node_ptr(n));

          for (unsigned int my_side_n=0;
               my_side_n < n_side_nodes;
               my_side_n++)
            {
              // We can have an FE type that supports an order
              // partially, such that sides do not support the same
              // order.  E.g. we say that a LAGRANGE PRISM21 supports
              // "third" order to distinguish its shape functions from
              // a PRISM18, but the QUAD9 sides will still only
              // support second order.
              FEType side_fe_type = fe_type;
              const int side_max_order =
                FEInterface::max_order(fe_type, my_side->type());

              if ((int)fe_type.order > side_max_order)
                side_fe_type.order = side_max_order;

              // Do not use the p_level(), if any, that is inherited by the side.
              libmesh_assert_less
                (my_side_n,
                 FEInterface::n_dofs(side_fe_type, /*extra_order=*/0,
                                     my_side.get()));

              const Node * my_node = my_nodes[my_side_n];

              // The support point of the DOF
              const Point & support_point = *my_node;

              // Figure out where my node lies on their reference element.
              const Point mapped_point = FEMap::inverse_map(Dim-1,
                                                            parent_side.get(),
                                                            support_point);

              // Compute the parent's side shape function values.
              for (unsigned int their_side_n=0;
                   their_side_n < n_side_nodes;
                   their_side_n++)
                {
                  // Do not use the p_level(), if any, that is inherited by the side.
                  libmesh_assert_less
                    (their_side_n,
                     FEInterface::n_dofs(side_fe_type,
                                         /*extra_order=*/0,
                                         parent_side.get()));

                  const Node * their_node = parent_nodes[their_side_n];
                  libmesh_assert(their_node);

                  // Do not use the p_level(), if any, that is inherited by the side.
                  const Real their_value = FEInterface::shape(side_fe_type,
                                                              /*extra_order=*/0,
                                                              parent_side.get(),
                                                              their_side_n,
                                                              mapped_point);

                  const Real their_mag = std::abs(their_value);
#ifdef DEBUG
                  // Protect for the case u_i ~= u_j,
                  // in which case i better equal j.
                  if (their_mag > 0.999)
                    {
                      libmesh_assert_equal_to (my_node, their_node);
                      libmesh_assert_less (std::abs(their_value - 1.), 0.001);
                    }
                  else
#endif
                    // To make nodal constraints useful for constructing
                    // sparsity patterns faster, we need to get EVERY
                    // POSSIBLE constraint coupling identified, even if
                    // there is no coupling in the isoparametric
                    // Lagrange case.
                    if (their_mag < 1.e-5)
                      {
                        // since we may be running this method concurrently
                        // on multiple threads we need to acquire a lock
                        // before modifying the shared constraint_row object.
                        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                        // A reference to the constraint row.
                        NodeConstraintRow & constraint_row = constraints[my_node].first;

                        constraint_row.emplace(their_node, 0.);
                      }
                  // To get nodal coordinate constraints right, only
                  // add non-zero and non-identity values for Lagrange
                  // basis functions.
                    else // (1.e-5 <= their_mag <= .999)
                      {
                        // since we may be running this method concurrently
                        // on multiple threads we need to acquire a lock
                        // before modifying the shared constraint_row object.
                        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                        // A reference to the constraint row.
                        NodeConstraintRow & constraint_row = constraints[my_node].first;

                        constraint_row.emplace(their_node, their_value);
                      }
                }
            }
        }
}

#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif // #ifdef LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_PERIODIC

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
void FEAbstract::compute_periodic_node_constraints (NodeConstraints & constraints,
                                                    const PeriodicBoundaries & boundaries,
                                                    const MeshBase & mesh,
                                                    const PointLocatorBase * point_locator,
                                                    const Elem * elem)
{
  // Only bother if we truly have periodic boundaries
  if (boundaries.empty())
    return;

  libmesh_assert(elem);

  // Only constrain active elements with this method
  if (!elem->active())
    return;

  const unsigned int Dim = elem->dim();

  const FEFamily mapping_family = FEMap::map_fe_type(*elem);
  const FEType fe_type(elem->default_side_order(), mapping_family);

  // Pull objects out of the loop to reduce heap operations
  std::vector<const Node *> my_nodes, neigh_nodes;
  std::unique_ptr<const Elem> my_side, neigh_side;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  std::vector<boundary_id_type> bc_ids;
  for (auto s : elem->side_index_range())
    {
      if (elem->neighbor_ptr(s))
        continue;

      mesh.get_boundary_info().boundary_ids (elem, s, bc_ids);
      for (const auto & boundary_id : bc_ids)
        {
          const PeriodicBoundaryBase * periodic = boundaries.boundary(boundary_id);
          if (periodic)
            {
              libmesh_assert(point_locator);

              // Get pointers to the element's neighbor.
              unsigned int s_neigh;
              const Elem * neigh = boundaries.neighbor(boundary_id, *point_locator, elem, s, &s_neigh);

              libmesh_error_msg_if
                (!neigh, "PeriodicBoundaries can't find a periodic neighbor for element " <<
                         elem->id() << " side " << s);

              // h refinement constraints:
              // constrain dofs shared between
              // this element and ones as coarse
              // as or coarser than this element.
              if (neigh->level() <= elem->level())
                {
#ifdef LIBMESH_ENABLE_AMR
                  libmesh_assert(neigh->active());
#endif // #ifdef LIBMESH_ENABLE_AMR

                  elem->build_side_ptr(my_side, s);
                  neigh->build_side_ptr(neigh_side, s_neigh);

                  const unsigned int n_side_nodes = my_side->n_nodes();

                  my_nodes.clear();
                  my_nodes.reserve (n_side_nodes);
                  neigh_nodes.clear();
                  neigh_nodes.reserve (n_side_nodes);

                  for (unsigned int n=0; n != n_side_nodes; ++n)
                    my_nodes.push_back(my_side->node_ptr(n));

                  for (unsigned int n=0; n != n_side_nodes; ++n)
                    neigh_nodes.push_back(neigh_side->node_ptr(n));

                  // Make sure we're not adding recursive constraints
                  // due to the redundancy in the way we add periodic
                  // boundary constraints, or adding constraints to
                  // nodes that already have AMR constraints
                  std::vector<bool> skip_constraint(n_side_nodes, false);

                  for (unsigned int my_side_n=0;
                       my_side_n < n_side_nodes;
                       my_side_n++)
                    {
                      // Do not use the p_level(), if any, that is inherited by the side.
                      libmesh_assert_less (my_side_n, FEInterface::n_dofs(fe_type, /*extra_order=*/0, my_side.get()));

                      const Node * my_node = my_nodes[my_side_n];

                      // If we've already got a constraint on this
                      // node, then the periodic constraint is
                      // redundant
                      {
                        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                        if (constraints.count(my_node))
                          {
                            skip_constraint[my_side_n] = true;
                            continue;
                          }
                      }

                      // Compute the neighbors's side shape function values.
                      for (unsigned int their_side_n=0;
                           their_side_n < n_side_nodes;
                           their_side_n++)
                        {
                          // Do not use the p_level(), if any, that is inherited by the side.
                          libmesh_assert_less (their_side_n, FEInterface::n_dofs(fe_type, /*extra_order=*/0, neigh_side.get()));

                          const Node * their_node = neigh_nodes[their_side_n];

                          // If there's a constraint on an opposing node,
                          // we need to see if it's constrained by
                          // *our side* making any periodic constraint
                          // on us recursive
                          {
                            Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                            if (!constraints.count(their_node))
                              continue;

                            const NodeConstraintRow & their_constraint_row =
                              constraints[their_node].first;

                            for (unsigned int orig_side_n=0;
                                 orig_side_n < n_side_nodes;
                                 orig_side_n++)
                              {
                                // Do not use the p_level(), if any, that is inherited by the side.
                                libmesh_assert_less (orig_side_n, FEInterface::n_dofs(fe_type, /*extra_order=*/0, my_side.get()));

                                const Node * orig_node = my_nodes[orig_side_n];

                                if (their_constraint_row.count(orig_node))
                                  skip_constraint[orig_side_n] = true;
                              }
                          }
                        }
                    }
                  for (unsigned int my_side_n=0;
                       my_side_n < n_side_nodes;
                       my_side_n++)
                    {
                      // Do not use the p_level(), if any, that is inherited by the side.
                      libmesh_assert_less (my_side_n, FEInterface::n_dofs(fe_type, /*extra_order=*/0, my_side.get()));

                      if (skip_constraint[my_side_n])
                        continue;

                      const Node * my_node = my_nodes[my_side_n];

                      // Figure out where my node lies on their reference element.
                      const Point neigh_point = periodic->get_corresponding_pos(*my_node);

                      // Figure out where my node lies on their reference element.
                      const Point mapped_point =
                        FEMap::inverse_map(Dim-1, neigh_side.get(),
                                           neigh_point);

                      for (unsigned int their_side_n=0;
                           their_side_n < n_side_nodes;
                           their_side_n++)
                        {
                          // Do not use the p_level(), if any, that is inherited by the side.
                          libmesh_assert_less (their_side_n, FEInterface::n_dofs(fe_type, /*extra_order=*/0, neigh_side.get()));

                          const Node * their_node = neigh_nodes[their_side_n];
                          libmesh_assert(their_node);

                          // Do not use the p_level(), if any, that is inherited by the side.
                          const Real their_value = FEInterface::shape(fe_type,
                                                                      /*extra_order=*/0,
                                                                      neigh_side.get(),
                                                                      their_side_n,
                                                                      mapped_point);

                          // since we may be running this method concurrently
                          // on multiple threads we need to acquire a lock
                          // before modifying the shared constraint_row object.
                          {
                            Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                            NodeConstraintRow & constraint_row =
                              constraints[my_node].first;

                            constraint_row.emplace(their_node, their_value);
                          }
                        }
                    }
                }
            }
        }
    }
}
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif // LIBMESH_ENABLE_PERIODIC


unsigned int FEAbstract::n_quadrature_points () const
{
  if (this->shapes_on_quadrature)
    {
      libmesh_assert(this->qrule);
      libmesh_assert_equal_to(this->qrule->n_points(),
                              this->_n_total_qp);
    }
  return this->_n_total_qp;
}

} // namespace libMesh
