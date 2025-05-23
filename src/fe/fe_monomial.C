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
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"

namespace libMesh
{

unsigned int monomial_n_dofs(const Elem * e, const Order o)
{
  libmesh_assert(e);
  return monomial_n_dofs(e->type(), o);
}


unsigned int monomial_n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {

      // constant shape functions
      // no matter what shape there is only one DOF.
    case CONSTANT:
      return (t != INVALID_ELEM) ? 1 : 0;


      // Discontinuous linear shape functions
      // expressed in the monomials.
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 2;

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 3;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
          case C0POLYHEDRON:
            return 4;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous quadratic shape functions
      // expressed in the monomials.
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 3;

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 6;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
          case C0POLYHEDRON:
            return 10;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous cubic shape functions
      // expressed in the monomials.
    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 4;

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 10;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
          case C0POLYHEDRON:
            return 20;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous quartic shape functions
      // expressed in the monomials.
    case FOURTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 5;

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 15;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case C0POLYHEDRON:
            return 35;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return (order+1);

          case C0POLYGON:
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case TRI7:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return (order+1)*(order+2)/2;

          case TET4:
          case TET10:
          case TET14:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case C0POLYHEDRON:
            return (order+1)*(order+2)*(order+3)/6;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
    }
} // monomial_n_dofs()


// Anonymous namespace for local helper functions
namespace {

void monomial_nodal_soln(const Elem * elem,
                         const Order order,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> & nodal_soln,
                         const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = order + add_p_level*elem->p_level();

  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
        libmesh_assert_equal_to (elem_soln.size(), 1);

        std::fill(nodal_soln.begin(), nodal_soln.end(), elem_soln[0]);

        return;
      }


      // For other orders, do interpolation at the nodes
      // explicitly.
    default:
      {
        // FEType object to be passed to various FEInterface functions below.
        FEType fe_type(order, MONOMIAL);

        const unsigned int n_sf =
          FEInterface::n_shape_functions(fe_type, elem);

        std::vector<Point> refspace_nodes;
        FEBase::get_refspace_nodes(elem_type,refspace_nodes);
        libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);
        libmesh_assert_equal_to (elem_soln.size(), n_sf);

        // Zero before summation
        std::fill(nodal_soln.begin(), nodal_soln.end(), 0);

        for (unsigned int n=0; n<n_nodes; n++)
          // u_i = Sum (alpha_i phi_i)
          for (unsigned int i=0; i<n_sf; i++)
            nodal_soln[n] += elem_soln[i] *
              FEInterface::shape(fe_type, elem, i, refspace_nodes[n]);

        return;
      } // default
    } // switch
} // monomial_nodal_soln()




} // anonymous namespace


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(MONOMIAL, monomial_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(MONOMIAL)


// Full specialization of n_dofs() function for every dimension
template <> unsigned int FE<0,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }

template <> unsigned int FE<0,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
// Monomials have no dofs at nodes, only element dofs.
template <> unsigned int FE<0,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <> unsigned int FE<0,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }

template <> unsigned int FE<0,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<1,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<2,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<3,MONOMIAL>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }


// Full specialization of get_continuity() function for every dimension.
template <> FEContinuity FE<0,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }

// Full specialization of is_hierarchic() function for every dimension.
// The monomials are hierarchic!
template <> bool FE<0,MONOMIAL>::is_hierarchic() const { return true; }
template <> bool FE<1,MONOMIAL>::is_hierarchic() const { return true; }
template <> bool FE<2,MONOMIAL>::is_hierarchic() const { return true; }
template <> bool FE<3,MONOMIAL>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR

// Full specialization of compute_constraints() function for 2D and
// 3D only.  There are no constraints for discontinuous elements, so
// we do nothing.
template <> void FE<2,MONOMIAL>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem *) {}
template <> void FE<3,MONOMIAL>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem *) {}

#endif // #ifdef LIBMESH_ENABLE_AMR

// Full specialization of shapes_need_reinit() function for every dimension.
template <> bool FE<0,MONOMIAL>::shapes_need_reinit() const { return false; }
template <> bool FE<1,MONOMIAL>::shapes_need_reinit() const { return false; }
template <> bool FE<2,MONOMIAL>::shapes_need_reinit() const { return false; }
template <> bool FE<3,MONOMIAL>::shapes_need_reinit() const { return false; }

} // namespace libMesh
