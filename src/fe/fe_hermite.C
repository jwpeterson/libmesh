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
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"

namespace libMesh
{

// ------------------------------------------------------------
// Hermite-specific implementations

// Anonymous namespace for local helper functions
namespace {

void hermite_nodal_soln(const Elem * elem,
                        const Order order,
                        const std::vector<Number> & elem_soln,
                        std::vector<Number> & nodal_soln,
                        const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(order, HERMITE);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(fe_type, elem, add_p_level);

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
        FEInterface::shape(fe_type, elem, i, refspace_nodes[n], add_p_level);

} // hermite_nodal_soln()



unsigned int HERMITE_n_dofs(const ElemType t, const Order o)
{
#ifdef DEBUG
  libmesh_error_msg_if(o < 3, "Error: Hermite elements require order>=3, but you asked for order=" << o);
#endif

  // Piecewise (bi/tri)cubic C1 Hermite splines
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case EDGE3:
      return (o+1);

    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case QUAD9:
    case QUADSHELL9:
      return ((o+1)*(o+1));

    case HEX8:
    case HEX20:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case HEX27:
      return ((o+1)*(o+1)*(o+1));

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }
} // HERMITE_n_dofs()



unsigned int HERMITE_n_dofs(const Elem * e, const Order o)
{
  libmesh_assert(e);
  return HERMITE_n_dofs(e->type(), o);
}



unsigned int HERMITE_n_dofs_at_node(const ElemType t,
                                    const Order o,
                                    const unsigned int n)
{
  libmesh_assert_greater (o, 2);
  // Piecewise (bi/tri)cubic C1 Hermite splines
  switch (t)
    {
    case NODEELEM:
      return 1;
    case EDGE2:
    case EDGE3:
      {
        switch (n)
          {
          case 0:
          case 1:
            return 2;
          case 2:
            //          Interior DoFs are carried on Elems
            //    return (o-3);
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
          }
      }

    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      {
        switch (n)
          {
            // Vertices
          case 0:
          case 1:
          case 2:
          case 3:
            return 4;
            // Edges
          case 4:
          case 5:
          case 6:
          case 7:
            return (2*(o-3));
          case 8:
            //          Interior DoFs are carried on Elems
            //    return ((o-3)*(o-3));
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD4/8/9!");
          }
      }

    case HEX8:
    case HEX20:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case HEX27:
      {
        switch (n)
          {
            // Vertices
          case 0:
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
          case 6:
          case 7:
            return 8;
            // Edges
          case 8:
          case 9:
          case 10:
          case 11:
          case 12:
          case 13:
          case 14:
          case 15:
          case 16:
          case 17:
          case 18:
          case 19:
            return (4*(o-3));
            // Faces
          case 20:
          case 21:
          case 22:
          case 23:
          case 24:
          case 25:
            return (2*(o-3)*(o-3));
          case 26:
            // Interior DoFs are carried on Elems
            //    return ((o-3)*(o-3)*(o-3));
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for HEX8/20/27!");
          }
      }

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }
} // HERMITE_n_dofs_at_node()



unsigned int HERMITE_n_dofs_at_node(const Elem & e,
                                    const Order o,
                                    const unsigned int n)
{
  return HERMITE_n_dofs_at_node(e.type(), o, n);
}



unsigned int HERMITE_n_dofs_per_elem(const ElemType t,
                                     const Order o)
{
  libmesh_assert_greater (o, 2);

  switch (t)
    {
    case NODEELEM:
      return 0;
    case EDGE2:
    case EDGE3:
      return (o-3);
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      return ((o-3)*(o-3));
    case HEX8:
      libmesh_assert_less (o, 4);
      libmesh_fallthrough();
    case HEX20:
    case HEX27:
      return ((o-3)*(o-3)*(o-3));

    case INVALID_ELEM:
      return 0;

    default:
      libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
    }
} // HERMITE_n_dofs_per_elem()



unsigned int HERMITE_n_dofs_per_elem(const Elem & e,
                                     const Order o)
{
  return HERMITE_n_dofs_per_elem(e.type(), o);
}


} // anonymous namespace


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(HERMITE, hermite_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(HERMITE)


// Instantiate n_dofs*() functions for every dimension
LIBMESH_DEFAULT_NDOFS(HERMITE)


// Hermite FEMs are C^1 continuous
template <> FEContinuity FE<0,HERMITE>::get_continuity() const { return C_ONE; }
template <> FEContinuity FE<1,HERMITE>::get_continuity() const { return C_ONE; }
template <> FEContinuity FE<2,HERMITE>::get_continuity() const { return C_ONE; }
template <> FEContinuity FE<3,HERMITE>::get_continuity() const { return C_ONE; }

// Hermite FEMs are hierarchic
template <> bool FE<0,HERMITE>::is_hierarchic() const { return true; }
template <> bool FE<1,HERMITE>::is_hierarchic() const { return true; }
template <> bool FE<2,HERMITE>::is_hierarchic() const { return true; }
template <> bool FE<3,HERMITE>::is_hierarchic() const { return true; }


#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <>
void FE<2,HERMITE>::compute_constraints (DofConstraints & constraints,
                                         DofMap & dof_map,
                                         const unsigned int variable_number,
                                         const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <>
void FE<3,HERMITE>::compute_constraints (DofConstraints & constraints,
                                         DofMap & dof_map,
                                         const unsigned int variable_number,
                                         const Elem * elem)
{ compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Hermite FEM shapes need reinit
template <> bool FE<0,HERMITE>::shapes_need_reinit() const { return true; }
template <> bool FE<1,HERMITE>::shapes_need_reinit() const { return true; }
template <> bool FE<2,HERMITE>::shapes_need_reinit() const { return true; }
template <> bool FE<3,HERMITE>::shapes_need_reinit() const { return true; }

} // namespace libMesh
