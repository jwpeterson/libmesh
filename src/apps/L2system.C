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

#include "L2system.h"

#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

using namespace libMesh;

HilbertSystem::~HilbertSystem () = default;

void HilbertSystem::init_data ()
{
  this->add_variable ("u", static_cast<Order>(_fe_order),
                      Utility::string_to_enum<FEFamily>(_fe_family));

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();
}



void HilbertSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * my_fe = nullptr;

  // Now make sure we have requested all the data
  // we need to build the L2 system.

  // We might have a multi-dimensional mesh
  const std::set<unsigned char> & elem_dims =
    c.elem_dimensions();

  for (const auto & dim : elem_dims)
    {
      c.get_element_fe( 0, my_fe, dim );

      my_fe->get_JxW();
      my_fe->get_phi();
      my_fe->get_xyz();

      if (this->_hilbert_order > 0)
        my_fe->get_dphi();

      c.get_side_fe( 0, my_fe, dim );
      my_fe->get_nothing();
    }

  // Build a corresponding context for the input system if we haven't
  // already
  auto & input_context = input_contexts[&c];
  if (input_system && !input_context)
    {
      input_context = std::make_unique<FEMContext>(*input_system);

      libmesh_assert(_goal_func.get());
      _goal_func->init_context(*input_context);
    }

  FEMSystem::init_context(context);
}


bool HilbertSystem::element_time_derivative (bool request_jacobian,
                                             DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  const Elem & elem = c.get_elem();

  if (!_subdomains_list.empty() &&
      !_subdomains_list.count(elem.subdomain_id()))
    return request_jacobian;

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = c.get_element_fe(0)->get_JxW();

  const std::vector<std::vector<Real>> & phi = c.get_element_fe(0)->get_phi();

  const std::vector<Point> & xyz = c.get_element_fe(0)->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.n_dof_indices(0);

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> & K = c.get_elem_jacobian(0, 0);
  DenseSubVector<Number> & F = c.get_elem_residual(0);

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  FEMContext & input_c = *libmesh_map_find(input_contexts, &c);
  if (input_system)
    {
      input_c.pre_fe_reinit(*input_system, &elem);
      input_c.elem_fe_reinit();
    }

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const Number u = c.interior_value(0, qp);
      const Number ufunc = (*_goal_func)(input_c, xyz[qp]);
      const Number err_u = u - ufunc;

      for (unsigned int i=0; i != n_u_dofs; i++)
        F(i) += JxW[qp] * (err_u * phi[i][qp]);

      if (_hilbert_order > 0)
        {
          const std::vector<std::vector<RealGradient>> & dphi =
            c.get_element_fe(0)->get_dphi();

          const Gradient grad_u = c.interior_gradient(0, qp);
          Gradient ufuncgrad = (*_goal_grad)(input_c, xyz[qp]);
          const Gradient err_grad_u = grad_u - ufuncgrad;

          for (unsigned int i=0; i != n_u_dofs; i++)
            F(i) += JxW[qp] * (err_grad_u * dphi[i][qp]);
        }

      if (request_jacobian)
        {
          const Number JxWxD = JxW[qp] *
            context.get_elem_solution_derivative();

          for (unsigned int i=0; i != n_u_dofs; i++)
            for (unsigned int j=0; j != n_u_dofs; ++j)
              K(i,j) += JxWxD * (phi[i][qp] * phi[j][qp]);

          if (_hilbert_order > 0)
            {
              const std::vector<std::vector<RealGradient>> & dphi =
                c.get_element_fe(0)->get_dphi();

              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; ++j)
                  K(i,j) += JxWxD * (dphi[i][qp] * dphi[j][qp]);
            }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}
