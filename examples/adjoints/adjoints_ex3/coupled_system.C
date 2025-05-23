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


// Example includes
#include "coupled_system.h"

// libMesh includes
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/parallel.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

// C++ includes
#include <functional> // std::reference_wrapper
#include <memory>

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function to set the Dirichlet boundary function for the inlet boundary velocities
class BdyFunction : public FunctionBase<Number>
{
public:
  BdyFunction (unsigned int u_var,
               unsigned int v_var,
               int sign) :
    _u_var(u_var),
    _v_var(v_var),
    _sign(sign)
  { this->_initialized = true; }

  virtual Number operator() (const Point &,
                             const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Number> & output)
  {
    output.resize(2);
    output.zero();
    const Real y=p(1);
    // Set the parabolic inflow boundary conditions at stations 0 & 1
    output(_u_var) = (_sign)*((y-2) * (y-3));
    output(_v_var) = 0;
  }

  virtual std::unique_ptr<FunctionBase<Number>> clone() const
  { return std::make_unique<BdyFunction>(_u_var, _v_var, int(_sign)); }

private:
  const unsigned int _u_var, _v_var;
  const Real _sign;
};


void CoupledSystem::init_data ()
{
  // Check the input file for Reynolds number, application type,
  // approximation type
  GetPot infile("coupled_system.in");
  Peclet = infile("Peclet", 1.);
  unsigned int pressure_p = infile("pressure_p", 1);
  std::string fe_family = infile("fe_family", std::string("LAGRANGE"));

  // LBB needs better-than-quadratic velocities for better-than-linear
  // pressures, and libMesh needs non-Lagrange elements for
  // better-than-quadratic velocities.
  libmesh_assert((pressure_p == 1) || (fe_family != "LAGRANGE"));

  FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);

  // Add the velocity components "u" & "v".  They
  // will be approximated using second-order approximation.
  u_var = this->add_variable ("u", static_cast<Order>(pressure_p+1),
                              fefamily);
  v_var = this->add_variable ("v", static_cast<Order>(pressure_p+1),
                              fefamily);

  // Add the pressure variable "p". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  p_var = this->add_variable ("p", static_cast<Order>(pressure_p),
                              fefamily);

  // Add the Concentration variable "C". They will
  // be approximated using second-order approximation, the same as the velocity components
  C_var = this->add_variable ("C", static_cast<Order>(pressure_p+1),
                              fefamily);

  // Tell the system to march velocity and concentration forward in
  // time, with first order time derivatives, but leave pressure as a
  // constraint only
  this->time_evolving(u_var, 1);
  this->time_evolving(v_var, 1);
  this->time_evolving(C_var, 1);

  // Useful debugging options
  this->verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);

  // Set Dirichlet boundary conditions
  const boundary_id_type left_inlet_id = 0;

  const boundary_id_type right_inlet_id = 1;
  std::set<boundary_id_type> right_inlet_bdy;
  right_inlet_bdy.insert(right_inlet_id);

  const boundary_id_type outlets_id = 2;
  std::set<boundary_id_type> outlets_bdy;
  outlets_bdy.insert(outlets_id);

  const boundary_id_type wall_id = 3;

  // The zero and constant functions
  ZeroFunction<Number> zero;
  ConstFunction<Number> one(1);

  // We need two boundary functions for the inlets, because the signs on the velocities
  // will be different
  int velocity_sign = 1;
  BdyFunction inflow_left(u_var, v_var, -velocity_sign);
  BdyFunction inflow_right(u_var, v_var, velocity_sign);

#ifdef LIBMESH_ENABLE_DIRICHLET
  // On the walls we will apply the no slip and no penetration boundary condition, u=0, v=0
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary ({wall_id}, {u_var,v_var}, &zero));

  // On the inlet (left), we apply parabolic inflow boundary conditions for the velocity, u = - (y-2)*(y-3), v=0
  // and set C = 1
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary ({left_inlet_id}, {u_var,v_var}, &inflow_left));
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary ({left_inlet_id}, {C_var}, &one));

  // On the inlet (right), we apply parabolic inflow boundary conditions for the velocity, u = (y-2)*(y-3), v=0
  // and set C = 0
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (right_inlet_bdy, {u_var,v_var}, &inflow_right));
  this->get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (right_inlet_bdy, {C_var}, &zero));
#endif // LIBMESH_ENABLE_DIRICHLET

  // Note that the remaining boundary conditions are the natural boundary conditions for the concentration
  // on the wall (grad(c) dot n = 0) and natural boundary conditions for the velocity and the concentration
  // on the outlets ((grad(velocity) dot n - p n) dot t = 0, grad(C) dot n = 0)

  // Do the parent's initialization after variables and boundary constraints are defined
  FEMSystem::init_data();
}

void CoupledSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // We should prerequest all the data
  // we will need to build the linear system.
  // Note that the concentration and velocity components
  // use the same basis.
  FEBase * u_elem_fe = nullptr;
  c.get_element_fe(u_var, u_elem_fe);
  u_elem_fe->get_JxW();
  u_elem_fe->get_phi();
  u_elem_fe->get_dphi();
  u_elem_fe->get_xyz();

  FEBase * p_elem_fe = nullptr;
  c.get_element_fe(p_var, p_elem_fe);
  p_elem_fe->get_phi();

  FEBase * side_fe = nullptr;
  c.get_side_fe(u_var, side_fe);

  side_fe->get_JxW();
  side_fe->get_phi();
  side_fe->get_xyz();

  // Don't waste time on side computations for p
  FEBase * p_side_fe = nullptr;
  c.get_side_fe(p_var, p_side_fe);
  p_side_fe->get_nothing();
}


bool CoupledSystem::element_time_derivative (bool request_jacobian,
                                             DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * u_elem_fe = nullptr;
  c.get_element_fe(u_var, u_elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = u_elem_fe->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real>> & phi = u_elem_fe->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = u_elem_fe->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  FEBase * p_elem_fe = nullptr;
  c.get_element_fe(p_var, p_elem_fe);

  const std::vector<std::vector<Real>> & psi = p_elem_fe->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_p_dofs = c.n_dof_indices(p_var);
  const unsigned int n_u_dofs = c.n_dof_indices(u_var);
  libmesh_assert_equal_to (n_u_dofs, c.n_dof_indices(v_var));

  // Stokes equation stiffness matrix diagonal blocks. There is no
  // uv/vu coupling for linear Stokes.
  std::reference_wrapper<DenseSubMatrix<Number>> Kdiag[2] =
    {
      c.get_elem_jacobian(u_var, u_var),
      c.get_elem_jacobian(v_var, v_var)
    };

  // Stokes equation velocity-pressure coupling blocks
  // Kup
  // Kvp
  std::reference_wrapper<DenseSubMatrix<Number>> B[2] =
    {
      c.get_elem_jacobian(u_var, p_var),
      c.get_elem_jacobian(v_var, p_var),
    };

  // Stokes equation residuals
  // Fu
  // Fv
  std::reference_wrapper<DenseSubVector<Number>> F[2] =
    {
      c.get_elem_residual(u_var),
      c.get_elem_residual(v_var)
    };

  // Concentration/velocity coupling blocks
  // KCu
  // KCv
  std::reference_wrapper<DenseSubMatrix<Number>> C[2] =
    {
      c.get_elem_jacobian(C_var, u_var),
      c.get_elem_jacobian(C_var, v_var)
    };

  // Concentration diagonal block and residual
  DenseSubMatrix<Number> & KCC = c.get_elem_jacobian(C_var, C_var);
  DenseSubVector<Number> & FC = c.get_elem_residual(C_var);

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Store (grad_u, grad_v) at current Newton iterate
  Gradient grad_uv[2];

  // Velocity at current Newton iterate
  NumberVectorValue U;

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the current Newton iterate
      Number p = c.interior_value(p_var, qp);
      Gradient grad_C = c.interior_gradient(C_var, qp);

      c.interior_value(u_var, qp, U(0));
      c.interior_value(v_var, qp, U(1));
      c.interior_gradient(u_var, qp, grad_uv[0]);
      c.interior_gradient(v_var, qp, grad_uv[1]);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          // Stokes equations residuals
          for (unsigned int d=0; d<2; ++d)
            F[d](i) += JxW[qp] *
              (p*dphi[i][qp](d) -                // pressure term
               (grad_uv[d]*dphi[i][qp]));        // diffusion term

          // Concentration Equation Residual
          FC(i) += JxW[qp] *
            ((U*grad_C)*phi[i][qp] +                // convection term
             (1./Peclet)*(grad_C*dphi[i][qp]));     // diffusion term

          // Note that the Fp block is identically zero unless we are using
          // some kind of artificial compressibility scheme...

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  for (unsigned int d=0; d<2; ++d)
                    {
                      // Matrix contributions for the uu and vv couplings.
                      Kdiag[d](i,j) += JxW[qp] * (-(dphi[i][qp]*dphi[j][qp]));

                      // Matrix contributions for the Cu and Cv couplings.
                      C[d](i,j) += JxW[qp] * ((phi[j][qp]*grad_C(d))*phi[i][qp]);
                    }

                  // Concentration equation Jacobian contribution
                  KCC(i,j) += JxW[qp]*
                    ((U*dphi[j][qp])*phi[i][qp] +            // nonlinear term (convection)
                     (1./Peclet)*(dphi[j][qp]*dphi[i][qp])); // diffusion term
                }

              // Matrix contributions for the up and vp couplings.
              for (unsigned int j=0; j != n_p_dofs; j++)
                for (unsigned int d=0; d<2; ++d)
                  B[d](i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](d);
            }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}


bool CoupledSystem::element_constraint (bool request_jacobian,
                                        DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * u_elem_fe = nullptr;
  c.get_element_fe(u_var, u_elem_fe);

  FEBase * p_elem_fe = nullptr;
  c.get_element_fe(p_var, p_elem_fe);

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> & JxW = u_elem_fe->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = u_elem_fe->get_dphi();

  // The pressure shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real>> & psi = p_elem_fe->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.n_dof_indices(u_var);
  const unsigned int n_p_dofs = c.n_dof_indices(p_var);

  // pressure-velocity coupling blocks
  // Kpu Kpv
  std::reference_wrapper<DenseSubMatrix<Number>> B[2] =
    {
      c.get_elem_jacobian(p_var, u_var),
      c.get_elem_jacobian(p_var, v_var)
    };

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Number> & Fp = c.get_elem_residual(p_var);

  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Store (grad_u, grad_v) at current Newton iterate
  Gradient grad_uv[2];

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate
      c.interior_gradient(u_var, qp, grad_uv[0]),
      c.interior_gradient(v_var, qp, grad_uv[1]);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          for (unsigned int d=0; d<2; ++d)
            Fp(i) += JxW[qp] * psi[i][qp] * grad_uv[d](d);

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                for (unsigned int d=0; d<2; ++d)
                  B[d](i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](d);
            }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

void CoupledSystem::postprocess()
{
  // We need to overload the postprocess function to set the
  // computed_QoI variable of the CoupledSystem class to the qoi value
  // stored in System::qoi[0]
  computed_QoI = this->get_qoi_value(0);
}

// These functions supply the nonlinear weighting for the adjoint residual error estimate which
// arise due to the convection term in the convection-diffusion equation:
// ||e((u_1)_h C,1)||_{L2} ||e(C^*)||_{L2} + ||e(u1 C,1_h)||_{L2} ||e(C^*)||_{L2}
// ||e((u_2)_h C,2)||_{L2} ||e(C^*)||_{L2} + ||e(u2 C,2_h)||_{L2} ||e(C^*)||_{L2}
// These functions compute (u_1)_h or C,1_h , and (u_2)_h or C,2_h , and supply it to the weighted patch recovery error estimator
// In CoupledFEMFunctionsx, the object built with var = 0, returns the (u_1)_h weight, while
// the object built with var = 1, returns the C,1_h weight. The switch statement
// distinguishes the behavior of the two objects
// Same thing for CoupledFEMFunctionsy
Number CoupledFEMFunctionsx::operator()(const FEMContext & c,
                                        const Point & p,
                                        const Real /* time */)
{
  Number weight = 0.0;

  switch(var)
    {
    case 0:
      {
        Gradient grad_C = c.point_gradient(3, p);

        weight = grad_C(0);
      }
      break;

    case 3:
      {
        Number u = c.point_value(0, p);

        weight = u;
      }
      break;

    default:
      libmesh_error_msg("Wrong variable number "                        \
                        << var                                          \
                        << " passed to CoupledFEMFunctionsx object! Quitting!");
    }

  return weight;
}

Number CoupledFEMFunctionsy::operator()(const FEMContext & c,
                                        const Point & p,
                                        const Real /* time */)
{
  Number weight = 0.0;

  switch(var)
    {
    case 1:
      {
        Gradient grad_C = c.point_gradient(3, p);

        weight = grad_C(1);
      }
      break;

    case 3:
      {
        Number v = c.point_value(1, p);

        weight = v;
      }
      break;

    default:
      libmesh_error_msg("Wrong variable number "                        \
                        << var                                          \
                        << " passed to CoupledFEMFunctionsy object! Quitting!");
    }

  return weight;
}
