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


// libMesh includes
#include "libmesh/diff_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/time_solver.h"
#include "libmesh/unsteady_solver.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/zero_function.h"

// C++ includes
#include <utility> // std::swap

namespace libMesh
{



DifferentiableSystem::DifferentiableSystem(EquationSystems & es,
                                           const std::string & name_in,
                                           const unsigned int number_in) :
  Parent      (es, name_in, number_in),
  time_solver (),
  deltat(1.),
  postprocess_sides(false),
  print_solution_norms(false),
  print_solutions(false),
  print_residual_norms(false),
  print_residuals(false),
  print_jacobian_norms(false),
  print_jacobians(false),
  print_element_solutions(false),
  print_element_residuals(false),
  print_element_jacobians(false),
  _constrain_in_solver(true),
  _diff_physics(),
  _diff_qoi()
{
}



DifferentiableSystem::~DifferentiableSystem () = default;



void DifferentiableSystem::clear ()
{
  // If we had no attached Physics object, clear our own Physics data
  if (this->_diff_physics.empty())
    this->clear_physics();

  this->_diff_physics = {}; // No stack::clear
  this->_diff_qoi = {};

  // If we had no attached QoI object, clear our own QoI data
  if (this->_diff_qoi.empty())
    this->clear_qoi();

  use_fixed_solution = false;
}



void DifferentiableSystem::reinit ()
{
  Parent::reinit();

  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);

  time_solver->reinit();
}



void DifferentiableSystem::init_data ()
{
  // If it isn't a separate initialized-upon-attachment object, do any
  // initialization our physics needs.
  if (this->_diff_physics.empty())
    this->init_physics(*this);

  // Do any initialization our solvers need
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);

  // Now check for second order variables and add their velocities to the System.
  if (!time_solver->is_steady())
    {
      const UnsteadySolver & unsteady_solver =
        cast_ref<const UnsteadySolver &>(*(time_solver.get()));

      if (unsteady_solver.time_order() == 1)
        this->add_second_order_dot_vars();
    }

  time_solver->init();

  // Next initialize ImplicitSystem data
  Parent::init_data();

  time_solver->init_data();
}

std::unique_ptr<DiffContext> DifferentiableSystem::build_context ()
{
  auto context = std::make_unique<DiffContext>(*this);
  context->set_deltat_pointer( &this->deltat );
  return context;
}


void DifferentiableSystem::assemble ()
{
  this->assembly(true, true);
}



void DifferentiableSystem::solve ()
{
  // Get the time solver object associated with the system, and tell it that
  // we are not solving the adjoint problem
  this->get_time_solver().set_is_adjoint(false);

  libmesh_assert_equal_to (&(time_solver->system()), this);
  time_solver->solve();
}



std::pair<unsigned int, Real> DifferentiableSystem::adjoint_solve (const QoISet & qoi_indices)
{
  // Get the time solver object associated with the system, and tell it that
  // we are solving the adjoint problem
  this->get_time_solver().set_is_adjoint(true);

  return time_solver->adjoint_solve(qoi_indices);

  //return this->ImplicitSystem::adjoint_solve(qoi_indices);
}



LinearSolver<Number> * DifferentiableSystem::get_linear_solver() const
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return this->time_solver->linear_solver().get();
}



std::pair<unsigned int, Real> DifferentiableSystem::get_linear_solve_parameters() const
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return std::make_pair(this->time_solver->diff_solver()->max_linear_iterations,
                        this->time_solver->diff_solver()->relative_residual_tolerance);
}



void DifferentiableSystem::add_second_order_dot_vars()
{
  const std::set<unsigned int> & second_order_vars = this->get_second_order_vars();
  if (!second_order_vars.empty())
    {
      for (const auto & var_id : second_order_vars)
        {
          const Variable & var = this->variable(var_id);
          std::string new_var_name = std::string("dot_")+var.name();

          unsigned int v_var_idx;

          if (var.active_subdomains().empty())
            v_var_idx = this->add_variable( new_var_name, var.type() );
          else
            v_var_idx = this->add_variable( new_var_name, var.type(), &var.active_subdomains() );

          _second_order_dot_vars.insert(std::pair<unsigned int, unsigned int>(var_id, v_var_idx));

          // The new velocities are time evolving variables of first order
          this->time_evolving( v_var_idx, 1 );

#ifdef LIBMESH_ENABLE_DIRICHLET
          // And if there are any boundary conditions set on the second order
          // variable, we also need to set it on its velocity variable.
          this->add_dot_var_dirichlet_bcs(var_id, v_var_idx);
#endif
        }
    }
}

#ifdef LIBMESH_ENABLE_DIRICHLET
void DifferentiableSystem::add_dot_var_dirichlet_bcs( unsigned int var_idx,
                                                      unsigned int dot_var_idx )
{
  // We're assuming that there could be a lot more variables than
  // boundary conditions, so we search each of the boundary conditions
  // for this variable rather than looping over boundary conditions
  // in a separate loop and searching through all the variables.
  const DirichletBoundaries * all_dbcs =
    this->get_dof_map().get_dirichlet_boundaries();

  if (all_dbcs)
    {
      // We need to cache the DBCs to be added so that we add them
      // after looping over the existing DBCs. Otherwise, we're polluting
      // the thing we're looping over.
      std::vector<DirichletBoundary> new_dbcs;

      for (const auto & dbc : *all_dbcs)
        {
          libmesh_assert(dbc);

          // Look for second order variable in the current
          // DirichletBoundary object
          std::vector<unsigned int>::const_iterator dbc_var_it =
            std::find( dbc->variables.begin(), dbc->variables.end(), var_idx );

          // If we found it, then we also need to add it's corresponding
          // "dot" variable to a DirichletBoundary
          std::vector<unsigned int> vars_to_add;
          if (dbc_var_it != dbc->variables.end())
            vars_to_add.push_back(dot_var_idx);

          if (!vars_to_add.empty())
            {
              // We need to check if the boundary condition is time-dependent.
              // Currently, we cannot automatically differentiate w.r.t. time
              // so if the user supplies a time-dependent Dirichlet BC, then
              // we can't automatically support the Dirichlet BC for the
              // "velocity" boundary condition, so we error. Otherwise,
              // the "velocity boundary condition will just be zero.
              bool is_time_evolving_bc = false;
              if (dbc->f)
                is_time_evolving_bc = dbc->f->is_time_dependent();
              else if (dbc->f_fem)
                // We it's a FEMFunctionBase object, it will be implicitly
                // time-dependent since it is assumed to depend on the solution.
                is_time_evolving_bc = true;
              else
                libmesh_error_msg("Could not find valid boundary function!");

              libmesh_error_msg_if(is_time_evolving_bc, "Cannot currently support time-dependent Dirichlet BC for dot variables!");
              libmesh_error_msg_if(!dbc->f, "Expected valid DirichletBoundary function");

              new_dbcs.emplace_back(dbc->b, vars_to_add, ZeroFunction<Number>());
            }
        }

      // Let the DofMap make its own deep copy of the DirichletBC objects
      for (const auto & dbc : new_dbcs)
        this->get_dof_map().add_dirichlet_boundary(dbc);

    } // if (all_dbcs)
}
#endif // LIBMESH_ENABLE_DIRICHLET

void DifferentiableSystem::attach_qoi( DifferentiableQoI * qoi_in )
{
  this->_diff_qoi = {};
  this->_diff_qoi.push(qoi_in->clone());

  auto & dq = this->_diff_qoi.top();
  // User needs to resize qoi system qoi accordingly
#ifdef LIBMESH_ENABLE_DEPRECATED
  // Call the old API for backwards compatibility
  dq->init_qoi( this->qoi );

  // Then the new API for forwards compatibility
  dq->init_qoi_count( *this );
#else
#ifndef NDEBUG
  // Make sure the user has updated their QoI subclass - call the old
  // API and make sure it does nothing
  std::vector<Number> deprecated_vector;
  dq->init_qoi( deprecated_vector );
  libmesh_assert(deprecated_vector.empty());
#endif

  // Then the new API
  dq->init_qoi_count( *this );
#endif
}

unsigned int DifferentiableSystem::get_second_order_dot_var( unsigned int var ) const
{
  // For SteadySolver or SecondOrderUnsteadySolvers, we just give back var
  unsigned int dot_var = var;

  if (!time_solver->is_steady())
    {
      const UnsteadySolver & unsteady_solver =
        cast_ref<const UnsteadySolver &>(*(time_solver.get()));

      if (unsteady_solver.time_order() == 1)
        dot_var = this->_second_order_dot_vars.find(var)->second;
    }

  return dot_var;
}

bool DifferentiableSystem::have_first_order_scalar_vars() const
{
  bool have_first_order_scalar_vars = false;

  if (this->have_first_order_vars())
    for (const auto & var : this->get_first_order_vars())
      if (this->variable(var).type().family == SCALAR)
        have_first_order_scalar_vars = true;

  return have_first_order_scalar_vars;
}

bool DifferentiableSystem::have_second_order_scalar_vars() const
{
  bool have_second_order_scalar_vars = false;

  if (this->have_second_order_vars())
    for (const auto & var : this->get_second_order_vars())
      if (this->variable(var).type().family == SCALAR)
        have_second_order_scalar_vars = true;

  return have_second_order_scalar_vars;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
void DifferentiableSystem::swap_physics ( DifferentiablePhysics * & swap_physics )
{
  // This isn't safe if users aren't very careful about memory
  // management and they don't (or aren't able to due to an exception)
  // swap back.
  libmesh_deprecated();

  // A mess of code for backwards compatibility
  if (this->_diff_physics.empty())
    {
      // Swap-something-else-for-self
      std::unique_ptr<DifferentiablePhysics> scary_hack(swap_physics);
      this->_diff_physics.push(std::move(scary_hack));
      swap_physics = this;
    }
  else if (swap_physics == this)
    {
      // The user must be cleaning up after a previous
      // swap-something-else-for-self
      libmesh_assert(!this->_diff_physics.empty());

      // So we don't want to delete what got swapped in, but we do
      // want to put it back into their pointer
      DifferentiablePhysics * old_p = this->_diff_physics.top().release();
      this->_diff_physics.pop();
      swap_physics = old_p;

      // And if the user is doing anything more sophisticated than
      // that then the user is sophisticated enough to upgrade to
      // push/pop.
      libmesh_assert(this->_diff_physics.empty());
    }
  else
    {
      // Swapping one external physics for another
      DifferentiablePhysics * old_p = this->_diff_physics.top().release();
      std::swap(old_p, swap_physics);
      this->_diff_physics.top().reset(old_p);
    }

  // If the physics has been swapped, we will reassemble
  // the matrix from scratch before doing an adjoint solve
  // rather than just transposing
  this->disable_cache();
}
#endif // LIBMESH_ENABLE_DEPRECATED



void DifferentiableSystem::push_physics ( DifferentiablePhysics & new_physics )
{
  this->_diff_physics.push(new_physics.clone_physics());

  // If the physics has been changed, we will reassemble
  // the matrix from scratch before doing an adjoint solve
  // rather than just transposing
  this->disable_cache();
}



void DifferentiableSystem::pop_physics ()
{
  libmesh_assert(!this->_diff_physics.empty());

  this->_diff_physics.pop();

  // If the physics has been changed, we will reassemble
  // the matrix from scratch before doing an adjoint solve
  // rather than just transposing
  this->disable_cache();
}


void DifferentiableSystem::set_constrain_in_solver(bool enable)
{
  _constrain_in_solver = enable;
  this->time_solver->diff_solver()->set_exact_constraint_enforcement(enable);
}


} // namespace libMesh
