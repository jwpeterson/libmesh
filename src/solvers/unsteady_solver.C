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


#include "libmesh/unsteady_solver.h"

#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/diff_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/error_vector.h"
#include "libmesh/int_range.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/sensitivity_data.h"
#include "libmesh/solution_history.h"

namespace libMesh
{



UnsteadySolver::UnsteadySolver (sys_type & s)
  : TimeSolver(s),
    old_local_nonlinear_solution (NumericVector<Number>::build(s.comm()).release()),
    first_solve                  (true),
    first_adjoint_step (true)
{
  old_adjoints.resize(s.n_qois());

  // Set the old adjoint pointers to nullptrs
  // We will use this nullness to skip the initial time instant,
  // when there is no older adjoint.
  for(auto j : make_range(s.n_qois()))
  {
    old_adjoints[j] = nullptr;
  }
}



UnsteadySolver::~UnsteadySolver () = default;



void UnsteadySolver::init ()
{
  TimeSolver::init();

  _system.add_vector("_old_nonlinear_solution");
}


void UnsteadySolver::init_adjoints ()
{
  TimeSolver::init_adjoints();

  // Add old adjoint solutions
  // To keep the number of vectors consistent between the primal and adjoint
  // time loops, we will also add the adjoint rhs vector during initialization
  for(auto i : make_range(_system.n_qois()))
  {
    std::string old_adjoint_solution_name = "_old_adjoint_solution";
    old_adjoint_solution_name+= std::to_string(i);
    _system.add_vector(old_adjoint_solution_name, false, GHOSTED);

    std::string adjoint_rhs_name = "adjoint_rhs";
    adjoint_rhs_name+= std::to_string(i);
    _system.add_vector(adjoint_rhs_name, false, GHOSTED);
  }

}


void UnsteadySolver::init_data()
{
  TimeSolver::init_data();

#ifdef LIBMESH_ENABLE_GHOSTED
  old_local_nonlinear_solution->init (_system.n_dofs(), _system.n_local_dofs(),
                                      _system.get_dof_map().get_send_list(), false,
                                      GHOSTED);
#else
  old_local_nonlinear_solution->init (_system.n_dofs(), false, SERIAL);
#endif
}



void UnsteadySolver::reinit ()
{
  TimeSolver::reinit();

#ifdef LIBMESH_ENABLE_GHOSTED
  old_local_nonlinear_solution->init (_system.n_dofs(), _system.n_local_dofs(),
                                      _system.get_dof_map().get_send_list(), false,
                                      GHOSTED);
#else
  old_local_nonlinear_solution->init (_system.n_dofs(), false, SERIAL);
#endif

  // localize the old solution
  NumericVector<Number> & old_nonlinear_soln =
    _system.get_vector("_old_nonlinear_solution");

  old_nonlinear_soln.localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}



void UnsteadySolver::solve ()
{
  if (first_solve)
    {
      advance_timestep();
      first_solve = false;
    }

  unsigned int solve_result = _diff_solver->solve();

  // If we requested the UnsteadySolver to attempt reducing dt after a
  // failed DiffSolver solve, check the results of the solve now.
  if (reduce_deltat_on_diffsolver_failure)
    {
      bool backtracking_failed =
        solve_result & DiffSolver::DIVERGED_BACKTRACKING_FAILURE;

      bool max_iterations =
        solve_result & DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;

      if (backtracking_failed || max_iterations)
        {
          // Cut timestep in half
          for (unsigned int nr=0; nr<reduce_deltat_on_diffsolver_failure; ++nr)
            {
              _system.deltat *= 0.5;
              libMesh::out << "Newton backtracking failed.  Trying with smaller timestep, dt="
                           << _system.deltat << std::endl;

              solve_result = _diff_solver->solve();

              // Check solve results with reduced timestep
              bool backtracking_still_failed =
                solve_result & DiffSolver::DIVERGED_BACKTRACKING_FAILURE;

              bool backtracking_max_iterations =
                solve_result & DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;

              if (!backtracking_still_failed && !backtracking_max_iterations)
                {
                  // Set the successful deltat as the last deltat
                  last_deltat = _system.deltat;

                  if (!quiet)
                    libMesh::out << "Reduced dt solve succeeded." << std::endl;
                  return;
                }
            }

          // If we made it here, we still couldn't converge the solve after
          // reducing deltat
          libMesh::out << "DiffSolver::solve() did not succeed after "
                       << reduce_deltat_on_diffsolver_failure
                       << " attempts." << std::endl;
          libmesh_convergence_failure();

        } // end if (backtracking_failed || max_iterations)
    } // end if (reduce_deltat_on_diffsolver_failure)

  // Set the successful deltat as the last deltat
  last_deltat = _system.deltat;
}



void UnsteadySolver::advance_timestep ()
{
  // The first access of advance_timestep happens via solve, not user code
  // It is used here to store any initial conditions data
  if (!first_solve)
    {
      // We call advance_timestep in user code after solve, so any solutions
      // we will be storing will be for the next time instance
      _system.time += _system.deltat;
    }
    else
    {
      // We are here because of a call to advance_timestep that happens
      // via solve, the very first solve. All we are doing here is storing
      // the initial condition. The actual solution computed via this solve
      // will be stored when we call advance_timestep in the user's timestep loop
      first_solve = false;
    }

  // If the user has attached a memory or file solution history object
  // to the solver, this will store the current solution indexed with
  // the current time
  solution_history->store(false, _system.time);

  NumericVector<Number> & old_nonlinear_soln =
    _system.get_vector("_old_nonlinear_solution");
  NumericVector<Number> & nonlinear_solution =
    *(_system.solution);

  old_nonlinear_soln = nonlinear_solution;

  old_nonlinear_soln.localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}

std::pair<unsigned int, Real> UnsteadySolver::adjoint_solve(const QoISet & qoi_indices)
{
  std::pair<unsigned int, Real> adjoint_output = _system.ImplicitSystem::adjoint_solve(qoi_indices);

  // Record the deltat we used for this adjoint timestep. This was determined completely
  // by SolutionHistory::retrieve methods. The adjoint_solve methods should never change deltat.
  last_deltat = _system.deltat;

  return adjoint_output;
}

void UnsteadySolver::adjoint_advance_timestep ()
{
  // Call the store function to store the adjoint we have computed (or
  // for first_adjoint_step, the adjoint initial condition) in this
  // time step for the time instance.
  solution_history->store(true, _system.time);

  // Before moving to the next time instant, copy over the current adjoint solutions into _old_adjoint_solutions
  for(auto i : make_range(_system.n_qois()))
  {
   std::string old_adjoint_solution_name = "_old_adjoint_solution";
   old_adjoint_solution_name+= std::to_string(i);
   NumericVector<Number> & old_adjoint_solution_i = _system.get_vector(old_adjoint_solution_name);
   NumericVector<Number> & adjoint_solution_i = _system.get_adjoint_solution(i);
   old_adjoint_solution_i = adjoint_solution_i;
  }

  _system.time -= _system.deltat;

  // Retrieve the primal solution vectors at this new (or for
  // first_adjoint_step, initial) time instance. These provide the
  // data to solve the adjoint problem for the next time instance.
  solution_history->retrieve(true, _system.time);

  // Dont forget to localize the old_nonlinear_solution !
  _system.get_vector("_old_nonlinear_solution").localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}

void UnsteadySolver::retrieve_timestep()
{
  // Retrieve all the stored vectors at the current time
  solution_history->retrieve(false, _system.time);

  // Dont forget to localize the old_nonlinear_solution !
  _system.get_vector("_old_nonlinear_solution").localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}

void UnsteadySolver::integrate_adjoint_sensitivity(const QoISet & qois, const ParameterVector & parameter_vector, SensitivityData & sensitivities)
{
  // CURRENTLY using the trapezoidal rule to integrate each timestep
  // (f(t_j) + f(t_j+1))/2 (t_j+1 - t_j)
  // Fix me: This function needs to be moved to the EulerSolver classes like the
  // other integrate_timestep functions, and use an integration rule consistent with
  // the theta method used for the time integration.

  // Get t_j
  Real time_left = _system.time;

  // Left side sensitivities to hold f(t_j)
  SensitivityData sensitivities_left(qois, _system, parameter_vector);

  // Get f(t_j)
  _system.adjoint_qoi_parameter_sensitivity(qois, parameter_vector, sensitivities_left);

  // Advance to t_j+1
  _system.time = _system.time + _system.deltat;

  // Get t_j+1
  Real time_right = _system.time;

  // Right side sensitivities f(t_j+1)
  SensitivityData sensitivities_right(qois, _system, parameter_vector);

  // Remove the sensitivity rhs vector from system since we did not write it to file and it cannot be retrieved
  _system.remove_vector("sensitivity_rhs0");

  // Retrieve the primal and adjoint solutions at the current timestep
  retrieve_timestep();

  // Get f(t_j+1)
  _system.adjoint_qoi_parameter_sensitivity(qois, parameter_vector, sensitivities_right);

  // Remove the sensitivity rhs vector from system since we did not write it to file and it cannot be retrieved
  _system.remove_vector("sensitivity_rhs0");

  // Get the contributions for each sensitivity from this timestep
  const auto pv_size = parameter_vector.size();
  for (auto i : make_range(qois.size(_system)))
    for (auto j : make_range(pv_size))
      sensitivities[i][j] = ( (sensitivities_left[i][j] + sensitivities_right[i][j])/2. )*(time_right - time_left);
}

void UnsteadySolver::integrate_qoi_timestep()
{
  libmesh_not_implemented();
}

#ifdef LIBMESH_ENABLE_AMR
void UnsteadySolver::integrate_adjoint_refinement_error_estimate(AdjointRefinementEstimator & /*adjoint_refinement_error_estimator*/, ErrorVector & /*QoI_elementwise_error*/)
{
  libmesh_not_implemented();
}
#endif // LIBMESH_ENABLE_AMR

Number UnsteadySolver::old_nonlinear_solution(const dof_id_type global_dof_number)
  const
{
  libmesh_assert_less (global_dof_number, _system.get_dof_map().n_dofs());
  libmesh_assert_less (global_dof_number, old_local_nonlinear_solution->size());

  return (*old_local_nonlinear_solution)(global_dof_number);
}



Real UnsteadySolver::du(const SystemNorm & norm) const
{

  std::unique_ptr<NumericVector<Number>> solution_copy =
    _system.solution->clone();

  solution_copy->add(-1., _system.get_vector("_old_nonlinear_solution"));

  solution_copy->close();

  return _system.calculate_norm(*solution_copy, norm);
}

void UnsteadySolver::update()
{
  // Dont forget to localize the old_nonlinear_solution !
  _system.get_vector("_old_nonlinear_solution").localize
    (*old_local_nonlinear_solution,
     _system.get_dof_map().get_send_list());
}

} // namespace libMesh
