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



#ifndef LIBMESH_DIFF_SOLVER_H
#define LIBMESH_DIFF_SOLVER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <vector>
#include <memory>

namespace libMesh
{

// Forward Declarations
class ImplicitSystem;
template <typename T> class NumericVector;

/**
 * Functor for use as callback in solve of nonlinear solver.
 */
class LinearSolutionMonitor
{
public:
  virtual void operator() (const NumericVector<Number> & delta_u,
                           const Real & norm_delta_u,
                           const NumericVector<Number> & u,
                           const Real & norm_u,
                           const NumericVector<Number> & res,
                           const Real & norm_res,
                           const unsigned int iteration) = 0;
  virtual ~LinearSolutionMonitor() = default;
};

/**
 * This is a generic class that defines a solver to handle
 * ImplicitSystem classes, including NonlinearImplicitSystem and
 * DifferentiableSystem   A user can define a solver by
 * deriving from this class and implementing certain functions.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2006-2010
 */
class DiffSolver : public ReferenceCountedObject<DiffSolver>,
                   public ParallelObject
{
public:
  /**
   * The type of system
   */
  typedef ImplicitSystem sys_type;

  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  DiffSolver (sys_type & s);

  /**
   * Factory method.  Requires a reference to the system
   * to be solved.
   *
   * \returns A NewtonSolver by default.
   */
  static std::unique_ptr<DiffSolver> build(sys_type & s);

  /**
   * Destructor.
   */
  virtual ~DiffSolver () = default;

  /**
   * The initialization function.  This method is used to
   * initialize internal data structures before a simulation begins.
   */
  virtual void init ();

  /**
   * The reinitialization function.  This method is used after
   * changes in the mesh.
   */
  virtual void reinit ();

  /**
   * This method performs a solve.  What occurs in
   * this method will depend on the type of solver.  See
   * the subclasses for more details.
   */
  virtual unsigned int solve () = 0;

  /**
   * \returns The number of "outer" (e.g. quasi-Newton) iterations
   * required by the last solve.
   */
  unsigned int total_outer_iterations() { return _outer_iterations; }

  /**
   * \returns The number of "inner" (e.g. Krylov) iterations
   * required by the last solve.
   */
  unsigned int total_inner_iterations() { return _inner_iterations; }

  /**
   * \returns The value of the SolveResult from the last solve.
   */
  unsigned int solve_result() { return _solve_result; }

  /**
   * \returns A constant reference to the system we are solving.
   */
  const sys_type & system () const { return _system; }

  /**
   * \returns A writable reference to the system we are solving.
   */
  sys_type & system () { return _system; }

  /**
   * Each linear solver step should exit after \p max_linear_iterations
   * is exceeded.
   */
  unsigned int max_linear_iterations;

  /**
   * The DiffSolver should exit in failure if \p max_nonlinear_iterations
   * is exceeded and \p continue_after_max_iterations is false, or should
   * end the nonlinear solve if \p max_nonlinear_iterations is exceeded and \p
   * continue_after_max_iterations is true.
   */
  unsigned int max_nonlinear_iterations;

  /**
   * The DiffSolver should not print anything to libMesh::out
   * unless quiet is set to false; default is true.
   */
  bool quiet;

  /**
   * The DiffSolver may print a lot more to libMesh::out
   * if verbose is set to true; default is false.
   */
  bool verbose;

  /**
   * Defaults to true, telling the DiffSolver to continue rather than exit when
   * a solve has reached its maximum number of nonlinear iterations.
   */
  bool continue_after_max_iterations;

  /**
   * Defaults to false, telling the DiffSolver to throw an error when
   * the backtracking scheme fails to find a descent direction.
   */
  bool continue_after_backtrack_failure;

  /**
   * The DiffSolver should exit after the residual is
   * reduced to either less than absolute_residual_tolerance
   * or less than relative_residual_tolerance times the
   * initial residual.
   *
   * Users should increase any of these tolerances that they want to use for a
   * stopping condition.
   */
  Real absolute_residual_tolerance;
  Real relative_residual_tolerance;

  /**
   * The DiffSolver should exit after the full nonlinear step norm is
   * reduced to either less than absolute_step_tolerance
   * or less than relative_step_tolerance times the largest
   * nonlinear solution which has been seen so far.
   *
   * Users should increase any of these tolerances that they want to use for a
   * stopping condition.
   */
  Real absolute_step_tolerance;
  Real relative_step_tolerance;

  /**
   * Any required linear solves will at first be done with this tolerance;
   * the DiffSolver may tighten the tolerance for later solves.
   */
  double initial_linear_tolerance;

  /**
   * The tolerance for linear solves is kept above this minimum
   */
  double minimum_linear_tolerance;

  /**
   * Enumeration return type for the solve() function.  Multiple SolveResults
   * may be combined (OR'd) in the single return.  To test which ones are present,
   * just AND the return value with any of the SolveResult flags defined below.
   */
  enum SolveResult {
    /**
     * A default or invalid solve result.  This usually means
     * no solve has occurred yet.
     */
    INVALID_SOLVE_RESULT = 0,

    /**
     * The solver converged but no
     * particular reason is specified.
     */
    CONVERGED_NO_REASON = 1,

    /**
     * The DiffSolver achieved the desired
     * absolute residual tolerance.
     */
    CONVERGED_ABSOLUTE_RESIDUAL = 2,

    /**
     * The DiffSolver achieved the desired
     * relative residual tolerance.
     */
    CONVERGED_RELATIVE_RESIDUAL = 4,

    /**
     * The DiffSolver achieved the desired
     * absolute step size tolerance.
     */
    CONVERGED_ABSOLUTE_STEP = 8,

    /**
     * The DiffSolver achieved the desired
     * relative step size tolerance.
     */
    CONVERGED_RELATIVE_STEP = 16,

    /**
     * The DiffSolver diverged but no
     * particular reason is specified.
     */
    DIVERGED_NO_REASON = 32,

    /**
     * The DiffSolver reached the maximum allowed
     * number of nonlinear iterations before satisfying
     * any convergence tests.
     */
    DIVERGED_MAX_NONLINEAR_ITERATIONS = 64,

    /**
     * The DiffSolver failed to find a descent direction
     * by backtracking (See newton_solver.C)
     */
    DIVERGED_BACKTRACKING_FAILURE = 128,

    /**
     * The linear solver used by the DiffSolver failed to
     * find a solution.
     */
    DIVERGED_LINEAR_SOLVER_FAILURE = 256
  };

  /**
   * Pointer to functor which is called right after each linear solve
   */
  std::unique_ptr<LinearSolutionMonitor> linear_solution_monitor;

  /**
   * Enable (or disable; it is \p true by default) exact enforcement
   * of constraints at the solver level, correcting any constrained
   * DoF coefficients in \p current_local_solution as well as applying
   * nonlinear residual and Jacobian terms based on constraint
   * equations.
   *
   * This is probably only safe to disable if user code is setting
   * nonlinear residual and Jacobian terms based on constraint
   * equations at an element-by-element level, by combining the
   * \p asymmetric_constraint_rows option with the
   * \p residual_constrain_element_vector processing option in
   * \p DofMap.
   */
  virtual void set_exact_constraint_enforcement(bool enable)
  {
    _exact_constraint_enforcement = enable;
  }

  bool exact_constraint_enforcement()
  {
    return _exact_constraint_enforcement;
  }

protected:

  /**
   * Whether we should enforce exact constraints globally during a
   * solve.
   */
  bool _exact_constraint_enforcement;

  /**
   * The largest solution norm which the DiffSolver has yet seen will be stored
   * here, to be used for stopping criteria based on relative_step_tolerance
   */
  Real max_solution_norm;

  /**
   * The largest nonlinear residual which the DiffSolver has yet seen will be
   * stored here, to be used for stopping criteria based on
   * relative_residual_tolerance
   */
  Real max_residual_norm;

  /**
   * The number of outer iterations used by the last solve
   */
  unsigned int _outer_iterations;

  /**
   * The number of inner iterations used by the last solve
   */
  unsigned int _inner_iterations;

  /**
   * A reference to the system we are solving.
   */
  sys_type & _system;

  /**
   * Initialized to zero.  solve_result is typically set internally in
   * the solve() function before it returns.  When non-zero,
   * solve_result tells the result of the latest solve.  See enum
   * definition for description.
   */
  unsigned int _solve_result;
};


} // namespace libMesh


#endif // LIBMESH_DIFF_SOLVER_H
