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


#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/petsc_diff_solver.h"
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_auto_fieldsplit.h"
#include "libmesh/boundary_info.h"

#ifdef LIBMESH_HAVE_PETSC

namespace libMesh
{

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed.
//
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  // Function to hand to PETSc's SNES,
  // which monitors convergence at X
  PetscErrorCode
  __libmesh_petsc_diff_solver_monitor (SNES snes,
                                       PetscInt its,
                                       PetscReal fnorm,
                                       void * ctx)
  {
    PetscFunctionBegin;

    PetscDiffSolver & solver =
      *(static_cast<PetscDiffSolver *> (ctx));

    if (solver.verbose)
      libMesh::out << "  PetscDiffSolver step " << its
                   << ", |residual|_2 = " << fnorm << std::endl;
    if (solver.linear_solution_monitor.get())
    {
      Vec petsc_delta_u;
      LibmeshPetscCall2(solver.comm(), SNESGetSolutionUpdate(snes, &petsc_delta_u));
      PetscVector<Number> delta_u(petsc_delta_u, solver.comm());
      delta_u.close();

      Vec petsc_u;
      LibmeshPetscCall2(solver.comm(), SNESGetSolution(snes, &petsc_u));
      PetscVector<Number> u(petsc_u, solver.comm());
      u.close();

      Vec petsc_res;
      LibmeshPetscCall2(solver.comm(), SNESGetFunction(snes, &petsc_res, nullptr, nullptr));
      PetscVector<Number> res(petsc_res, solver.comm());
      res.close();

      (*solver.linear_solution_monitor)(
                                        delta_u, delta_u.l2_norm(),
                                        u, u.l2_norm(),
                                        res, res.l2_norm(), its);
    }
    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

  // Functions to hand to PETSc's SNES,
  // which compute the residual or jacobian at X
  PetscErrorCode
  __libmesh_petsc_diff_solver_residual (SNES, Vec x, Vec r, void * ctx)
  {
    PetscFunctionBegin;

    libmesh_assert(x);
    libmesh_assert(r);
    libmesh_assert(ctx);

    PetscDiffSolver & solver =
      *(static_cast<PetscDiffSolver*> (ctx));
    ImplicitSystem & sys = solver.system();

    if (solver.verbose)
      libMesh::out << "Assembling the residual" << std::endl;

    PetscVector<Number> & X_system =
      *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> & R_system =
      *cast_ptr<PetscVector<Number> *>(sys.rhs);
    PetscVector<Number> X_input(x, sys.comm()), R_input(r, sys.comm());

    // DiffSystem assembles from the solution and into the rhs, so swap
    // those with our input vectors before assembling.  They'll probably
    // already be references to the same vectors, but PETSc might do
    // something tricky.
    X_input.swap(X_system);
    R_input.swap(R_system);

    // We may need to localize a parallel solution
    sys.update();

    // We may need to correct a non-conforming solution
    if (solver.exact_constraint_enforcement())
      sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

    // Do DiffSystem assembly
    sys.assembly(true, false, !solver.exact_constraint_enforcement());
    R_system.close();

    // Swap back
    X_input.swap(X_system);
    R_input.swap(R_system);

    // No errors, we hope
    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }


  PetscErrorCode
  __libmesh_petsc_diff_solver_jacobian (SNES,
                                        Vec x,
                                        Mat libmesh_dbg_var(j),
                                        Mat libmesh_dbg_var(pc),
                                        void * ctx)
  {
    PetscFunctionBegin;

    libmesh_assert(x);
    libmesh_assert(j);
    //  libmesh_assert_equal_to (pc, j);  // We don't use separate preconditioners yet
    libmesh_assert(ctx);

    PetscDiffSolver & solver =
      *(static_cast<PetscDiffSolver*> (ctx));
    ImplicitSystem & sys = solver.system();

    if (solver.verbose)
      libMesh::out << "Assembling the Jacobian" << std::endl;

    PetscVector<Number> & X_system =
      *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X_input(x, sys.comm());

    PetscMatrixBase<Number> & J_system =
      *cast_ptr<PetscMatrixBase<Number> *>(sys.matrix);
    libmesh_assert(J_system.mat() == pc);

    // DiffSystem assembles from the solution and into the jacobian, so
    // swap those with our input vectors before assembling.  They'll
    // probably already be references to the same vectors, but PETSc
    // might do something tricky.
    X_input.swap(X_system);

    // We may need to localize a parallel solution
    sys.update();

    // We may need to correct a non-conforming solution
    if (solver.exact_constraint_enforcement())
      sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

    // Do DiffSystem assembly
    sys.assembly(false, true, !solver.exact_constraint_enforcement());
    J_system.close();

    // Swap back
    X_input.swap(X_system);

    // No errors, we hope
    PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  }

} // extern "C"


PetscDiffSolver::PetscDiffSolver (sys_type & s)
  : Parent(s)
{
}


void PetscDiffSolver::init ()
{
  LOG_SCOPE("init()", "PetscDiffSolver");

  Parent::init();

  this->setup_petsc_data();
}



PetscDiffSolver::~PetscDiffSolver () = default;



void PetscDiffSolver::clear()
{
  LOG_SCOPE("clear()", "PetscDiffSolver");

  // calls custom deleter
  _snes.destroy();

#if !PETSC_VERSION_LESS_THAN(3,7,3)
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)
  _dm_wrapper.clear();
#endif
#endif
}



void PetscDiffSolver::reinit()
{
  LOG_SCOPE("reinit()", "PetscDiffSolver");

  // We need to wipe out all the old PETSc data
  // if we are reinit'ing, since we'll need to build
  // it all back up again.
  this->clear();

  Parent::reinit();

  this->setup_petsc_data();
}



DiffSolver::SolveResult convert_solve_result(SNESConvergedReason r)
{
  switch (r)
    {
    case SNES_CONVERGED_FNORM_ABS:
      return DiffSolver::CONVERGED_ABSOLUTE_RESIDUAL;
    case SNES_CONVERGED_FNORM_RELATIVE:
      return DiffSolver::CONVERGED_RELATIVE_RESIDUAL;
    case SNES_CONVERGED_SNORM_RELATIVE:
      return DiffSolver::CONVERGED_RELATIVE_STEP;
    case SNES_CONVERGED_ITS:
      // SNES_CONVERGED_TR_DELTA was changed to a diverged condition,
      // SNES_DIVERGED_TR_DELTA, in PETSc 1c6b2ff8df. This change will
      // likely be in 3.12 and later releases.
#if PETSC_VERSION_LESS_THAN(3,12,0)
    case SNES_CONVERGED_TR_DELTA:
#endif
      return DiffSolver::CONVERGED_NO_REASON;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
    case SNES_DIVERGED_FUNCTION_COUNT:
    case SNES_DIVERGED_FNORM_NAN:
    case SNES_DIVERGED_INNER:
    case SNES_DIVERGED_LINEAR_SOLVE:
    case SNES_DIVERGED_LOCAL_MIN:
      return DiffSolver::DIVERGED_NO_REASON;
    case SNES_DIVERGED_MAX_IT:
      return DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS;
    case SNES_DIVERGED_LINE_SEARCH:
      return DiffSolver::DIVERGED_BACKTRACKING_FAILURE;
      // In PETSc, SNES_CONVERGED_ITERATING means
      // the solve is still iterating, but by the
      // time we get here, we must have either
      // converged or diverged, so
      // SNES_CONVERGED_ITERATING is invalid.
    case SNES_CONVERGED_ITERATING:
      return DiffSolver::INVALID_SOLVE_RESULT;
    default:
      break;
    }
  return DiffSolver::INVALID_SOLVE_RESULT;
}



unsigned int PetscDiffSolver::solve()
{
  LOG_SCOPE("solve()", "PetscDiffSolver");

#if !PETSC_VERSION_LESS_THAN(3,7,3)
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)
  // GMG is currently not supported if we enable children to be associated with
  // boundary sides
  libmesh_assert(!_system.get_mesh().get_boundary_info().is_children_on_boundary_side());
#endif
#endif

  PetscVector<Number> & x =
    *(cast_ptr<PetscVector<Number> *>(_system.solution.get()));
  PetscMatrixBase<Number> & jac =
    *(cast_ptr<PetscMatrixBase<Number> *>(_system.matrix));
  PetscVector<Number> & r =
    *(cast_ptr<PetscVector<Number> *>(_system.rhs));

  LibmeshPetscCall(SNESSetFunction (_snes, r.vec(),
                                    __libmesh_petsc_diff_solver_residual, this));

  LibmeshPetscCall(SNESSetJacobian (_snes, jac.mat(), jac.mat(),
                                    __libmesh_petsc_diff_solver_jacobian, this));

  LibmeshPetscCall(SNESSetFromOptions(_snes));

  LibmeshPetscCall(SNESSolve (_snes, LIBMESH_PETSC_NULLPTR, x.vec()));

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  if (this->_exact_constraint_enforcement)
    _system.get_dof_map().enforce_constraints_exactly(_system);
#endif

  SNESConvergedReason reason;
  LibmeshPetscCall(SNESGetConvergedReason(_snes, &reason));

  PetscInt l_its, nl_its;
  LibmeshPetscCall(SNESGetLinearSolveIterations(_snes, &l_its));
  this->_inner_iterations = l_its;

  LibmeshPetscCall(SNESGetIterationNumber(_snes, &nl_its));
  this->_outer_iterations = nl_its;

  return convert_solve_result(reason);
}

void PetscDiffSolver::setup_petsc_data()
{
  LibmeshPetscCall(SNESCreate(this->comm().get(), _snes.get()));

  LibmeshPetscCall(SNESMonitorSet (_snes, __libmesh_petsc_diff_solver_monitor,
                                   this, LIBMESH_PETSC_NULLPTR));

  if (libMesh::on_command_line("--solver-system-names"))
    LibmeshPetscCall(SNESSetOptionsPrefix(_snes, (_system.name()+"_").c_str()));

  bool use_petsc_dm = libMesh::on_command_line("--use_petsc_dm");

  // This needs to be called before SNESSetFromOptions
#if !PETSC_VERSION_LESS_THAN(3,7,3)
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)
  if (use_petsc_dm)
    this->_dm_wrapper.init_and_attach_petscdm(_system, _snes);
#endif
#endif

  // If we're not using PETSc DM, let's keep around
  // the old style for fieldsplit
  if (!use_petsc_dm)
    {
      LibmeshPetscCall(SNESSetFromOptions(_snes));

      KSP my_ksp;
      LibmeshPetscCall(SNESGetKSP(_snes, &my_ksp));

      PC my_pc;
      LibmeshPetscCall(KSPGetPC(my_ksp, &my_pc));

      petsc_auto_fieldsplit(my_pc, _system);
    }
}

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
