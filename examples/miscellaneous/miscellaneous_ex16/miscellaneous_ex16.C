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

// <h1>Miscellaneous Example 16 - Static Condensation</h1>
// \author Alexander D. Lindsay
// \date 2024
//
// This example is equivalent to the third introductory example except that
// we add and use a StaticCondensation object. The StaticCondensation class
// forward eliminates internal degrees of freedom and then performs a global
// solve only with the trace degrees of freedom. We then perform backward
// substitution to recover the solution for the internal degrees of freedom.
// We verify that this solution matches the solution when the global solve
// is performed on both interior and trace degrees of freedom

#include "libmesh/libmesh_config.h"

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/static_condensation.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// For mesh refinement
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

// I/O utilities.
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"

// For the solver for the system with static condensation
#include "libmesh/petsc_linear_solver.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

#if defined(LIBMESH_HAVE_EIGEN_DENSE) && defined(LIBMESH_HAVE_PETSC)
// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void assemble_poisson(EquationSystems & es, const std::string & system_name);
#endif

// Function prototype for the exact solution.
Real exact_solution(const Real x, const Real y, const Real z = 0.);

int
main(int argc, char ** argv)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init(argc, argv);

#if !defined(LIBMESH_HAVE_EIGEN_DENSE)
  libmesh_example_requires(false, "--enable-eigen");
#elif !defined(LIBMESH_HAVE_PETSC)
  libmesh_example_requires(false, "--enable-petsc");
#else

  // Brief message to the user regarding the program name
  // and command line arguments.
  libMesh::out << "Running " << argv[0];

  for (int i = 1; i < argc; i++)
    libMesh::out << " " << argv[i];
  libMesh::out << std::endl << std::endl;

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Parse the input file.
  GetPot infile("miscellaneous_ex16.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  const unsigned int grid_size = infile("grid_size", 5);

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
  // elements instead of the default QUAD4's we used in example 2
  // allow us to use higher-order approximation.
  MeshTools::Generation::build_square(mesh, grid_size, grid_size, -1., 1., -1., 1., QUAD9);

  // Print information about the mesh to the screen.
  // Note that 5x5 QUAD9 elements actually has 11x11 nodes,
  // so this mesh is significantly larger than the one in example 2.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // Declare the Poisson system and its variables.
  // The Poisson system is another example of a steady system.
  auto & sys = equation_systems.add_system<LinearImplicitSystem>("Poisson");

  // Adds the variable "u" to "Poisson".  "u"
  // will be approximated using second-order approximation.
  sys.add_variable("u", SECOND);

  // Give the system a pointer to the matrix assembly
  // function.  This will be called when needed by the
  // library.
  sys.attach_assemble_function(assemble_poisson);

  // Now perform same steps for system with static condensation enabled
  auto & sc_sys = equation_systems.add_system<LinearImplicitSystem>("SC_Poisson");
  sc_sys.add_variable("v", SECOND);
  sc_sys.attach_assemble_function(assemble_poisson);
  sc_sys.create_static_condensation();

#ifdef LIBMESH_ENABLE_AMR
  // Define the mesh refinement object that takes care of adaptively
  // refining the mesh.
  MeshRefinement mesh_refinement(mesh);

  // These parameters determine the proportion of elements that will
  // be refined and coarsened. Any element within 30% of the maximum
  // error on any element will be refined, and any element within 30%
  // of the minimum error on any element might be coarsened
  mesh_refinement.refine_fraction()  = 0.7;
  mesh_refinement.coarsen_fraction() = 0.3;
  // We won't refine any element more than 2 times in total
  mesh_refinement.max_h_level()      = 2;
#endif

  // Initialize the data structures for the equation system.
  equation_systems.init();

#ifdef LIBMESH_ENABLE_AMR
  // Refinement parameters
  const unsigned int max_r_steps = 2; // Refine the mesh 2 times

  for (const auto r_step : make_range(max_r_steps + 1))
    {
#endif
      // Prints information about the system to the screen.
      equation_systems.print_info();

      // Solve
      sys.solve();
      auto * sc_solver = dynamic_cast<PetscLinearSolver<Number> *>(sc_sys.get_linear_solver());
      libmesh_assert(sc_solver);
      KSP sc_ksp = sc_solver->ksp();
      LibmeshPetscCall2(sc_solver->comm(), KSPSetType(sc_ksp, KSPPREONLY));
      LibmeshPetscCall2(sc_solver->comm(), KSPSetInitialGuessNonzero(sc_ksp, PETSC_FALSE));
      sc_sys.solve();

      libmesh_error_msg_if(!libMesh::relative_fuzzy_equals(*sys.solution, *sc_sys.solution, 1e-4),
                           "mismatching solution");
      libMesh::out << "Static condensation reduced problem size to "
                   << sc_sys.get_static_condensation().get_condensed_mat().m() << std::endl << std::endl;

#if defined(LIBMESH_HAVE_EXODUS_API) && !defined(LIBMESH_ENABLE_PARMESH)
      // After solving the system write the solution
      // to an Exodus-formatted file.
      ExodusII_IO exii_io(mesh);
      const std::string file_name =
#ifdef LIBMESH_ENABLE_AMR
          "out_" + std::to_string(r_step) + ".e";
#else
          "out.e";
#endif
      exii_io.write_equation_systems(file_name, equation_systems);
#endif

#ifdef LIBMESH_ENABLE_AMR
      // We need to ensure that the mesh is not refined on the last iteration
      // of this loop, since we do not want to refine the mesh unless we are
      // going to solve the equation system for that refined mesh.
      if (r_step != max_r_steps)
        {
          // Error estimation objects, see Adaptivity Example 2 for details
          ErrorVector error;
          KellyErrorEstimator error_estimator;
          error_estimator.use_unweighted_quadrature_rules = true;

          // Compute the error for each active element
          error_estimator.estimate_error(sys, error);

          // Output error estimate magnitude
          libMesh::out << "Error estimate\nl2 norm = " << error.l2_norm()
                       << "\nmaximum = " << error.maximum() << std::endl << std::endl;

          // Flag elements to be refined and coarsened
          mesh_refinement.flag_elements_by_error_fraction(error);

          // Perform refinement and coarsening
          mesh_refinement.refine_and_coarsen_elements();

          // Reinitialize the equation_systems object for the newly refined
          // mesh. One of the steps in this is project the solution onto the
          // new mesh
          equation_systems.reinit();
        }
    }
#endif // LIBMESH_ENABLE_AMR
#endif // defined(LIBMESH_HAVE_EIGEN_DENSE) && defined(LIBMESH_HAVE_PETSC)

  // All done.
  return 0;
}

#if defined(LIBMESH_HAVE_EIGEN_DENSE) && defined(LIBMESH_HAVE_PETSC)

// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void
assemble_poisson(EquationSystems & es, const std::string & system_name)
{
  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>(system_name);

  // Get a pointer to the StaticCondensation class if it exists
  StaticCondensation * sc = nullptr;
  if (system.has_static_condensation())
    sc = &system.get_static_condensation();

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.  Introduction Example 4
  // describes some advantages of  std::unique_ptr's in the context of
  // quadrature rules.
  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule(dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule(&qrule);

  // Declare a special finite element object for
  // boundary integration.
  std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));

  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim - 1, FIFTH);

  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point> & q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // The global system matrix
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element ranges are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // range will iterate from the first to the last element on
  // the local processor.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the
  // active_local_element_ptr_range.
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    // Cache the number of degrees of freedom on this element, for
    // use as a loop bound later.  We use cast_int to explicitly
    // convert from size() (which may be 64-bit) to unsigned int
    // (which may be 32-bit but which is definitely enough to count
    // *local* degrees of freedom.
    const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit(elem);

    // With one variable, we should have the same number of degrees
    // of freedom as shape functions.
    libmesh_assert_equal_to(n_dofs, phi.size());

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).

    // The  DenseMatrix::resize() and the  DenseVector::resize()
    // members will automatically zero out the matrix  and vector.
    Ke.resize(n_dofs, n_dofs);

    Fe.resize(n_dofs);

    // Now loop over the quadrature points.  This handles
    // the numeric integration.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      // Now we will build the element matrix.  This involves
      // a double loop to integrate the test functions (i) against
      // the trial functions (j).
      for (unsigned int i = 0; i != n_dofs; i++)
        for (unsigned int j = 0; j != n_dofs; j++)
        {
          Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
        }

      // This is the end of the matrix summation loop
      // Now we build the element right-hand-side contribution.
      // This involves a single loop in which we integrate the
      // "forcing function" in the PDE against the test functions.
      {
        const Real x = q_point[qp](0);
        const Real y = q_point[qp](1);
        const Real eps = 1.e-3;

        // "fxy" is the forcing function for the Poisson equation.
        // In this case we set fxy to be a finite difference
        // Laplacian approximation to the (known) exact solution.
        //
        // We will use the second-order accurate FD Laplacian
        // approximation, which in 2D is
        //
        // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
        //                u(i-1,j) + u(i+1,j) +
        //                -4*u(i,j))/h^2
        //
        // Since the value of the forcing function depends only
        // on the location of the quadrature point (q_point[qp])
        // we will compute it here, outside of the i-loop
        const Real fxy =
            -(exact_solution(x, y - eps) + exact_solution(x, y + eps) + exact_solution(x - eps, y) +
              exact_solution(x + eps, y) - 4. * exact_solution(x, y)) /
            eps / eps;

        for (unsigned int i = 0; i != n_dofs; i++)
          Fe(i) += JxW[qp] * fxy * phi[i][qp];
      }
    }

    // We have now reached the end of the RHS summation,
    // and the end of quadrature point loop, so
    // the interior element integration has
    // been completed.  However, we have not yet addressed
    // boundary conditions.  For this example we will only
    // consider simple Dirichlet boundary conditions.
    //
    // There are several ways Dirichlet boundary conditions
    // can be imposed.  A simple approach, which works for
    // interpolary bases like the standard Lagrange polynomials,
    // is to assign function values to the
    // degrees of freedom living on the domain boundary. This
    // works well for interpolary bases, but is more difficult
    // when non-interpolary (e.g Legendre or Hierarchic) bases
    // are used.
    //
    // Dirichlet boundary conditions can also be imposed with a
    // "penalty" method.  In this case essentially the L2 projection
    // of the boundary values are added to the matrix. The
    // projection is multiplied by some large factor so that, in
    // floating point arithmetic, the existing (smaller) entries
    // in the matrix and right-hand-side are effectively ignored.
    //
    // This amounts to adding a term of the form (in latex notation)
    //
    // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta
    // \Omega} u \phi_i
    //
    // where
    //
    // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
    {

      // The following loop is over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (auto side : elem->side_index_range())
        if (elem->neighbor_ptr(side) == nullptr)
        {
          // The value of the shape functions at the quadrature
          // points.
          const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();

          // The Jacobian * Quadrature Weight at the quadrature
          // points on the face.
          const std::vector<Real> & JxW_face = fe_face->get_JxW();

          // The XYZ locations (in physical space) of the
          // quadrature points on the face.  This is where
          // we will interpolate the boundary value function.
          const std::vector<Point> & qface_point = fe_face->get_xyz();

          // Compute the shape function values on the element
          // face.
          fe_face->reinit(elem, side);

          // Some shape functions will be 0 on the face, but for
          // ease of indexing and generality of code we loop over
          // them anyway
          libmesh_assert_equal_to(n_dofs, phi_face.size());

          // Loop over the face quadrature points for integration.
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            // The location on the boundary of the current
            // face quadrature point.
            const Real xf = qface_point[qp](0);
            const Real yf = qface_point[qp](1);

            // The penalty value.  \frac{1}{\epsilon}
            // in the discussion above.
            const Real penalty = 1.e10;

            // The boundary value.
            const Real value = exact_solution(xf, yf);

            // Matrix contribution of the L2 projection.
            for (unsigned int i = 0; i != n_dofs; i++)
              for (unsigned int j = 0; j != n_dofs; j++)
                Ke(i, j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];

            // Right-hand-side contribution of the L2
            // projection.
            for (unsigned int i = 0; i != n_dofs; i++)
              Fe(i) += JxW_face[qp] * penalty * value * phi_face[i][qp];
          }
        }
    }

    // We have now finished the quadrature point loop,
    // and have therefore applied all the boundary conditions.

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

    if (sc)
      sc->set_current_elem(*elem);

    // The element matrix and right-hand-side are now built
    // for this element.  Add them to the global matrix and
    // right-hand-side vector.  The  SparseMatrix::add_matrix()
    // and  NumericVector::add_vector() members do this for us.
    matrix.add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
  }

  matrix.close();
}

#endif // defined(LIBMESH_HAVE_EIGEN_DENSE) && defined(LIBMESH_HAVE_PETSC)
