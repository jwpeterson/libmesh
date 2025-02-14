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



#ifndef LIBMESH_UNIFORM_REFINEMENT_ESTIMATOR_H
#define LIBMESH_UNIFORM_REFINEMENT_ESTIMATOR_H

// Local Includes
#include "libmesh/error_estimator.h"
#include "libmesh/libmesh.h"

// C++ includes
#include <cstddef>
#include <vector>

#ifdef LIBMESH_ENABLE_AMR

namespace libMesh
{

/**
 * This class implements a ``brute force'' error estimator
 * which integrates differences between the current solution
 * and the solution on a uniformly refined (in h and/or p,
 * for an arbitrary number of levels) grid.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class UniformRefinementEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Sets the most common default parameter values.
   */
  UniformRefinementEstimator();

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  UniformRefinementEstimator (const UniformRefinementEstimator &) = default;
  UniformRefinementEstimator (UniformRefinementEstimator &&) = default;
  UniformRefinementEstimator & operator= (const UniformRefinementEstimator &) = default;
  UniformRefinementEstimator & operator= (UniformRefinementEstimator &&) = default;
  virtual ~UniformRefinementEstimator() = default;

  /**
   * This function does uniform refinements and a
   * solve to get an improved solution on each cell,
   * then estimates the error by integrating differences
   * between the coarse and fine solutions.
   *
   * system.solve() must be called, and so should
   * have no side effects.
   *
   * Only the provided system is solved on the refined mesh;
   * for problems decoupled into multiple systems, use of
   * estimate_errors() should be more reliable.
   *
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System & system,
                               ErrorVector & error_per_cell,
                               const NumericVector<Number> * solution_vector = nullptr,
                               bool estimate_parent_error = false) override;

  /**
   * Currently this function ignores the error_norm member variable,
   * and uses the function argument error_norms instead.
   *
   * This function is named estimate_errors instead of estimate_error
   * because otherwise C++ can get confused.
   */
  virtual void estimate_errors (const EquationSystems & equation_systems,
                                ErrorVector & error_per_cell,
                                const std::map<const System *, SystemNorm> & error_norms,
                                const std::map<const System *, const NumericVector<Number> *> * solution_vectors = nullptr,
                                bool estimate_parent_error = false) override;

  /**
   * Currently this function ignores the component_scale member variable,
   * because it calculates each error individually, unscaled.
   *
   * The user selects which errors get computed by filling a map with error
   * vectors: If errors_per_cell[&system][v] exists, it will be filled with the
   * error values in variable \p v of \p system
   */
  virtual void estimate_errors (const EquationSystems & equation_systems,
                                ErrorMap & errors_per_cell,
                                const std::map<const System *, const NumericVector<Number> *> * solution_vectors = nullptr,
                                bool estimate_parent_error = false) override;

  /**
   * Increases or decreases the order of the quadrature rule used for numerical
   * integration.  The default \p extraorder is 1, because properly
   * integrating L2 error requires integrating the squares of terms
   * with order p+1, and 2p+2 is 1 higher than what we default to
   * using for reasonable mass matrix integration.
   */
  void extra_quadrature_order (const int extraorder)
  { _extra_order = extraorder; }

  virtual ErrorEstimatorType type() const override;

  /**
   * How many h refinements to perform to get the fine grid
   */
  unsigned char number_h_refinements;

  /**
   * How many p refinements to perform to get the fine grid
   */
  unsigned char number_p_refinements;

protected:

  /**
   * Extra order to use for quadrature rule
   */
  int _extra_order;

  /**
   * The code for estimate_error and both estimate_errors versions is very
   * similar, so we use the same function for all three
   */
  virtual void _estimate_error (const EquationSystems * equation_systems,
                                const System * system,
                                ErrorVector * error_per_cell,
                                ErrorMap * errors_per_cell,
                                const std::map<const System *, SystemNorm > * error_norms,
                                const std::map<const System *, const NumericVector<Number> *> * solution_vectors = nullptr,
                                bool estimate_parent_error = false);
};

} // namespace libMesh

#endif // #ifdef LIBMESH_ENABLE_AMR

#endif // LIBMESH_UNIFORM_REFINEMENT_ESTIMATOR_H
