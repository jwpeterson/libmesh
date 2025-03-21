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



#ifndef LIBMESH_SYSTEM_NORM_H
#define LIBMESH_SYSTEM_NORM_H

// Local includes
#include "libmesh/libmesh_common.h" // for Real

// C++ includes
#include <vector>

namespace libMesh
{

// Forward declarations
enum FEMNormType : int;

/**
 * This class defines a norm/seminorm to be applied to a NumericVector which
 * contains coefficients in a finite element space.
 *
 * Discrete vector norms and weighted l2 combinations of Sobolev norms and
 * seminorms are representable.
 *
 * For ease of use, if a norm or weight is queried for variable n but
 * has not been set for variable n, then the norm or weight set for
 * the highest-numbered variable below n is returned.
 *
 * \author Roy H. Stogner
 * \date 2008
 */
class SystemNorm
{
public:

  /**
   * Constructor, defaults to DISCRETE_L2 with weight of 1.0.
   */
  SystemNorm();

  /**
   * Constructor, for discrete vector norms, systems with one variable,
   * and systems for which the same norm type should be used with a
   * weight of one on each variable.
   *
   * This is deliberately an implicit constructor; we want user code
   * to be able to include lines like "error_norm = L2"
   */
  SystemNorm(const FEMNormType & t);

  /**
   * Constructor, for unweighted sobolev norms on systems with multiple
   * variables.
   *
   * For a system with n variables, the final norm will be the l2 norm of the
   * n-vector of the norms in each variable.
   */
  explicit
  SystemNorm(std::vector<FEMNormType> norms);

  /**
   * Constructor, for weighted sobolev norms on systems with multiple
   * variables.
   *
   * For a system with n variables, the final norm will be the l2 norm of the
   * n-vector of the norms in each variable, each multiplied by weight.
   */
  SystemNorm(std::vector<FEMNormType> norms,
             std::vector<Real> & weights);

  /**
   * Constructor, for weighted sobolev norms on systems with multiple
   * variables and their adjoints
   *
   * For a system with n variables, the final norm computed will be of the form
   * norm_u^T*R*norm_z where R is a scaling matrix
   */
  SystemNorm(std::vector<FEMNormType> norms,
             std::vector<std::vector<Real>> & weights);

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  SystemNorm (const SystemNorm &) = default;
  SystemNorm (SystemNorm &&) = default;
  SystemNorm & operator= (const SystemNorm &) = default;
  SystemNorm & operator= (SystemNorm &&) = default;
  virtual ~SystemNorm() = default;

  /**
   * \returns \p true if this is purely a discrete norm
   */
  bool is_discrete() const;

  /**
   * \returns The weighted norm v^T*W*v where W represents our
   * weights matrix or weights vector times identity matrix.
   */
  Real calculate_norm(const std::vector<Real> & v);

  /**
   * \returns The weighted inner product v1^T*W*v2 where R is our weights
   */
  Real calculate_norm(const std::vector<Real> & v1,
                      const std::vector<Real> & v2);

  /**
   * \returns \p true if no weight matrix W is specified or an identity matrix is specified, otherwise returns false
   */
  bool is_identity();

  /**
   * \returns The type of the norm in variable \p var
   *
   * If no norm has been explicitly set for \p var, then the
   * highest-index norm explicitly set is returned, or DISCRETE_L2
   * is returned if no norms have been explicitly set.
   */
  FEMNormType type(unsigned int var) const;

  /**
   * Sets the type of the norm in variable \p var, as well as for any
   * unset variables with index less than \p var
   */
  void set_type(unsigned int var, const FEMNormType & t);

  /**
   * \returns The weight corresponding to the norm in variable \p var
   *
   * If no weight has been explicitly set for \p var, then the
   * highest-index weight explicitly set is returned, or 1.0 is
   * returned if no weights have been explicitly set.
   */
  Real weight(unsigned int var) const;

  /**
   * Sets the weight corresponding to the norm in variable \p var, as
   * well as for any unset variables with index less than \p var.
   */
  void set_weight(unsigned int var, Real w);

  /**
   * Sets the weight corresponding to the norm from the variable pair v1(var1) coming from v2(var2). See calculate_norm
   */
  void set_off_diagonal_weight(unsigned int i, unsigned int j, Real w);

  /**
   * \returns The squared weight corresponding to the norm in variable
   * \p var.  We cache that at construction time to save a few flops.
   */
  Real weight_sq(unsigned int var) const;



private:
  std::vector<FEMNormType> _norms;

  std::vector<Real> _weights;
  std::vector<Real> _weights_sq;

  /**
   * One more data structure needed to store the off diagonal
   * components for the generalize SystemNorm case
   */
  std::vector<std::vector<Real>> _off_diagonal_weights;
};

} // namespace libMesh

#endif // LIBMESH_SYSTEM_NORM_H
