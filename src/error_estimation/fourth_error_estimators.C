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


#include "libmesh/libmesh_config.h"

#include "libmesh/fourth_error_estimators.h"

#include "libmesh/enum_error_estimator_type.h"

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/error_vector.h"
#include "libmesh/fe_base.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/system.h"
#include "libmesh/dense_vector.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/enum_norm_type.h"

// C++ includes
#include <algorithm> // for std::fill
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>    // for sqrt

namespace libMesh
{


LaplacianErrorEstimator::LaplacianErrorEstimator() :
  JumpErrorEstimator()
{
  error_norm = H2_SEMINORM;
}



ErrorEstimatorType
LaplacianErrorEstimator::type() const
{
  return LAPLACIAN;
}



void
LaplacianErrorEstimator::init_context(FEMContext & c)
{
  const unsigned int n_vars = c.n_vars();
  for (unsigned int v=0; v<n_vars; v++)
    {
      // Possibly skip this variable
      if (error_norm.weight(v) == 0.0) continue;

      // FIXME: Need to generalize this to vector-valued elements. [PB]
      FEBase * side_fe = nullptr;

      const std::set<unsigned char> & elem_dims =
        c.elem_dimensions();

      for (const auto & dim : elem_dims)
        {
          c.get_side_fe( v, side_fe, dim );

          // We'll need hessians on both sides for flux jump computation
          side_fe->get_d2phi();
        }
    }
}



void
LaplacianErrorEstimator::internal_side_integration ()
{
  const Elem & coarse_elem = coarse_context->get_elem();
  const Elem & fine_elem   = fine_context->get_elem();

  const DenseVector<Number> & Ucoarse = coarse_context->get_elem_solution();
  const DenseVector<Number> & Ufine   = fine_context->get_elem_solution();

  unsigned short dim = fine_elem.dim();

  FEBase * fe_fine = nullptr;
  fine_context->get_side_fe( var, fe_fine, dim );

  FEBase * fe_coarse = nullptr;
  coarse_context->get_side_fe( var, fe_coarse, dim );

  Real error = 1.e-30;
  unsigned int n_qp = fe_fine->n_quadrature_points();

  const std::vector<std::vector<RealTensor>> & d2phi_coarse = fe_coarse->get_d2phi();
  const std::vector<std::vector<RealTensor>> & d2phi_fine = fe_fine->get_d2phi();
  const std::vector<Real> & JxW_face = fe_fine->get_JxW();

  for (unsigned int qp=0; qp != n_qp; ++qp)
    {
      // Calculate solution gradients on fine and coarse elements
      // at this quadrature point
      Number laplacian_fine = 0., laplacian_coarse = 0.;

      const unsigned int n_coarse_dofs = Ucoarse.size();
      for (unsigned int i=0; i != n_coarse_dofs; ++i)
        {
          laplacian_coarse += d2phi_coarse[i][qp](0,0) * Ucoarse(i);
          if (dim > 1)
            laplacian_coarse += d2phi_coarse[i][qp](1,1) * Ucoarse(i);
          if (dim > 2)
            laplacian_coarse += d2phi_coarse[i][qp](2,2) * Ucoarse(i);
        }

      const unsigned int n_fine_dofs = Ufine.size();
      for (unsigned int i=0; i != n_fine_dofs; ++i)
        {
          laplacian_fine += d2phi_fine[i][qp](0,0) * Ufine(i);
          if (dim > 1)
            laplacian_fine += d2phi_fine[i][qp](1,1) * Ufine(i);
          if (dim > 2)
            laplacian_fine += d2phi_fine[i][qp](2,2) * Ufine(i);
        }


      // Find the jump in the Laplacian
      // at this quadrature point
      const Number jump = laplacian_fine - laplacian_coarse;
      const Real jump2 = TensorTools::norm_sq(jump);

      // Accumulate the jump integral
      error += JxW_face[qp] * jump2;
    }

  // Add the h-weighted jump integral to each error term
  fine_error =
    error * fine_elem.hmax() * error_norm.weight(var);
  coarse_error =
    error * coarse_elem.hmax() * error_norm.weight(var);
}

} // namespace libMesh

#else // defined (LIBMESH_ENABLE_SECOND_DERIVATIVES)

namespace libMesh
{

LaplacianErrorEstimator::LaplacianErrorEstimator() :
  JumpErrorEstimator()
{
  libmesh_error_msg("Error: LaplacianErrorEstimator requires second " \
                    << "derivative support; try configuring libmesh with " \
                    << "--enable-second");
}

ErrorEstimatorType
LaplacianErrorEstimator::type() const
{
  return LAPLACIAN;
}



void LaplacianErrorEstimator::init_context (FEMContext &) {}

void LaplacianErrorEstimator::internal_side_integration () {}

} // namespace libMesh

#endif // defined (LIBMESH_ENABLE_SECOND_DERIVATIVES)
