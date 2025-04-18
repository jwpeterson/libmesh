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

#include "libmesh/fe_transformation_base.h"
#include "libmesh/h1_fe_transformation.h"
#include "libmesh/hcurl_fe_transformation.h"
#include "libmesh/hdiv_fe_transformation.h"
#include "libmesh/fe_type.h"
#include "libmesh/enum_to_string.h"

// C++ Includes
#include <memory>

namespace libMesh
{

template<typename OutputShape>
std::unique_ptr<FETransformationBase<OutputShape>> FETransformationBase<OutputShape>::build( const FEType & fe_type )
{
  switch (fe_type.family)
    {
      // H1 Conforming Elements
    case LAGRANGE:
    case HIERARCHIC:
    case BERNSTEIN:
    case SZABAB:
    case CLOUGH: // PB: Really H2
    case HERMITE: // PB: Really H2
    case SUBDIVISION:
    case LAGRANGE_VEC:
    case HIERARCHIC_VEC:
    case MONOMIAL: // PB: Shouldn't this be L2 conforming?
    case MONOMIAL_VEC: // PB: Shouldn't this be L2 conforming?
    case XYZ: // PB: Shouldn't this be L2 conforming?
    case RATIONAL_BERNSTEIN:
    case L2_HIERARCHIC:
    case L2_HIERARCHIC_VEC:
    case SIDE_HIERARCHIC:
    case L2_LAGRANGE: // PB: Shouldn't this be L2 conforming?
    case L2_LAGRANGE_VEC: // PB: Shouldn't this be L2 conforming?
    case JACOBI_20_00: // PB: For infinite elements...
    case JACOBI_30_00: // PB: For infinite elements...
      return std::make_unique<H1FETransformation<OutputShape>>();

      // HCurl Conforming Elements
    case NEDELEC_ONE:
      return std::make_unique<HCurlFETransformation<OutputShape>>();

      // HDiv Conforming Elements
    case RAVIART_THOMAS:
    case L2_RAVIART_THOMAS:
      return std::make_unique<HDivFETransformation<OutputShape>>();

      // L2 Conforming Elements

      // Other...
    case SCALAR:
      // Should never need this for SCALARs
      return std::make_unique<H1FETransformation<OutputShape>>();

    default:
      libmesh_error_msg("Unknown family = " << Utility::enum_to_string(fe_type.family));
    }
}

template class LIBMESH_EXPORT FETransformationBase<Real>;
template class LIBMESH_EXPORT FETransformationBase<RealGradient>;

} // namespace libMesh
