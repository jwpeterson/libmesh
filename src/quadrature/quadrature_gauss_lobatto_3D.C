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



// Local includes
#include "libmesh/quadrature_gauss_lobatto.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/quadrature_gauss.h"

namespace libMesh
{

void QGaussLobatto::init_3D()
{
  switch (_type)
    {
    case HEX8:
    case HEX20:
    case HEX27:
      {
        // We compute the 3D quadrature rule as a tensor
        // product of the 1D quadrature rule.
        QGaussLobatto q1D(1, get_order());
        tensor_product_hex(q1D);
        return;
      }

      // We fall back on a Gauss type rule for other types of elements,
      // but we warn the user (once) that we are doing this instead of
      // silently switching out quadrature rules on them.
    default:
      {
        libmesh_warning("Warning: QGaussLobatto falling back on QGauss rule "
                        "for unsupported Elem type: " << Utility::enum_to_string(_type));

        QGauss gauss_rule(_dim, _order);
        gauss_rule.init(*this);

        // Swap points and weights with the about-to-be destroyed rule.
        _points.swap (gauss_rule.get_points() );
        _weights.swap(gauss_rule.get_weights());
      }
    }
}

} // namespace libMesh
