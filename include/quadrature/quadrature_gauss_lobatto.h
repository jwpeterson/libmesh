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



#ifndef LIBMESH_QUADRATURE_GAUSS_LOBATTO_H
#define LIBMESH_QUADRATURE_GAUSS_LOBATTO_H

// Local includes
#include "libmesh/quadrature.h"

namespace libMesh
{

/**
 * This class implements Gauss-Lobatto quadrature for 1D elements and 2D/3D
 * tensor product elements.  Properties of Gauss-Lobatto quadrature rules:
 * .) Include the "end-points" of the domain (have points on edges/faces in 2D/3D).
 * .) Rules with n points can exactly integrate polynomials of degree 2n-3.
 *
 * http://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Lobatto_rules
 *
 * \author John W. Peterson
 * \date 2014
 * \brief Implements 1D and 2/3D tensor product Gauss-Lobatto quadrature rules.
 */
class QGaussLobatto final : public QBase
{
public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGaussLobatto (unsigned int dim,
                 Order order=INVALID_ORDER);

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  QGaussLobatto (const QGaussLobatto &) = default;
  QGaussLobatto (QGaussLobatto &&) = default;
  QGaussLobatto & operator= (const QGaussLobatto &) = default;
  QGaussLobatto & operator= (QGaussLobatto &&) = default;
  virtual ~QGaussLobatto() = default;

  /**
   * \returns \p QGAUSS_LOBATTO.
   */
  virtual QuadratureType type() const override;

  virtual std::unique_ptr<QBase> clone() const override;

private:

  virtual void init_1D () override;
  virtual void init_2D () override;
  virtual void init_3D () override;
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_GAUSS_LOBATTO_H
