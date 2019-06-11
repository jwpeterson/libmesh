// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/quadrature_nodal.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/quadrature_simpson.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

void QNodal::init_3D(const ElemType, unsigned int)
{
#if LIBMESH_DIM == 3

  switch (_type)
    {
    case TET4:
    case PRISM6:
    case HEX8:
      {
        QTrap rule(/*dim=*/3, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }

    case PYRAMID5:
      {
        // A rule with 5 points which is exact only for constants.
        // Note: it is not possible to construct a rule with 5 points
        // at the vertices which is exact for linears, so the weights
        // are instead chosen as described below.  The quadrature
        // points are numbered the same way as the reference element
        // nodes.
        _points =
          {
            Point(-1,-1,0.), Point(+1,-1,0.), Point(+1,+1,0.), Point(-1,+1,0.),
            Point(0.,0.,+1)
          };

        // base vertex (wb), and apex vertex (wa) weights are obtained by:
        // 1.) Requiring that they sum to the reference element volume.
        // 2.) Minimizing the Frobenius norm of the difference between
        //     the resulting nodal quadrature (diagonal) mass matrix
        //     and the true mass matrix for the reference element.
        Real wb = Real(58) / 225;
        Real wa = Real(68) / 225;

        _weights = {wb, wb, wb, wb, wa};

        return;
      }

    case PRISM15:
      {
        // A rule with 15 points which is exact for linears, and
        // naturally produces a lumped approximation to the mass
        // matrix. The quadrature points are numbered the same way as
        // the reference element nodes.
        _points =
          {
            Point(0.,0.,-1), Point(+1,0.,-1), Point(0.,+1,-1),
            Point(0.,0.,+1), Point(+1,0.,+1), Point(0.,+1,+1),
            Point(.5,0.,-1), Point(.5,.5,-1), Point(0.,.5,-1),
            Point(0.,0.,0.), Point(+1,0.,0.), Point(0.,+1,0.),
            Point(.5,0.,+1), Point(.5,.5,+1), Point(0.,.5,+1),
          };

        // vertex (wv), tri edge (wt), and quad edge (wq) weights are
        // obtained by:
        // 1.) Requiring that they sum to the reference element volume.
        // 2.) Minimizing the Frobenius norm of the difference between
        //     the resulting nodal quadrature (diagonal) mass matrix
        //     and the true mass matrix for the reference element.
        Real wv = Real(26) / 675;
        Real wt = Real(17) / 225;
        Real wq = Real(71) / 675;

        _weights = {wv, wv, wv, wv, wv, wv,
                    wt, wt, wt,
                    wq, wq, wq,
                    wt, wt, wt};

        return;
      }

    case PYRAMID14:
      {
        // A rule with 14 points which is exact for linears, and
        // naturally produces a lumped approximation to the mass
        // matrix. The quadrature points are numbered the same way as
        // the reference element nodes.
        _points =
          {
            Point(-1.,-1., 0.), Point( 1.,-1., 0.), Point( 1., 1., 0.), Point(-1., 1., 0.),
            Point( 0., 0., 1.),
            Point( 0.,-1., 0.), Point( 1., 0., 0.), Point( 0., 1., 0.), Point(-1., 0., 0.),
            Point(-.5,-.5, .5), Point( .5,-.5, .5), Point( .5, .5, .5), Point(-.5, .5, .5),
            Point( 0., 0., 0.)
          };

        // wb (4) base (vertex) weights
        // wa (1) apex (vertex) weight
        // wq (4) quad (edge) weights
        // we (4) tri (edge) weights
        // wc (1) base centroid weight
        Real wb = Real(4726)  / 89775; // ~0.05264
        Real wa = Real(7783)  / 89775; // ~0.08669
        Real wq = Real(2519)  / 29925; // ~0.08418
        Real we = Real(11071) / 89775; // ~0.12332
        Real wc = Real(881)   / 4275;  // ~0.20608

        _weights = {wb, wb, wb, wb,
                    wa,
                    wq, wq, wq, wq,
                    we, we, we, we,
                    wc};

        return;
      }

    case HEX20:
      {
        // A rule with 20 points which is exact for linears, and
        // naturally produces a lumped approximation to the mass
        // matrix. The quadrature points are numbered the same way as
        // the reference element nodes.
        _points =
          {
            Point(-1,-1,-1), Point(+1,-1,-1), Point(+1,+1,-1), Point(-1,+1,-1),
            Point(-1,-1,+1), Point(+1,-1,+1), Point(+1,+1,+1), Point(-1,+1,+1),
            Point(0.,-1,-1), Point(+1,0.,-1), Point(0.,+1,-1), Point(-1,0.,-1),
            Point(-1,-1,0.), Point(+1,-1,0.), Point(+1,+1,0.), Point(-1,+1,0.),
            Point(0.,-1,+1), Point(+1,0.,+1), Point(0.,+1,+1), Point(-1,0.,+1)
          };

        // vertex (wv), and edge (we) weights are obtained by:
        // 1.) Requiring that they sum to the reference element volume.
        // 2.) Minimizing the Frobenius norm of the difference between
        //     the resulting nodal quadrature (diagonal) mass matrix
        //     and the true mass matrix for the reference element.
        Real wv = Real(136) / 585;
        Real we = Real(898) / 1755;

        _weights = {wv, wv, wv, wv, wv, wv, wv, wv,
                    we, we, we, we, we, we, we, we, we, we, we, we};

        return;
      }

    case TET10:
    case PRISM18:
    case HEX27:
      {
        QSimpson rule(/*dim=*/3, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }

    default:
      libmesh_error_msg("ERROR: Unsupported type: " << Utility::enum_to_string(_type));
    }
#endif
}

} // namespace libMesh
