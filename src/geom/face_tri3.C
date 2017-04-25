// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes

// Local includes
#include "libmesh/side.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"

namespace {
/**
 * The following function is from:
 *
 * Routines for Arbitrary Precision Floating-point Arithmetic
 * and Fast Robust Geometric Predicates
 * (predicates.c)
 *
 * May 18, 1996
 *
 * Placed in the public domain by
 * Jonathan Richard Shewchuk
 * School of Computer Science
 * Carnegie Mellon University
 * 5000 Forbes Avenue
 * Pittsburgh, Pennsylvania  15213-3891
 * jrs@cs.cmu.edu
 */
double orient3d(const double * pa,
                const double * pb,
                const double * pc,
                const double * pd);
}

namespace libMesh
{



// ------------------------------------------------------------
// Tri3 class static member initializations
const unsigned int Tri3::side_nodes_map[3][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 0}  // Side 2
  };


#ifdef LIBMESH_ENABLE_AMR

const float Tri3::_embedding_matrix[4][3][3] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2
      {1.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0}, // 1
      {0.5, 0.0, 0.5}  // 2
    },

    // embedding matrix for child 1
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 1.0, 0.0}, // 1
      {0.0, 0.5, 0.5}  // 2
    },

    // embedding matrix for child 2
    {
      // 0    1    2
      {0.5, 0.0, 0.5}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.0, 0.0, 1.0}  // 2
    },

    // embedding matrix for child 3
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.5, 0.0, 0.5}  // 2
    }
  };

#endif



// ------------------------------------------------------------
// Tri3 class member functions

bool Tri3::is_vertex(const unsigned int) const
{
  return true;
}

bool Tri3::is_edge(const unsigned int) const
{
  return false;
}

bool Tri3::is_face(const unsigned int) const
{
  return false;
}

bool Tri3::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 2; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

UniquePtr<Elem> Tri3::build_side_ptr (const unsigned int i,
                                      bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return UniquePtr<Elem>(new Side<Edge2,Tri3>(this,i));

  else
    {
      Elem * edge = new Edge2;
      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<edge->n_nodes(); ++n)
        edge->set_node(n) = this->node_ptr(Tri3::side_nodes_map[i][n]);

      return UniquePtr<Elem>(edge);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}


void Tri3::connectivity(const unsigned int libmesh_dbg_var(sf),
                        const IOPackage iop,
                        std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(4);
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(2)+1;
        return;
      }

    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}






Real Tri3::volume () const
{
  // 3-node triangles have the following formula for computing the area
  return 0.5 * cross_norm(point(1) - point(0),
                          point(2) - point(0));
}



std::pair<Real, Real> Tri3::min_and_max_angle() const
{
  Point v10 ( this->point(1) - this->point(0) );
  Point v20 ( this->point(2) - this->point(0) );
  Point v21 ( this->point(2) - this->point(1) );

  const Real
    len_10=v10.norm(),
    len_20=v20.norm(),
    len_21=v21.norm();

  const Real
    theta0=std::acos(( v10*v20)/len_10/len_20),
    theta1=std::acos((-v10*v21)/len_10/len_21),
    theta2=libMesh::pi - theta0 - theta1
    ;

  libmesh_assert_greater (theta0, 0.);
  libmesh_assert_greater (theta1, 0.);
  libmesh_assert_greater (theta2, 0.);

  return std::make_pair(std::min(theta0, std::min(theta1,theta2)),
                        std::max(theta0, std::max(theta1,theta2)));
}

bool Tri3::contains_point (const Point & p, Real tol) const
{
  // See "Barycentric Technique" section at
  // http://www.blackpawn.com/texts/pointinpoly for details.

  // Compute vectors
  Point v0 = this->point(1) - this->point(0);
  Point v1 = this->point(2) - this->point(0);
  Point v2 = p - this->point(0);

  // Compute dot products
  Real dot00 = v0 * v0;
  Real dot01 = v0 * v1;
  Real dot02 = v0 * v2;
  Real dot11 = v1 * v1;
  Real dot12 = v1 * v2;

  // Out of plane check
  if (std::abs(triple_product(v2, v0, v1)) / std::max(dot00, dot11) > tol)
    return false;

  // Compute barycentric coordinates
  Real invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  Real u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  Real v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // Check if point is in triangle
  return (u > -tol) && (v > -tol) && (u + v < 1 + tol);
}

} // namespace libMesh



namespace {

/**
 * Original documentation for these functions:
 *
 * orient3d() Adaptive exact 3D orientation test.  Robust.
 *
 *            Return a positive value if the point pd lies below the
 *            plane passing through pa, pb, and pc; "below" is defined so
 *            that pa, pb, and pc appear in counterclockwise order when
 *            viewed from above the plane.  Returns a negative value if
 *            pd lies above the plane.  Returns zero if the points are
 *            coplanar.  The result is also a rough approximation of six
 *            times the signed volume of the tetrahedron defined by the
 *            four points.
 *
 *  Uses exact arithmetic to ensure a correct answer.  The
 *  result returned is the determinant of a matrix.  In orient3d() only,
 *  this determinant is computed adaptively, in the sense that exact
 *  arithmetic is used only to the degree it is needed to ensure that the
 *  returned value has the correct sign.  Hence, orient3d() is usually quite
 *  fast, but will run more slowly when the input points are coplanar or
 *  nearly so.
 *
 * Notes on adapting the original function for use in libmesh
 * - Copied orient3dadapt(), orient3d() into an anonymous namespace.
 * - Deleted INEXACT macro.
 * - Replaced double macro with hardcoded double.
 * - Copied the Split, Two_Product, etc. macros from predicates.c.
 * - Copied constants from constants.C.
 * - Include <float.h>
 * - Reindented everything.
 */

// DBL_MANT_DIG, etc.
#include <float.h>

const double splitter = ((DBL_MANT_DIG + 1) >> 1) + 1.0;

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (double) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Split(a, ahi, alo)     \
  c = (double) (splitter * a); \
  abig = (double) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y)            \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (double) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y)                    \
  x = (double) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Two_Product(a, b, x, y) \
  x = (double) (a * b); \
  Two_Product_Tail(a, b, x, y)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0)    \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)


double orient3dadapt(const double * pa, const double * pb, const double * pc, const double * pd, double permanent)
{
  double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
  double det, errbound;

  double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
  double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
  double bc[4], ca[4], ab[4];
  double bc3, ca3, ab3;
  double adet[8], bdet[8], cdet[8];
  int alen, blen, clen;
  double abdet[16];
  int ablen;
  double *finnow, *finother, *finswap;
  double fin1[192], fin2[192];
  int finlength;

  double adxtail, bdxtail, cdxtail;
  double adytail, bdytail, cdytail;
  double adztail, bdztail, cdztail;
  double at_blarge, at_clarge;
  double bt_clarge, bt_alarge;
  double ct_alarge, ct_blarge;
  double at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
  int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
  double bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
  double adxt_cdy1, adxt_bdy1, bdxt_ady1;
  double bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
  double adxt_cdy0, adxt_bdy0, bdxt_ady0;
  double bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
  double adyt_cdx1, adyt_bdx1, bdyt_adx1;
  double bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
  double adyt_cdx0, adyt_bdx0, bdyt_adx0;
  double bct[8], cat[8], abt[8];
  int bctlen, catlen, abtlen;
  double bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
  double adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
  double bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
  double adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
  double u[4], v[12], w[16];
  double u3;
  int vlength, wlength;
  double negate;

  double bvirt;
  double avirt, bround, around;
  double c;
  double abig;
  double ahi, alo, bhi, blo;
  double err1, err2, err3;
  double _i, _j, _k;
  double _0;

  adx = (double) (pa[0] - pd[0]);
  bdx = (double) (pb[0] - pd[0]);
  cdx = (double) (pc[0] - pd[0]);
  ady = (double) (pa[1] - pd[1]);
  bdy = (double) (pb[1] - pd[1]);
  cdy = (double) (pc[1] - pd[1]);
  adz = (double) (pa[2] - pd[2]);
  bdz = (double) (pb[2] - pd[2]);
  cdz = (double) (pc[2] - pd[2]);

  Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
  Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
  Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
  bc[3] = bc3;
  alen = scale_expansion_zeroelim(4, bc, adz, adet);

  Two_Product(cdx, ady, cdxady1, cdxady0);
  Two_Product(adx, cdy, adxcdy1, adxcdy0);
  Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
  ca[3] = ca3;
  blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

  Two_Product(adx, bdy, adxbdy1, adxbdy0);
  Two_Product(bdx, ady, bdxady1, bdxady0);
  Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
  ab[3] = ab3;
  clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

  det = estimate(finlength, fin1);
  errbound = o3derrboundB * permanent;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
  Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
  Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
  Two_Diff_Tail(pa[1], pd[1], ady, adytail);
  Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
  Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
  Two_Diff_Tail(pa[2], pd[2], adz, adztail);
  Two_Diff_Tail(pb[2], pd[2], bdz, bdztail);
  Two_Diff_Tail(pc[2], pd[2], cdz, cdztail);

  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
      && (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
    return det;
  }

  errbound = o3derrboundC * permanent + resulterrbound * Absolute(det);
  det += (adz * ((bdx * cdytail + cdy * bdxtail)
                 - (bdy * cdxtail + cdx * bdytail))
          + adztail * (bdx * cdy - bdy * cdx))
    + (bdz * ((cdx * adytail + ady * cdxtail)
              - (cdy * adxtail + adx * cdytail))
       + bdztail * (cdx * ady - cdy * adx))
    + (cdz * ((adx * bdytail + bdy * adxtail)
              - (ady * bdxtail + bdx * adytail))
       + cdztail * (adx * bdy - ady * bdx));
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  finnow = fin1;
  finother = fin2;

  if (adxtail == 0.0) {
    if (adytail == 0.0) {
      at_b[0] = 0.0;
      at_blen = 1;
      at_c[0] = 0.0;
      at_clen = 1;
    } else {
      negate = -adytail;
      Two_Product(negate, bdx, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      Two_Product(adytail, cdx, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    }
  } else {
    if (adytail == 0.0) {
      Two_Product(adxtail, bdy, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      negate = -adxtail;
      Two_Product(negate, cdy, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    } else {
      Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
      Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
      Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
                   at_blarge, at_b[2], at_b[1], at_b[0]);
      at_b[3] = at_blarge;
      at_blen = 4;
      Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
      Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
      Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
                   at_clarge, at_c[2], at_c[1], at_c[0]);
      at_c[3] = at_clarge;
      at_clen = 4;
    }
  }
  if (bdxtail == 0.0) {
    if (bdytail == 0.0) {
      bt_c[0] = 0.0;
      bt_clen = 1;
      bt_a[0] = 0.0;
      bt_alen = 1;
    } else {
      negate = -bdytail;
      Two_Product(negate, cdx, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    }
  } else {
    if (bdytail == 0.0) {
      Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      negate = -bdxtail;
      Two_Product(negate, ady, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    } else {
      Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
      Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
      Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
                   bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
      bt_c[3] = bt_clarge;
      bt_clen = 4;
      Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
      Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
      Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
                   bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
      bt_a[3] = bt_alarge;
      bt_alen = 4;
    }
  }
  if (cdxtail == 0.0) {
    if (cdytail == 0.0) {
      ct_a[0] = 0.0;
      ct_alen = 1;
      ct_b[0] = 0.0;
      ct_blen = 1;
    } else {
      negate = -cdytail;
      Two_Product(negate, adx, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    }
  } else {
    if (cdytail == 0.0) {
      Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      negate = -cdxtail;
      Two_Product(negate, bdy, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    } else {
      Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
      Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
      Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
                   ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
      ct_a[3] = ct_alarge;
      ct_alen = 4;
      Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
      Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
      Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
                   ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
      ct_b[3] = ct_blarge;
      ct_blen = 4;
    }
  }

  bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
  wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
  wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
  wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  if (adztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, bc, adztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  if (adxtail != 0.0) {
    if (bdytail != 0.0) {
      Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
      Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (cdytail != 0.0) {
      negate = -adxtail;
      Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
      Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (bdxtail != 0.0) {
    if (cdytail != 0.0) {
      Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
      Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (adytail != 0.0) {
      negate = -bdxtail;
      Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
      Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (cdxtail != 0.0) {
    if (adytail != 0.0) {
      Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
      Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (bdytail != 0.0) {
      negate = -cdxtail;
      Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
      Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }

  if (adztail != 0.0) {
    wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  return finnow[finlength - 1];
} // orient3dadapt()


double orient3d(const double * pa, const double * pb, const double * pc, const double * pd)
{
  double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
  double bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
  double det;
  double permanent, errbound;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];
  adz = pa[2] - pd[2];
  bdz = pb[2] - pd[2];
  cdz = pc[2] - pd[2];

  bdxcdy = bdx * cdy;
  cdxbdy = cdx * bdy;

  cdxady = cdx * ady;
  adxcdy = adx * cdy;

  adxbdy = adx * bdy;
  bdxady = bdx * ady;

  det = adz * (bdxcdy - cdxbdy)
    + bdz * (cdxady - adxcdy)
    + cdz * (adxbdy - bdxady);

  permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz)
    + (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz)
    + (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz);
  errbound = o3derrboundA * permanent;
  if ((det > errbound) || (-det > errbound)) {
    return det;
  }

  return orient3dadapt(pa, pb, pc, pd, permanent);
}

} // anonymous namespace
