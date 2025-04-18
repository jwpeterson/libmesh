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



#ifndef LIBMESH_FE_TYPE_H
#define LIBMESH_FE_TYPE_H

// Local includes
#include "libmesh/compare_types.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/enum_order.h" // return Order objects
#include "libmesh/enum_fe_family.h" // LAGRANGE
#include "libmesh/enum_inf_map_type.h" // CARTESIAN
#include "libmesh/hashing.h"

// C++ includes
#include <memory>
#include <functional>

namespace libMesh
{

// Forward declarations
class QBase;
class FEType;

/**
 * This provides a shim class that wraps the Order enum.
 * The purpose of this is to store the order as an int
 * instead of an enum (to enable higher orders) while
 * retaining backwards compatibility.
 */
class OrderWrapper
{
public:

  /**
   * Constructor. Enables implicit conversion from an Order
   * enum to an OrderWrapper.
   */
  OrderWrapper(Order order) :
    _order(static_cast<int>(order))
  {}

  /**
   * Constructor. Enables implicit conversion from an int
   * to an OrderWrapper.
   */
  OrderWrapper(int order) :
    _order(order)
  {}

  /**
   * Operator that enables implicit conversion to
   * an Order enum.
   */
  operator Order() const
  {
    return static_cast<Order>(_order);
  }

  /**
   * Explicitly request the order as an int.
   */
  int get_order() const
  {
    return _order;
  }

  template <typename T>
  inline OrderWrapper & operator+=(T p)
  {
    _order += int(p);
    return *this;
  }

  template <typename T>
  inline OrderWrapper & operator-=(T p)
  {
    _order -= int(p);
    return *this;
  }

private:

  /**
   * The approximation order of the element.
   */
  int _order;

  friend struct std::hash<FEType>;
};

/**
 * Overload comparison operators for OrderWrapper.
 */
inline bool operator==(const OrderWrapper & lhs, const OrderWrapper & rhs){ return lhs.get_order() == rhs.get_order(); }
inline bool operator!=(const OrderWrapper & lhs, const OrderWrapper & rhs){ return !(lhs == rhs); }
inline bool operator< (const OrderWrapper & lhs, const OrderWrapper & rhs){ return lhs.get_order() < rhs.get_order(); }
inline bool operator> (const OrderWrapper & lhs, const OrderWrapper & rhs){ return rhs < lhs; }
inline bool operator<=(const OrderWrapper & lhs, const OrderWrapper & rhs){ return !(lhs > rhs); }
inline bool operator>=(const OrderWrapper & lhs, const OrderWrapper & rhs){ return !(lhs < rhs); }

// First disambiguate everything that would be ambiguated by the
// subsequent disambiguations
#define OrderWrapperOperators(comparisontype)           \
  inline bool operator==(comparisontype lhs, Order rhs) \
  { return lhs == static_cast<comparisontype>(rhs); }   \
  inline bool operator==(Order lhs, comparisontype rhs) \
  { return static_cast<comparisontype>(lhs) == rhs; }   \
  inline bool operator!=(comparisontype lhs, Order rhs) \
  { return !(lhs == rhs); }                             \
  inline bool operator!=(Order lhs, comparisontype rhs) \
  { return !(lhs == rhs); }                             \
  inline bool operator< (comparisontype lhs, Order rhs) \
  { return lhs < static_cast<comparisontype>(rhs); }    \
  inline bool operator< (Order lhs, comparisontype rhs) \
  { return static_cast<comparisontype>(lhs) < rhs; }    \
  inline bool operator> (comparisontype lhs, Order rhs) \
  { return rhs < lhs; }                                 \
  inline bool operator> (Order lhs, comparisontype rhs) \
  { return rhs < lhs; }                                 \
  inline bool operator<=(comparisontype lhs, Order rhs) \
  { return !(lhs > rhs); }                              \
  inline bool operator<=(Order lhs, comparisontype rhs) \
  { return !(lhs > rhs); }                              \
  inline bool operator>=(comparisontype lhs, Order rhs) \
  { return !(lhs < rhs); }                              \
  inline bool operator>=(Order lhs, comparisontype rhs) \
  { return !(lhs < rhs); }

OrderWrapperOperators(int)
OrderWrapperOperators(unsigned int)
#if LIBMESH_SIZEOF_SIZE_T != LIBMESH_SIZEOF_UNSIGNED_INT
OrderWrapperOperators(std::size_t)
#endif

// Now disambiguate all the things
inline bool operator==(int lhs, const OrderWrapper & rhs){ return lhs == rhs.get_order(); }
inline bool operator==(const OrderWrapper & lhs, int rhs){ return lhs.get_order() == rhs; }
inline bool operator==(Order lhs, const OrderWrapper & rhs){ return lhs == rhs.get_order(); }
inline bool operator==(const OrderWrapper & lhs, Order rhs){ return lhs.get_order() == rhs; }
inline bool operator!=(int lhs, const OrderWrapper & rhs){ return !(lhs == rhs); }
inline bool operator!=(const OrderWrapper & lhs, int rhs){ return !(lhs == rhs); }
inline bool operator!=(Order lhs, const OrderWrapper & rhs){ return !(lhs == rhs); }
inline bool operator!=(const OrderWrapper & lhs, Order rhs){ return !(lhs == rhs); }
inline bool operator< (int lhs, const OrderWrapper & rhs){ return lhs < rhs.get_order(); }
inline bool operator< (const OrderWrapper & lhs, int rhs){ return lhs.get_order() < rhs; }
inline bool operator< (Order lhs, const OrderWrapper & rhs){ return lhs < rhs.get_order(); }
inline bool operator< (const OrderWrapper & lhs, Order rhs){ return lhs.get_order() < rhs; }
inline bool operator> (int lhs, const OrderWrapper & rhs){ return rhs < lhs; }
inline bool operator> (const OrderWrapper & lhs, int rhs){ return rhs < lhs; }
inline bool operator> (Order lhs, const OrderWrapper & rhs){ return rhs < lhs; }
inline bool operator> (const OrderWrapper & lhs, Order rhs){ return rhs < lhs; }
inline bool operator<=(int lhs, const OrderWrapper & rhs){ return !(lhs > rhs); }
inline bool operator<=(const OrderWrapper & lhs, int rhs){ return !(lhs > rhs); }
inline bool operator<=(Order lhs, const OrderWrapper & rhs){ return !(lhs > rhs); }
inline bool operator<=(const OrderWrapper & lhs, Order rhs){ return !(lhs > rhs); }
inline bool operator>=(int lhs, const OrderWrapper & rhs){ return !(lhs < rhs); }
inline bool operator>=(const OrderWrapper & lhs, int rhs){ return !(lhs < rhs); }
inline bool operator>=(Order lhs, const OrderWrapper & rhs){ return !(lhs < rhs); }
inline bool operator>=(const OrderWrapper & lhs, Order rhs){ return !(lhs < rhs); }

/**
 * Overload stream operators.
 */
inline std::ostream & operator << (std::ostream & os, const OrderWrapper & order)
{
  os << order.get_order();
  return os;
}

/**
 * class FEType hides (possibly multiple) FEFamily and approximation
 * orders, thereby enabling specialized finite element families.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Manages the family, order, etc. parameters for a given FE.
 */
class FEType
{
public:

#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * Constructor.  Optionally takes the approximation \p Order
   * and the finite element family \p FEFamily
   */
  FEType(const int      o = 1,
         const FEFamily f = LAGRANGE) :
    order(o),
    family(f)
  {}

  /**
   * The approximation order of the element.
   */
  OrderWrapper order;

  /**
   * The type of finite element.  Valid types are \p LAGRANGE,
   * \p HIERARCHIC, etc...
   */
  FEFamily family;

#else

  /**
   * Constructor.  Optionally takes the approximation \p Order
   * and the finite element family \p FEFamily.
   *
   * \note For non-infinite elements, the \p order and \p base order
   * are the same, as with the \p family and \p base_family.  It must
   * be so, otherwise what we switch on would change when infinite
   * elements are not compiled in.
   */
  FEType(const int        o  = 1,
         const FEFamily   f  = LAGRANGE,
         const int        ro = THIRD,
         const FEFamily   rf = JACOBI_20_00,
         const InfMapType im = CARTESIAN) :
    order(o),
    radial_order(ro),
    family(f),
    radial_family(rf),
    inf_map(im)
  {}

  /**
   * The approximation order in the base of the infinite element.
   */
  OrderWrapper order;

  /**
   * The approximation order in radial direction of the infinite element.
   */
  OrderWrapper radial_order;

  /**
   * For InfFE, \p family contains the radial shape family, while
   * \p base_family contains the approximation type in circumferential
   * direction.  Valid types are \p LAGRANGE, \p HIERARCHIC, etc...
   */
  FEFamily family;

  /**
   * The type of approximation in radial direction.  Valid types are
   * \p JACOBI_20_00, \p JACOBI_30_00, etc...
   */
  FEFamily radial_family;

  /**
   * The coordinate mapping type of the infinite element.
   * When the infinite elements are defined over a surface with
   * a separable coordinate system (sphere, spheroid, ellipsoid),
   * the infinite elements may take advantage of this fact.
   */
  InfMapType inf_map;

#endif // ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * Tests equality
   */
  bool operator== (const FEType & f2) const
  {
    return (order == f2.order
            && family == f2.family
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
            && radial_order == f2.radial_order
            && radial_family == f2.radial_family
            && inf_map == f2.inf_map
#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
            );
  }

  /**
   * Tests inequality
   */
  bool operator!= (const FEType & f2) const
  {
    return !(*this == f2);
  }

  /**
   * An ordering to make FEType useful as a std::map key
   */
  bool operator< (const FEType & f2) const
  {
    if (order != f2.order)
      return (order < f2.order);
    if (family != f2.family)
      return (family < f2.family);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    if (radial_order != f2.radial_order)
      return (radial_order < f2.radial_order);
    if (radial_family != f2.radial_family)
      return (radial_family < f2.radial_family);
    if (inf_map != f2.inf_map)
      return (inf_map < f2.inf_map);
#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    return false;
  }

  /**
   * \returns The default quadrature order for this \p FEType.  The
   * default quadrature order is calculated assuming a polynomial of
   * degree \p order and is based on integrating the mass matrix for
   * such an element exactly on affine elements.
   */
  Order default_quadrature_order () const;

  /**
   * \returns A quadrature rule of appropriate type and order for this \p
   * FEType.  The default quadrature rule is based on integrating the mass
   * matrix for such an element exactly, with an additional power on
   * the basis order to help account for nonlinearities and/or
   * nonuniform coefficients.  Higher or lower degree rules can be
   * chosen by changing the extraorder parameter.
   */
  std::unique_ptr<QBase> default_quadrature_rule (const unsigned int dim,
                                                  const int extraorder=0) const;

  /**
   * \returns The default quadrature order for integrating unweighted
   * basis functions of this \p FEType.
   * The unweighted quadrature order is calculated assuming a
   * polynomial of degree \p order and is based on integrating the
   * shape functions for such an element exactly on affine elements.
   */
  Order unweighted_quadrature_order () const;

  /**
   * \returns A quadrature rule of appropriate type and order for
   * unweighted integration of this \p FEType.  The default quadrature
   * rule is based on integrating the shape functions on an affine
   * element exactly.  Higher or lower degree rules can be chosen by
   * changing the extraorder parameter.
   */
  std::unique_ptr<QBase> unweighted_quadrature_rule (const unsigned int dim,
                                                     const int extraorder=0) const;


private:

};



//-------------------------------------------------------------------
// FEType inline methods
inline
Order FEType::default_quadrature_order () const
{
  return static_cast<Order>(2*static_cast<unsigned int>(order.get_order()) + 1);
}

inline
Order FEType::unweighted_quadrature_order () const
{
  return order;
}

} // namespace libMesh

namespace std
{
template <>
struct hash<libMesh::FEType>
{
  std::size_t operator()(const libMesh::FEType & fe_type) const
    {
      std::size_t seed = 0;
      // Old compiler versions seem to need the static_cast
      libMesh::boostcopy::hash_combine(seed, static_cast<int>(fe_type.family));
      libMesh::boostcopy::hash_combine(seed, fe_type.order._order);
      return seed;
    }
};
}

#endif // LIBMESH_FE_TYPE_H
