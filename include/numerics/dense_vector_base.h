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



#ifndef LIBMESH_DENSE_VECTOR_BASE_H
#define LIBMESH_DENSE_VECTOR_BASE_H


// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/int_range.h"

// C++ includes

namespace libMesh
{

/**
 * Defines an abstract dense vector base class for use in
 * Finite Element-type computations. Specialized dense vectors,
 * for example DenseSubVectors, can be derived from this class.
 *
 * \author John W. Peterson
 * \date 2003
 */
template<typename T>
class DenseVectorBase
{
public:
  /**
   * Constructor.
   */
  DenseVectorBase() = default;

  /**
   * The 5 special functions can be defaulted for this class, as it
   * does not manage any memory itself.
   */
  DenseVectorBase (DenseVectorBase &&) = default;
  DenseVectorBase (const DenseVectorBase &) = default;
  DenseVectorBase & operator= (const DenseVectorBase &) = default;
  DenseVectorBase & operator= (DenseVectorBase &&) = default;
  virtual ~DenseVectorBase() = default;

  /**
   * Set every element in the vector to 0.  Needs to
   * be pure virtual since the storage method may
   * be different in derived classes.
   */
  virtual void zero() = 0;

  /**
   * \returns The \p (i) element of the vector.
   */
  virtual T el(const unsigned int i) const = 0;

  /**
   * \returns The \p (i) element of the vector as a writable reference.
   */
  virtual T & el(const unsigned int i) = 0;

  /**
   * \returns The size of the vector.
   */
  virtual unsigned int size() const = 0;

  /**
   * \returns \p true iff size() is 0.
   */
  virtual bool empty() const { return (this->size() == 0); }

  /**
   * Pretty-print the vector to \p stdout.
   */
  void print(std::ostream & os) const;

  /**
   * Same as above, but allows you to print using the
   * usual stream syntax.
   */
  friend std::ostream & operator << (std::ostream & os, const DenseVectorBase<T> & v)
  {
    v.print(os);
    return os;
  }

  /**
   * Prints the entries of the vector with additional
   * decimal places in scientific notation.
   */
  void print_scientific(std::ostream & os, unsigned precision=8) const;
};

template<typename T>
void DenseVectorBase<T>::print (std::ostream & os) const
{
  for (auto i : make_range(this->size()))
    os << std::setw(8)
       << this->el(i)
       << std::endl;
}

} // namespace libMesh



#endif // LIBMESH_DENSE_VECTOR_BASE_H
