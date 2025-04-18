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


// Local Includes
#include "libmesh/numeric_vector.h"
#include "libmesh/distributed_vector.h"
#include "libmesh/laspack_vector.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/trilinos_epetra_vector.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/int_range.h"


// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::abs
#include <limits>
#include <memory>


namespace libMesh
{



//------------------------------------------------------------------
// NumericVector methods

// Full specialization for Real datatypes
template <typename T>
std::unique_ptr<NumericVector<T>>
NumericVector<T>::build(const Parallel::Communicator & comm,
                        SolverPackage solver_package,
                        ParallelType parallel_type)
{
  // Build the appropriate vector
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      return std::make_unique<LaspackVector<T>>(comm, parallel_type);
#endif

#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      return std::make_unique<PetscVector<T>>(comm, parallel_type);
#endif

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA
    case TRILINOS_SOLVERS:
      return std::make_unique<EpetraVector<T>>(comm, parallel_type);
#endif

#ifdef LIBMESH_HAVE_EIGEN
    case EIGEN_SOLVERS:
      return std::make_unique<EigenSparseVector<T>>(comm, parallel_type);
#endif

    default:
      return std::make_unique<DistributedVector<T>>(comm, parallel_type);
    }
}



template <typename T>
void NumericVector<T>::set_type(ParallelType t)
{
  // Check for no-op
  if (_type == t)
    return;

  // If the NumericVector is not yet initialized, then it is generally
  // safe to change the ParallelType, with minor restrictions.
  if (!this->initialized())
    {
      // If ghosted vectors are not enabled and the user requested a
      // GHOSTED vector, fall back on SERIAL.
#ifndef LIBMESH_ENABLE_GHOSTED
      if (t == GHOSTED)
        {
          _type = SERIAL;
          return;
        }
#endif

      _type = t;
      return;
    }

  // If we made it here, then the NumericVector was already
  // initialized and we don't currently allow the ParallelType to be
  // changed, although this could potentially be added later.
  libmesh_not_implemented();
}

template <typename T>
void NumericVector<T>::insert (const T * v,
                               const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert (v);

  for (auto i : index_range(dof_indices))
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void NumericVector<T>::insert (const NumericVector<T> & V,
                               const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());
  libmesh_assert (V.readable());

  for (auto i : index_range(dof_indices))
    this->set (dof_indices[i], V(i));
}



template <typename T>
int NumericVector<T>::compare (const NumericVector<T> & other_vector,
                               const Real threshold) const
{
  libmesh_assert(this->compatible(other_vector));

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index())
  {
    if (std::abs((*this)(i) - other_vector(i)) > threshold)
      first_different_i = i;
    else
      i++;
  }

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}


template <typename T>
int NumericVector<T>::local_relative_compare (const NumericVector<T> & other_vector,
                                              const Real threshold) const
{
  libmesh_assert(this->compatible(other_vector));

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  do
    {
      if (std::abs((*this)(i) - other_vector(i)) > threshold *
          std::max(std::abs((*this)(i)), std::abs(other_vector(i))))
        first_different_i = i;
      else
        i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index());

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}


template <typename T>
int NumericVector<T>::global_relative_compare (const NumericVector<T> & other_vector,
                                               const Real threshold) const
{
  libmesh_assert(this->compatible(other_vector));

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  const Real my_norm = this->linfty_norm();
  const Real other_norm = other_vector.linfty_norm();
  const Real abs_threshold = std::max(my_norm, other_norm) * threshold;

  do
    {
      if (std::abs((*this)(i) - other_vector(i) ) > abs_threshold)
        first_different_i = i;
      else
        i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index());

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}

/*
// Full specialization for float datatypes (DistributedVector wants this)

template <>
int NumericVector<float>::compare (const NumericVector<float> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if (std::abs((*this)(i) - other_vector(i) ) > threshold)
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<last_local_index());

return rvalue;
}

// Full specialization for double datatypes
template <>
int NumericVector<double>::compare (const NumericVector<double> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if (std::abs((*this)(i) - other_vector(i) ) > threshold)
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<last_local_index());

return rvalue;
}

#ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
// Full specialization for long double datatypes
template <>
int NumericVector<long double>::compare (const NumericVector<long double> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if (std::abs((*this)(i) - other_vector(i) ) > threshold)
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<last_local_index());

return rvalue;
}
#endif


// Full specialization for Complex datatypes
template <>
int NumericVector<Complex>::compare (const NumericVector<Complex> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if ((std::abs((*this)(i).real() - other_vector(i).real()) > threshold) || (std::abs((*this)(i).imag() - other_vector(i).imag()) > threshold))
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<this->last_local_index());

return rvalue;
}
*/


template <class T>
Real NumericVector<T>::subset_l1_norm (const std::set<numeric_index_type> & indices) const
{
  libmesh_assert (this->readable());

  const NumericVector<T> & v = *this;

  Real norm = 0;

  for (const auto & index : indices)
    norm += std::abs(v(index));

  this->comm().sum(norm);

  return norm;
}

template <class T>
Real NumericVector<T>::subset_l2_norm (const std::set<numeric_index_type> & indices) const
{
  libmesh_assert (this->readable());

  const NumericVector<T> & v = *this;

  Real norm = 0;

  for (const auto & index : indices)
    norm += TensorTools::norm_sq(v(index));

  this->comm().sum(norm);

  return std::sqrt(norm);
}

template <class T>
Real NumericVector<T>::subset_linfty_norm (const std::set<numeric_index_type> & indices) const
{
  libmesh_assert (this->readable());

  const NumericVector<T> & v = *this;

  Real norm = 0;

  for (const auto & index : indices)
    {
      Real value = std::abs(v(index));
      if (value > norm)
        norm = value;
    }

  this->comm().max(norm);

  return norm;
}



template <class T>
Real NumericVector<T>::l2_norm_diff (const NumericVector<T> & v) const
{
  libmesh_assert(this->compatible(v));

  Real norm = 0;
  for (const auto i : make_range(this->first_local_index(), this->last_local_index()))
    norm += TensorTools::norm_sq((*this)(i) - v(i));

  this->comm().sum(norm);

  return std::sqrt(norm);
}



template <class T>
Real NumericVector<T>::l1_norm_diff (const NumericVector<T> & v) const
{
  libmesh_assert(this->compatible(v));

  Real norm = 0;
  for (const auto i : make_range(this->first_local_index(), this->last_local_index()))
    norm += libMesh::l1_norm_diff((*this)(i), v(i));

  this->comm().sum(norm);

  return norm;
}



template <typename T>
void NumericVector<T>::add_vector (const T * v,
                                   const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v);

  for (auto i : index_range(dof_indices))
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void NumericVector<T>::add_vector (const NumericVector<T> & v,
                                   const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v.readable());

  const std::size_t n = dof_indices.size();
  libmesh_assert_equal_to(v.size(), n);
  for (numeric_index_type i=0; i != n; i++)
    this->add (dof_indices[i], v(i));
}



template <typename T>
void NumericVector<T>::add_vector (const NumericVector<T> & v,
                                   const ShellMatrix<T> & a)
{
  libmesh_assert(this->compatible(v));

  a.vector_mult_add(*this,v);
}



template <typename T>
bool NumericVector<T>::readable () const
{
  return this->initialized() && this->closed();
}


template <typename T>
bool NumericVector<T>::compatible (const NumericVector<T> & v) const
{
  return this->readable() && v.readable() &&
         this->size() == v.size() &&
         this->local_size() == v.local_size() &&
         this->first_local_index() == v.first_local_index() &&
         this->last_local_index() == v.last_local_index();
}


//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT NumericVector<Number>;

} // namespace libMesh
