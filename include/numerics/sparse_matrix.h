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



#ifndef LIBMESH_SPARSE_MATRIX_H
#define LIBMESH_SPARSE_MATRIX_H


// Local includes
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"
#include "libmesh/enum_parallel_type.h" // PARALLEL
#include "libmesh/enum_matrix_build_type.h" // AUTOMATIC
#include "libmesh/enum_solver_package.h"

// C++ includes
#include <cstddef>
#include <iomanip>
#include <vector>
#include <memory>

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class DenseMatrix;
class DofMap;
namespace SparsityPattern {
  class Build;
  class Graph;
}
template <typename T> class NumericVector;

// This template helper function must be declared before it
// can be defined below.
template <typename T>
std::ostream & operator << (std::ostream & os, const SparseMatrix<T> & m);


/**
 * Generic sparse matrix. This class contains pure virtual members
 * that must be overridden in derived classes.  Using a common base
 * class allows for uniform access to sparse matrices from various
 * different solver packages in different formats.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
template <typename T>
class SparseMatrix : public ReferenceCountedObject<SparseMatrix<T>>,
                     public ParallelObject
{
public:
  /**
   * Constructor; initializes the matrix to be empty, without any
   * structure, i.e.  the matrix is not usable at all. This
   * constructor is therefore only useful for matrices which are
   * members of a class. All other matrices should be created at a
   * point in the data flow where all necessary information is
   * available.
   *
   * You have to initialize the matrix before usage with
   * \p init(...).
   */
  explicit
  SparseMatrix (const Parallel::Communicator & comm);

  /**
   * This _looks_ like a copy assignment operator, but note that,
   * unlike normal copy assignment operators, it is pure virtual. This
   * function should be overridden in derived classes so that they can
   * be copied correctly via references to the base class. This design
   * usually isn't a good idea in general, but in this context it
   * works because we usually don't have a mix of different kinds of
   * SparseMatrix active in the library at a single time.
   *
   * \returns A reference to *this as the base type.
   */
  virtual SparseMatrix<T> & operator= (const SparseMatrix<T> &)
  {
    libmesh_not_implemented();
    return *this;
  }

  /**
   * These 3 special functions can be defaulted for this class, as it
   * does not manage any memory itself.
   */
  SparseMatrix (SparseMatrix &&) = default;
  SparseMatrix (const SparseMatrix &) = default;
  SparseMatrix & operator= (SparseMatrix &&) = default;

  /**
   * While this class doesn't manage any memory, the derived class might and
   * users may be deleting through a pointer to this base class
   */
  virtual ~SparseMatrix() = default;

  /**
   * Builds a \p SparseMatrix<T> using the linear solver package specified by
   * \p solver_package
   */
  static std::unique_ptr<SparseMatrix<T>>
  build(const Parallel::Communicator & comm,
        const SolverPackage solver_package = libMesh::default_solver_package(),
        const MatrixBuildType matrix_build_type = MatrixBuildType::AUTOMATIC);

  virtual SolverPackage solver_package() = 0;

  /**
   * \returns \p true if the matrix has been initialized,
   * \p false otherwise.
   */
  virtual bool initialized() const { return _is_initialized; }

  /**
   * Set a pointer to the \p DofMap to use.  If a separate sparsity
   * pattern is not being used, use the one from the DofMap.
   *
   * The lifetime of \p dof_map must exceed the lifetime of \p this.
   */
  void attach_dof_map (const DofMap & dof_map);

  /**
   * Set a pointer to a sparsity pattern to use.  Useful in cases
   * where a matrix requires a wider (or for efficiency narrower)
   * pattern than most matrices in the system, or in cases where no
   * system sparsity pattern is being calculated by the DofMap.
   *
   * The lifetime of \p sp must exceed the lifetime of \p this.
   */
  void attach_sparsity_pattern (const SparsityPattern::Build & sp);

  /**
   * \returns \p true if this sparse matrix format needs to be fed the
   * graph of the sparse matrix.
   *
   * This is true for \p LaspackMatrix, but not \p PetscMatrixBase subclasses.  In
   * the case where the full graph is not required, we can efficiently
   * approximate it to provide a good estimate of the required size of
   * the sparse matrix.
   */
  virtual bool need_full_sparsity_pattern() const
  { return false; }

  /**
   * @returns Whether this matrix needs the sparsity pattern computed by the \p DofMap
   */
  virtual bool require_sparsity_pattern() const { return !this->use_hash_table(); }

  /**
   * Updates the matrix sparsity pattern. When your \p SparseMatrix<T>
   * implementation does not need this data, simply do not override
   * this method.
   */
  virtual void update_sparsity_pattern (const SparsityPattern::Graph &) {}

  /**
   * Initialize SparseMatrix with the specified sizes.
   *
   * \param m The global number of rows.
   * \param n The global number of columns.
   * \param m_l The local number of rows.
   * \param n_l The local number of columns.
   * \param nnz The number of on-diagonal nonzeros per row (defaults to 30).
   * \param noz The number of off-diagonal nonzeros per row (defaults to 10).
   * \param blocksize Optional value indicating dense coupled blocks
   * for systems with multiple variables all of the same type.
   */
  virtual void init (const numeric_index_type m,
                     const numeric_index_type n,
                     const numeric_index_type m_l,
                     const numeric_index_type n_l,
                     const numeric_index_type nnz=30,
                     const numeric_index_type noz=10,
                     const numeric_index_type blocksize=1) = 0;

  /**
   * Initialize this matrix using the sparsity structure computed by \p dof_map.
   * @param type The serial/parallel/ghosted type of the matrix
   */
  virtual void init (ParallelType type = PARALLEL) = 0;

  /**
   * Restores the \p SparseMatrix<T> to a pristine state.
   */
  virtual void clear () = 0;

  /**
   * Set all entries to 0.
   */
  virtual void zero () = 0;

  /**
   * \returns A smart pointer to a copy of this matrix with the same
   * type, size, and partitioning, but with all zero entries.
   *
   * \note This must be overridden in the derived classes.
   */
  virtual std::unique_ptr<SparseMatrix<T>> zero_clone () const = 0;

  /**
   * \returns A smart pointer to a copy of this matrix.
   *
   * \note This must be overridden in the derived classes.
   */
  virtual std::unique_ptr<SparseMatrix<T>> clone () const = 0;

  /**
   * Sets all row entries to 0 then puts \p diag_value in the diagonal entry.
   */
  virtual void zero_rows (std::vector<numeric_index_type> & rows, T diag_value = 0.0);

  /**
   * Calls the SparseMatrix's internal assembly routines, ensuring
   * that the values are consistent across processors.
   */
  virtual void close () = 0;

  /**
   *  For PETSc matrix , this function is similar to close but without shrinking memory.
   *  This is useful when we want to switch between ADD_VALUES and INSERT_VALUES.
   *  close should be called before using the matrix.
   */
  virtual void flush () { close(); }

  /**
   * \returns The row-dimension of the matrix.
   */
  virtual numeric_index_type m () const = 0;

  /**
   * Get the number of rows owned by this process
   */
  virtual numeric_index_type local_m () const { return row_stop() - row_start(); }

  /**
   * Get the number of columns owned by this process
   */
  virtual numeric_index_type local_n () const { return col_stop() - col_start(); }

  /**
   * \returns The column-dimension of the matrix.
   */
  virtual numeric_index_type n () const = 0;

  /**
   * \returns The index of the first matrix row stored on this
   * processor.
   */
  virtual numeric_index_type row_start () const = 0;

  /**
   * \returns The index of the last matrix row (+1) stored on this
   * processor.
   */
  virtual numeric_index_type row_stop () const = 0;

  /**
   * \returns The index of the first matrix column owned by this
   * processor.
   */
  virtual numeric_index_type col_start () const = 0;

  /**
   * \returns The index of the last matrix column (+1) owned by this
   * processor.
   */
  virtual numeric_index_type col_stop () const = 0;

  /**
   * Set the element \p (i,j) to \p value.  Throws an error if the
   * entry does not exist. Zero values can be "stored" in non-existent
   * fields.
   */
  virtual void set (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) = 0;

  /**
   * Add \p value to the element \p (i,j).  Throws an error if the
   * entry does not exist. Zero values can be "added" to non-existent
   * entries.
   */
  virtual void add (const numeric_index_type i,
                    const numeric_index_type j,
                    const T value) = 0;

  /**
   * Add the full matrix \p dm to the SparseMatrix.  This is useful
   * for adding an element matrix at assembly time.
   */
  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & rows,
                           const std::vector<numeric_index_type> & cols) = 0;

  /**
   * Same as \p add_matrix, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  virtual void add_matrix (const DenseMatrix<T> & dm,
                           const std::vector<numeric_index_type> & dof_indices) = 0;

  /**
   * Add the full matrix \p dm to the SparseMatrix.  This is useful
   * for adding an element matrix at assembly time.  The matrix is
   * assumed blocked, and \p brow, \p bcol correspond to the *block*
   * row and column indices.
   */
  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & brows,
                                 const std::vector<numeric_index_type> & bcols);

  /**
   * Same as \p add_block_matrix(), but assumes the row and column
   * maps are the same.  Thus the matrix \p dm must be square.
   */
  virtual void add_block_matrix (const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & dof_indices)
  { this->add_block_matrix (dm, dof_indices, dof_indices); }

  /**
   * Compute \f$ A \leftarrow A + a*X \f$ for scalar \p a, matrix \p X.
   */
  virtual void add (const T a, const SparseMatrix<T> & X) = 0;

  /**
   * Compute Y = A*X for matrix \p X.
   */
  virtual void matrix_matrix_mult (SparseMatrix<T> & /*X*/, SparseMatrix<T> & /*Y*/, bool /*reuse*/)
  { libmesh_not_implemented(); }

  /**
   * Add \p scalar* \p spm to the rows and cols of this matrix (A):
   * A(rows[i], cols[j]) += scalar * spm(i,j)
   */
  virtual void add_sparse_matrix (const SparseMatrix<T> & /*spm*/,
                                  const std::map<numeric_index_type, numeric_index_type> & /*rows*/,
                                  const std::map<numeric_index_type, numeric_index_type> & /*cols*/,
                                  const T /*scalar*/)
  { libmesh_not_implemented(); }

  /**
   * \returns A copy of matrix entry \p (i,j).
   *
   * \note This may be an expensive operation, and you should always
   * be careful where you call this function.
   */
  virtual T operator () (const numeric_index_type i,
                         const numeric_index_type j) const = 0;

  /**
   * \returns The \f$ \ell_1 \f$-norm of the matrix, that is the max column sum:
   * \f$ |M|_1 = \max_{j} \sum_{i} |M_{ij}| \f$
   *
   * This is the natural matrix norm that is compatible with the
   * \f$ \ell_1 \f$-norm for vectors, i.e. \f$ |Mv|_1 \leq |M|_1 |v|_1 \f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  virtual Real l1_norm () const = 0;

  /**
   * \returns The \f$ \ell_{\infty} \f$-norm of the matrix, that is the max row sum:
   *
   * \f$ |M|_{\infty} = \max_{i} \sum_{j} |M_{ij}| \f$
   *
   * This is the natural matrix norm that is compatible to the
   * \f$ \ell_{\infty} \f$-norm of vectors, i.e.
   * \f$ |Mv|_{\infty} \leq |M|_{\infty} |v|_{\infty} \f$.
   * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
   */
  virtual Real linfty_norm () const = 0;

  /**
   * \returns The l1_norm() of the difference of \p this and \p other_mat
   */
  Real l1_norm_diff (const SparseMatrix<T> & other_mat) const;

  /**
   * \returns \p true if the matrix has been assembled.
   */
  virtual bool closed() const = 0;

  /**
   * \returns the global number of non-zero entries in the matrix
   * sparsity pattern
   */
  virtual std::size_t n_nonzeros() const;

  /**
   * Print the contents of the matrix to the screen in a uniform
   * style, regardless of matrix/solver package being used.
   */
  void print(std::ostream & os=libMesh::out, const bool sparse=false) const;

  /**
   * Same as the print method above, but allows you to print to a
   * stream in the standard syntax.
   *
   * \code
   * template <typename U>
   * friend std::ostream & operator << (std::ostream & os, const SparseMatrix<U> & m);
   * \endcode
   *
   * \note The above syntax, which does not require any
   * prior declaration of operator<<, declares *any* instantiation of
   * SparseMatrix<X> is friend to *any* instantiation of
   * operator<<(ostream &, SparseMatrix<Y> &).  It would not happen in
   * practice, but in principle it means that SparseMatrix<Complex>
   * would be friend to operator<<(ostream &, SparseMatrix<Real>).
   *
   * \note The form below, which requires a previous
   * declaration of the operator<<(stream &, SparseMatrix<T> &) function
   * (see top of this file), means that any instantiation of
   * SparseMatrix<T> is friend to the specialization
   * operator<<(ostream &, SparseMatrix<T> &), but e.g. SparseMatrix<U>
   * is *not* friend to the same function.  So this is slightly
   * different to the form above...
   *
   * This method seems to be the "preferred" technique, see
   * http://www.parashift.com/c++-faq-lite/template-friends.html
   */
  friend std::ostream & operator << <>(std::ostream & os, const SparseMatrix<T> & m);

  /**
   * Print the contents of the matrix to the screen in a
   * package-personalized style, if available.
   */
  virtual void print_personal(std::ostream & os=libMesh::out) const = 0;

  /**
   * Print the contents of the matrix in Matlab's sparse matrix
   * format. Optionally prints the matrix to the file named \p name.
   * If \p name is not specified it is dumped to the screen.
   */
  virtual void print_matlab(const std::string & /*name*/ = "") const;

  /**
   * Write the contents of the matrix to a file in PETSc's binary
   * sparse matrix format.
   */
  virtual void print_petsc_binary(const std::string & filename);

  /**
   * Write the contents of the matrix to a file in PETSc's HDF5
   * sparse matrix format.
   */
  virtual void print_petsc_hdf5(const std::string & filename);


  /**
   * Read the contents of the matrix from a file, with the file format
   * inferred from the extension of \p filename.
   */
  virtual void read(const std::string & filename);

  /**
   * Read the contents of the matrix from the Matlab-script sparse
   * matrix format used by PETSc.
   *
   * If the size and sparsity of the matrix in \p filename appears
   * consistent with the existing sparsity of \p this then the
   * existing parallel decomposition and sparsity will be retained.
   * If not, then \p this will be initialized with the sparsity from
   * the file, linearly partitioned onto the number of processors
   * available.
   */
  virtual void read_matlab(const std::string & filename);

  /**
   * Read the contents of the matrix from a file in PETSc's binary
   * sparse matrix format.
   */
  virtual void read_petsc_binary(const std::string & filename);

  /**
   * Read the contents of the matrix from a file in PETSc's HDF5
   * sparse matrix format.
   */
  virtual void read_petsc_hdf5(const std::string & filename);

  /**
   * This function creates a matrix called "submatrix" which is defined
   * by the row and column indices given in the "rows" and "cols" entries.
   * Currently this operation is only defined for the PetscMatrixBase subclasses.
   * Note: The \p rows and \p cols vectors need to be sorted;
   *       Use the nosort version below if \p rows and \p cols vectors are not sorted;
   *       The \p rows and \p cols only contain indices that are owned by this processor.
   */
  virtual void create_submatrix(SparseMatrix<T> & submatrix,
                                const std::vector<numeric_index_type> & rows,
                                const std::vector<numeric_index_type> & cols) const
  {
    this->_get_submatrix(submatrix,
                         rows,
                         cols,
                         false); // false means DO NOT REUSE submatrix
  }

  /**
   * Similar to the above function, this function creates a \p
   * submatrix which is defined by the indices given in the \p rows
   * and \p cols vectors.
   * Note: Both \p rows and \p cols can be unsorted;
   *       Use the above function for better efficiency if your indices are sorted;
   *       \p rows and \p cols can contain indices that are owned by other processors.
   */
  virtual void create_submatrix_nosort(SparseMatrix<T> & /*submatrix*/,
                                         const std::vector<numeric_index_type> & /*rows*/,
                                         const std::vector<numeric_index_type> & /*cols*/) const
  {
    libmesh_not_implemented();
  }

  /**
   * This function is similar to the one above, but it allows you to reuse
   * the existing sparsity pattern of "submatrix" instead of reallocating
   * it again.  This should hopefully be more efficient if you are frequently
   * extracting submatrices of the same size.
   */
  virtual void reinit_submatrix(SparseMatrix<T> & submatrix,
                                const std::vector<numeric_index_type> & rows,
                                const std::vector<numeric_index_type> & cols) const
  {
    this->_get_submatrix(submatrix,
                         rows,
                         cols,
                         true); // true means REUSE submatrix
  }

  /**
   * Multiplies the matrix by the NumericVector \p arg and stores the
   * result in NumericVector \p dest.
   */
  void vector_mult (NumericVector<T> & dest,
                    const NumericVector<T> & arg) const;

  /**
   * Multiplies the matrix by the NumericVector \p arg and adds the
   * result to the NumericVector \p dest.
   */
  void vector_mult_add (NumericVector<T> & dest,
                        const NumericVector<T> & arg) const;

  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVector<T> & dest) const = 0;

  /**
   * Copies the transpose of the matrix into \p dest, which may be
   * *this.
   */
  virtual void get_transpose (SparseMatrix<T> & dest) const = 0;

  /**
   * Get a row from the matrix
   * @param i The matrix row to get
   * @param indices A container that will be filled with the column indices
   *                corresponding to (possibly) non-zero values
   * @param values A container holding the column values
   */
  virtual void get_row(numeric_index_type i,
                       std::vector<numeric_index_type> & indices,
                       std::vector<T> & values) const = 0;

  /**
   * Scales all elements of this matrix by \p scale
   */
  virtual void scale(const T scale);

  /**
   * @returns Whether the matrix supports hash table assembly
   */
  virtual bool supports_hash_table() const { return false; }

  /**
   * Sets whether to use hash table assembly. This will error if the passed-in value is true and the
   * matrix type does not support hash tables. Hash table or hash map assembly means storing maps
   * from i-j locations in the matrix to values. Because it is a hash map as opposed to a contiguous
   * array of data, no preallocation is required to use it
   */
  void use_hash_table(bool use_hash);

  /**
   * @returns Whether this matrix is using hash table assembly. Hash table or hash map assembly
   * means storing maps from i-j locations in the matrix to values. Because it is a hash map as
   * opposed to a contiguous array of data, no preallocation is required to use it
   */
  bool use_hash_table() const { return _use_hash_table; }

  /**
   * Reset the memory storage of the matrix. Unlike \p clear(), this does not destroy the matrix but
   * rather will reset the matrix to use the original preallocation or when using hash table matrix
   * assembly (see \p use_hash_table()) will reset (clear) the hash table used for assembly. In the
   * words of the \p MatResetPreallocation documentation in PETSc, 'current values in the matrix are
   * lost in this call', so a user can expect to have back their original sparsity pattern in a
   * zeroed state
   */
  virtual void restore_original_nonzero_pattern() { libmesh_not_implemented(); }

protected:
  /**
   * Protected implementation of the create_submatrix and reinit_submatrix
   * routines.
   *
   * \note This function must be overridden in derived classes for it
   * to work properly!
   */
  virtual void _get_submatrix(SparseMatrix<T> & /*submatrix*/,
                              const std::vector<numeric_index_type> & /*rows*/,
                              const std::vector<numeric_index_type> & /*cols*/,
                              const bool /*reuse_submatrix*/) const
  {
    libmesh_not_implemented();
  }

  /**
   * The \p DofMap object associated with this object.  May be queried
   * for degree-of-freedom counts on processors.
   */
  DofMap const * _dof_map;

  /**
   * The \p sparsity pattern associated with this object.  Should be
   * queried for entry counts (or with need_full_sparsity_pattern,
   * patterns) when needed.
   */
  SparsityPattern::Build const * _sp;

  /**
   * Flag indicating whether or not the matrix has been initialized.
   */
  bool _is_initialized;

  /**
   * Flag indicating whether the matrix is assembled using a hash table
   */
  bool _use_hash_table;
};



//-----------------------------------------------------------------------
// SparseMatrix inline members
template <typename T>
void
SparseMatrix<T>::use_hash_table(const bool use_hash)
{
  libmesh_error_msg_if(use_hash && !this->supports_hash_table(),
                       "This matrix class does not support hash table assembly");
  this->_use_hash_table = use_hash;
}

// For SGI MIPSpro this implementation must occur after
// the full specialization of the print() member.
//
// It's generally easier to define these friend functions in the header
// file.
template <typename T>
std::ostream & operator << (std::ostream & os, const SparseMatrix<T> & m)
{
  m.print(os);
  return os;
}

template <typename T>
auto
l1_norm(const SparseMatrix<T> & mat)
{
  return mat.l1_norm();
}

template <typename T>
auto
l1_norm_diff(const SparseMatrix<T> & mat1, const SparseMatrix<T> & mat2)
{
  return mat1.l1_norm_diff(mat2);
}

} // namespace libMesh


#endif // LIBMESH_SPARSE_MATRIX_H
