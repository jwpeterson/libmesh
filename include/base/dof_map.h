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



#ifndef LIBMESH_DOF_MAP_H
#define LIBMESH_DOF_MAP_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/variable.h"
#include "libmesh/threads.h"
#include "libmesh/threads_allocators.h"
#include "libmesh/elem_range.h"
#include "libmesh/ghosting_functor.h"
#include "libmesh/sparsity_pattern.h"
#include "libmesh/parallel_object.h"
#include "libmesh/point.h"
#include "libmesh/utility.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/mesh_subdivision_support.h"
#include "libmesh/dof_map_base.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"
#include "timpi/parallel_sync.h"

// C++ Includes
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <map>
#include <string>
#include <vector>
#include <memory>

namespace libMesh
{

// Forward Declarations
class CouplingMatrix;
class DefaultCoupling;
class DirichletBoundary;
class DirichletBoundaries;
class DofMap;
class DofObject;
class FEType;
class MeshBase;
class PeriodicBoundaryBase;
class PeriodicBoundaries;
class System;
class NonlinearImplicitSystem;
class StaticCondensationDofMap;
template <typename T> class DenseVectorBase;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
enum Order : int;



// ------------------------------------------------------------
// Do we need constraints for anything?

#if defined(LIBMESH_ENABLE_AMR) ||              \
  defined(LIBMESH_ENABLE_PERIODIC) ||           \
  defined(LIBMESH_ENABLE_DIRICHLET)
#  define LIBMESH_ENABLE_CONSTRAINTS 1
#endif

// ------------------------------------------------------------
// AMR constraint matrix types

#ifdef LIBMESH_ENABLE_CONSTRAINTS
/**
 * A row of the Dof constraint matrix.
 */
typedef std::map<dof_id_type, Real,
                 std::less<dof_id_type>,
                 Threads::scalable_allocator<std::pair<const dof_id_type, Real>>> DofConstraintRow;

/**
 * The constraint matrix storage format.
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Don't delete this from
 * a pointer-to-std-map; the destructor isn't virtual!
 */
class DofConstraints : public std::map<dof_id_type,
                                       DofConstraintRow,
                                       std::less<dof_id_type>,
                                       Threads::scalable_allocator<std::pair<const dof_id_type, DofConstraintRow>>>
{
};

/**
 * Storage for DofConstraint right hand sides for a particular
 * problem.  Each dof id with a non-zero constraint offset
 * stores it in such a structure.
 */
class DofConstraintValueMap :
    public std::map<dof_id_type, Number,
                    std::less<dof_id_type>,
                    Threads::scalable_allocator<std::pair<const dof_id_type, Number>>>
{
};

/**
 * Storage for DofConstraint right hand sides for all adjoint
 * problems.
 */
class AdjointDofConstraintValues :
    public std::map<unsigned int, DofConstraintValueMap,
                    std::less<unsigned int>,
                    Threads::scalable_allocator
                    <std::pair<const unsigned int, DofConstraintValueMap>>>
{
};

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
/**
 * A row of the Node constraint mapping.  Currently this just
 * stores the topology of the constrained Nodes, but for forward
 * compatibility we also include coefficients, so we could add
 * Lagrange-positioned-node constraints later.
 */
typedef std::map<const Node *, Real,
                 std::less<const Node *>,
                 Threads::scalable_allocator<std::pair<const Node * const, Real>>> NodeConstraintRow;

/**
 * The Node constraint storage format.
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Don't delete this from
 * a pointer-to-std-map; the destructor isn't virtual!
 */
class NodeConstraints : public std::map<const Node *,
                                        std::pair<NodeConstraintRow,Point>,
                                        std::less<const Node *>,
                                        Threads::scalable_allocator<std::pair<const Node * const, std::pair<NodeConstraintRow,Point>>>>
{
};
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif // LIBMESH_ENABLE_CONSTRAINTS



/**
 * This class handles the numbering of degrees of freedom on a mesh.
 * For systems of equations the class supports a fixed number of variables.
 * The degrees of freedom are numbered such that sequential, contiguous blocks
 * belong to distinct processors.  This is so that the resulting data
 * structures will work well with parallel linear algebra packages.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 * \brief Manages the degrees of freedom (DOFs) in a simulation.
 */
class DofMap : public DofMapBase,
               public ReferenceCountedObject<DofMap>
{
public:

  /**
   * Constructor.  Requires the number of the system for which we
   * will be numbering degrees of freedom & the parent object
   * we are contained in, which defines our communication space.
   */
  explicit
  DofMap(const unsigned int sys_number,
         MeshBase & mesh);

  /**
   * Destructor.
   */
  ~DofMap();

  /**
   * Backwards compatibility for prior AugmentSparsityPattern users.
   */
  class AugmentSparsityPattern : public SparsityPattern::AugmentSparsityPattern
  {};

  /**
   * Abstract base class to be used to add user-defined parallel
   * degree of freedom couplings.
   */
  class AugmentSendList
  {
  public:
    virtual ~AugmentSendList () = default;

    /**
     * User-defined function to augment the send list.
     */
    virtual void augment_send_list (std::vector<dof_id_type> & send_list) = 0;
  };

  /**
   * Additional matrices may be attached to this \p DofMap.
   * They are initialized to the same sparsity structure as
   * the major matrix.
   */
  void attach_matrix (SparseMatrix<Number> & matrix);

  /**
   * Additional matrices may be be temporarily initialized by this \p
   * DofMap.
   * They are initialized to the same sparsity structure as
   * the major matrix.
   */
  void update_sparsity_pattern(SparseMatrix<Number> & matrix) const;

  /**
   * Matrices should not be attached more than once.  We can test for
   * an already-attached matrix if necessary using \p is_attached
   */
  bool is_attached (SparseMatrix<Number> & matrix);

  /**
   * Distribute dofs on the current mesh.  Also builds the send list for
   * processor \p proc_id, which defaults to 0 for ease of use in serial
   * applications.
   * \returns The total number of DOFs for the System, summed across all procs.
   */
  std::size_t distribute_dofs (MeshBase &);

  /**
   * Computes the sparsity pattern for the matrices corresponding to
   * \p proc_id and sends that data to Linear Algebra packages for
   * preallocation of sparse matrices.
   */
  void compute_sparsity (const MeshBase &);

  /**
   * Returns true iff a sparsity pattern has already been computed.
   */
  bool computed_sparsity_already () const;

  /**
   * Sets the current policy for constructing sparsity patterns: if
   * \p use_constraints is true (for robustness), we explicitly
   * account for sparsity entries created by constraint matrix pre-
   * and post- application.  If \p use_constraints is false (for
   * speed), we calculate only the sparsity pattern of an
   * unconstrained matrix.  This is false by default, because in
   * nearly all applications our constraints do not increase the
   * number of non-zeros required in a sparse matrix.
   */
  void set_constrained_sparsity_construction(bool use_constraints);

  /**
   * Sets need_full_sparsity_pattern to true regardless of the requirements by matrices
   */
  void full_sparsity_pattern_needed();

  /**
   * Returns true iff the current policy when constructing sparsity
   * patterns is to explicitly account for sparsity entries created by
   * constraint matrix pre- and post- application.
   */
  bool constrained_sparsity_construction();

  /**
   * Clears the sparsity pattern
   */
  void clear_sparsity();

  /**
   * Remove any default ghosting functor(s).  User-added ghosting
   * functors will be unaffected.
   *
   * Unless user-added equivalent ghosting functors exist, removing
   * the default coupling functor is only safe for explicit solves,
   * and removing the default algebraic ghosting functor is only safe
   * for codes where no evaluations on neighbor cells (e.g. no jump
   * error estimators) are done.
   *
   * Defaults can be restored manually via add_default_ghosting(), or
   * automatically if clear() returns the DofMap to a default state.
   */
  void remove_default_ghosting();

  /**
   * Add the default functor(s) for coupling and algebraic ghosting.
   * User-added ghosting functors will be unaffected.
   */
  void add_default_ghosting();

  /**
   * Iterator type for coupling and algebraic ghosting functor ranges.
   * This has changed in the past and may change again; code should
   * use auto or the type here.
   */
  typedef std::vector<GhostingFunctor *>::const_iterator GhostingFunctorIterator;

  /**
   * Adds a functor which can specify coupling requirements for
   * creation of sparse matrices.
   * Degree of freedom pairs which match the elements and variables
   * returned by these functors will be added to the sparsity pattern,
   * and the degrees of freedom which live on other processors will be
   * added to the send_list for use on ghosted vectors, and the
   * elements which live on other processors will be ghosted on a
   * distributed mesh.
   *
   * GhostingFunctor memory must be managed by the code which calls
   * this function; the GhostingFunctor lifetime is expected to extend
   * until either the functor is removed or the DofMap is destructed.
   *
   * When \p to_mesh is true, the \p coupling_functor is also added to
   * our associated mesh, to ensure that coupled elements do not get
   * lost during mesh distribution.  (if coupled elements were
   * *already* lost there's no getting them back after the fact,
   * sorry)
   *
   * If \p to_mesh is false, no change to mesh ghosting is made;
   * the Mesh must already have ghosting functor(s) specifying a
   * superset of \p coupling_functor or this is a horrible bug.
   */
  void add_coupling_functor(GhostingFunctor & coupling_functor,
                            bool to_mesh = true);

  /**
   * Adds a functor which can specify coupling requirements for
   * creation of sparse matrices.
   *
   * GhostingFunctor memory when using this method is managed by the
   * shared_ptr mechanism.
   */
  void add_coupling_functor(std::shared_ptr<GhostingFunctor> coupling_functor,
                            bool to_mesh = true)
  { _shared_functors[coupling_functor.get()] = coupling_functor;
    this->add_coupling_functor(*coupling_functor, to_mesh); }

  /**
   * Removes a functor which was previously added to the set of
   * coupling functors, from both this DofMap and from the underlying
   * mesh.
   */
  void remove_coupling_functor(GhostingFunctor & coupling_functor);

  /**
   * Beginning of range of coupling functors
   */
  GhostingFunctorIterator coupling_functors_begin() const
  { return _coupling_functors.begin(); }

  /**
   * End of range of coupling functors
   */
  GhostingFunctorIterator coupling_functors_end() const
  { return _coupling_functors.end(); }

  /**
   * Default coupling functor
   */
  DefaultCoupling & default_coupling() { return *_default_coupling; }

  /**
   * Adds a functor which can specify algebraic ghosting requirements
   * for use with distributed vectors.  Degrees of freedom on other
   * processors which match the elements and variables returned by
   * these functors will be added to the send_list, and the elements
   * on other processors will be ghosted on a distributed mesh, so
   * that the elements can always be found and the solutions on them
   * will always be evaluable.
   *
   * GhostingFunctor memory must be managed by the code which calls
   * this function; the GhostingFunctor lifetime is expected to extend
   * until either the functor is removed or the DofMap is destructed.
   *
   * When \p to_mesh is true, the \p coupling_functor is also added to
   * our associated mesh, to ensure that evaluable elements do not get
   * lost during mesh distribution.  (if evaluable elements were
   * *already* lost there's no getting them back after the fact,
   * sorry)
   *
   * If \p to_mesh is false, no change to mesh ghosting is made;
   * the Mesh must already have ghosting functor(s) specifying a
   * superset of \p evaluable_functor or this is a horrible bug.
   */
  void add_algebraic_ghosting_functor(GhostingFunctor & evaluable_functor,
                                      bool to_mesh = true);

  /**
   * Adds a functor which can specify algebraic ghosting requirements
   * for use with distributed vectors.
   *
   * GhostingFunctor memory when using this method is managed by the
   * shared_ptr mechanism.
   */
  void add_algebraic_ghosting_functor(std::shared_ptr<GhostingFunctor> evaluable_functor,
                                      bool to_mesh = true)
  { _shared_functors[evaluable_functor.get()] = evaluable_functor;
    this->add_algebraic_ghosting_functor(*evaluable_functor, to_mesh); }

  /**
   * Removes a functor which was previously added to the set of
   * algebraic ghosting functors, from both this DofMap and from the
   * underlying mesh.
   */
  void remove_algebraic_ghosting_functor(GhostingFunctor & evaluable_functor);

  /**
   * Beginning of range of algebraic ghosting functors
   */
  GhostingFunctorIterator algebraic_ghosting_functors_begin() const
  { return _algebraic_ghosting_functors.begin(); }

  /**
   * End of range of algebraic ghosting functors
   */
  GhostingFunctorIterator algebraic_ghosting_functors_end() const
  { return _algebraic_ghosting_functors.end(); }

  /**
   * Default algebraic ghosting functor
   */
  DefaultCoupling & default_algebraic_ghosting() { return *_default_evaluating; }

  /**
   * Attach an object to use to populate the
   * sparsity pattern with extra entries.
   *
   * Care must be taken that when adding entries they are sorted into the Rows
   *
   * Further, you _must_ modify n_nz and n_oz properly!
   *
   * This is an advanced function... use at your own peril!
   */
  void attach_extra_sparsity_object (SparsityPattern::AugmentSparsityPattern & asp)
  {
    _augment_sparsity_pattern = &asp;
  }

  /**
   * Attach a function pointer to use as a callback to populate the
   * sparsity pattern with extra entries.
   *
   * Care must be taken that when adding entries they are sorted into the Rows
   *
   * Further, you _must_ modify n_nz and n_oz properly!
   *
   * This is an advanced function... use at your own peril!
   */
  void attach_extra_sparsity_function(void (*func)(SparsityPattern::Graph & sparsity,
                                                   std::vector<dof_id_type> & n_nz,
                                                   std::vector<dof_id_type> & n_oz,
                                                   void *),
                                      void * context = nullptr)
  { _extra_sparsity_function = func; _extra_sparsity_context = context; }

  /**
   * Attach an object to populate the send_list with extra entries.
   * This should only add to the send list, but no checking is done
   * to enforce this behavior.
   *
   * This is an advanced function... use at your own peril!
   */
  void attach_extra_send_list_object (DofMap::AugmentSendList & asl)
  {
    _augment_send_list = &asl;
  }

  /**
   * Attach a function pointer to use as a callback to populate the
   * send_list with extra entries.
   */
  void attach_extra_send_list_function(void (*func)(std::vector<dof_id_type> &, void *),
                                       void * context = nullptr)
  { _extra_send_list_function = func; _extra_send_list_context = context; }

  /**
   * Takes the \p _send_list vector (which may have duplicate entries)
   * and sorts it.  The duplicate entries are then removed, resulting in
   * a sorted \p _send_list with unique entries.  Also calls any user-provided
   * methods for adding to the send list.
   */
  void prepare_send_list ();

  /**
   * Clears the \p _send_list vector. This should be done in order to completely
   * rebuild the send_list from scratch rather than merely adding to the existing
   * send_list.
   */
  void clear_send_list ()
  {
    _send_list.clear();
  }

  /**
   * Clears the \p _send_list vector and then rebuilds it. This may be needed
   * in special situations, for example when an algebraic coupling functor cannot
   * be added to the \p DofMap until after it is completely setup. Then this method
   * can be used to rebuild the send_list once the algebraic coupling functor is
   * added. Note that while this will recommunicate constraints with the updated
   * send_list, this does assume no new constraints have been added since the previous
   * reinit_constraints call.
   */
  void reinit_send_list (MeshBase & mesh);


  /**
   * \returns A constant reference to the \p _send_list for this processor.
   *
   * The \p _send_list contains the global indices of all the
   * variables in the global solution vector that influence the
   * current processor.  This information can be used for gathers at
   * each solution step to retrieve solution values needed for
   * computation.
   */
  const std::vector<dof_id_type> & get_send_list() const { return _send_list; }

  /**
   * \returns A constant reference to the \p _n_nz list for this processor.
   *
   * The vector contains the bandwidth of the on-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<dof_id_type> & get_n_nz() const
  {
    libmesh_assert(_sp);
    return _sp->get_n_nz();
  }

  /**
   * \returns A constant reference to the \p _n_oz list for this processor.
   *
   * The vector contains the bandwidth of the off-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<dof_id_type> & get_n_oz() const
  {
    libmesh_assert(_sp);
    return _sp->get_n_oz();
  }


  /**
   * \returns A constant pointer to the sparsity pattern stored here,
   * once that has been computed.  Returns null if no sparsity pattern
   * has yet been computed.
   *
   * If need_full_sparsity_pattern is false, the "sparsity pattern"
   * may only own n_nz and n_oz lists.
   */
  const SparsityPattern::Build * get_sparsity_pattern() const
  {
    return _sp.get();
  }

  // /**
  //  * Add an unknown of order \p order and finite element type
  //  * \p type to the system of equations.
  //  */
  // void add_variable (const Variable & var);

  /**
   * Add a group of unknowns of order \p order and finite element type
   * \p type to the system of equations.
   */
  void add_variable_group (VariableGroup var_group);

  /**
   * Specify whether or not we perform an extra (opt-mode enabled) check
   * for constraint loops. If a constraint loop is present then
   * the system constraints are not valid, so if \p error_on_constraint_loop
   * is true we will throw an error in this case.
   *
   * \note We previously referred to these types of constraints as
   * "cyclic" but that has now been deprecated, and these will now
   * instead be referred to as "constraint loops" in libMesh.
   */
  void set_error_on_cyclic_constraint(bool error_on_cyclic_constraint);
  void set_error_on_constraint_loop(bool error_on_constraint_loop);

  /**
   * \returns The \p VariableGroup description object for group \p g.
   */
  const VariableGroup & variable_group (const unsigned int c) const;

  const Variable & variable (const unsigned int c) const override;

  /**
   * \returns The approximation order for variable \p c.
   */
  Order variable_order (const unsigned int c) const;

  /**
   * \returns The approximation order for \p VariableGroup \p vg.
   */
  Order variable_group_order (const unsigned int vg) const;

  /**
   * \returns The finite element type for variable \p c.
   */
  const FEType & variable_type (const unsigned int c) const;

  /**
   * \returns The finite element type for \p VariableGroup \p vg.
   */
  const FEType & variable_group_type (const unsigned int vg) const;

  /**
   * \returns The number of variables in the global solution vector. Defaults
   * to 1, should be 1 for a scalar equation, 3 for 2D incompressible Navier
   * Stokes (u,v,p), etc...
   */
  unsigned int n_variable_groups() const
  { return cast_int<unsigned int>(_variable_groups.size()); }

  unsigned int n_variables() const override
  { return cast_int<unsigned int>(_variables.size()); }

  /**
   * \returns The variable group number that the provided variable number belongs to
   */
  unsigned int var_group_from_var_number(unsigned int var_num) const;

  /**
   * \returns \p true if the variables are capable of being stored in a blocked
   * form.  Presently, this means that there can only be one variable group,
   * and that the group has more than one variable.
   */
  bool has_blocked_representation() const
  {
    return ((this->n_variable_groups() == 1) && (this->n_variables() > 1));
  }

  /**
   * \returns The block size, if the variables are amenable to block storage.
   * Otherwise 1.
   * This routine was originally designed to enable a blocked storage, but
   * it turns out this information is still super useful for solvers even when
   * we do not use the blocked storage (e.g., MATMPIBAIJ in PETSc). For example (in PCHMG),
   * for a system of PDEs, to construct an efficient multilevel preconditioner, we coarsen
   * the matrix of one single PDE instead of the entire huge matrix. In order to
   * accomplish this, we need to know many PDEs we have. Another use case,
   * the fieldsplit preconditioner can be constructed in place with this info without
   * involving any user efforts.
   */
  unsigned int block_size() const
  {
    return (this->has_blocked_representation() ? this->n_variables() : 1);
  }

  using DofMapBase::n_dofs;
  /**
   * \returns The total number of degrees of freedom for a particular
   * variable \p vn.
   */
  dof_id_type n_dofs(const unsigned int vn) const
  {
    dof_id_type n = this->n_local_dofs(vn);
    this->comm().sum(n);
    return n;
  }

  /**
   * \returns The number of SCALAR dofs.
   */
  dof_id_type n_SCALAR_dofs() const { return _n_SCALAR_dofs; }

  using DofMapBase::n_local_dofs;
  /**
   * \returns The number of degrees of freedom on this processor for a
   * particular variable \p vn.  This is an O(N) operation on serialized or
   * O(N/Nproc) operation on distributed meshes.
   */
  dof_id_type n_local_dofs(const unsigned int vn) const
  {
    dof_id_type n;
    this->local_variable_indices(n, _mesh, vn);
    return n;
  }

  /**
   * \returns The number of degrees of freedom on each partition for a
   * particular variable \p vn.
   */
  std::vector<dof_id_type> n_dofs_per_processor(const unsigned int vn) const
  {
    std::vector<dof_id_type> n_local_dofs(this->n_processors(), 0);
    this->comm().allgather(this->n_local_dofs(vn), n_local_dofs);
    return n_local_dofs;
  }

  /**
   * \returns The processor id that owns the dof index \p dof
   */
  processor_id_type dof_owner(const dof_id_type dof) const
  { std::vector<dof_id_type>::const_iterator ub =
      std::upper_bound(_end_df.begin(), _end_df.end(), dof);
    libmesh_assert (ub != _end_df.end());
    return cast_int<processor_id_type>(ub - _end_df.begin());
  }

  void dof_indices (const Elem * const elem,
                    std::vector<dof_id_type> & di) const;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element.  For one variable, and potentially for a
   * non-default element p refinement level
   */
  void dof_indices (const Elem * const elem,
                    std::vector<dof_id_type> & di,
                    const unsigned int vn,
                    int p_level = -12345) const override;

  /**
   * Retrieves degree of freedom indices for a given \p elem and then performs actions for these
   * indices defined by the user-provided functors \p scalar_dofs_functor and \p field_dofs_functor.
   * This API is useful when a user wants to do more than simply fill a degree of freedom container
   * @param elem The element to get degrees of freedom for
   * @param di A container for degrees of freedom. It is up to the provided functors how this gets
   *           filled
   * @param vn The variable number to retrieve degrees of freedom for
   * @param scalar_dofs_functor The functor that acts on scalar degrees of freedom. This functor has
   *                            the interface:
   *                            void scalar_dofs_functor(const Elem & elem,
   *                                                    std::vector<dof_id_type> & di,
   *                                                    const std::vector<dof_id_type> & scalar_dof_indices)
   *                            where \p di is the degree of freedom container described above and
   *                            \p scalar_dof_indices are the scalar dof indices available to
   *                            \p elem
   * @param field_dofs_functor The functor that acts on "field" (e.g. non-scalar, non-global)
   *                           degrees of freedom. This functor has
   *                           the interface:
   *                           void field_dofs_functor(const Elem & elem,
   *                                                   const unsigned int node_num,
   *                                                   const unsigned int var_num,
   *                                                   std::vector<dof_id_type> & di,
   *                                                   const dof_id_type field_dof)
   *                           where \p field_dof represents a field degree of freedom to act on and
   *                           is associated with \p node_num and \p var_num. If the degree of
   *                           freedom is elemental than \p node_num will be \p invalid_uint. \p di
   *                           is again the degree of freedom container provided above
   */
  template <typename ScalarDofsFunctor, typename FieldDofsFunctor>
  void dof_indices(const Elem * const elem,
                   std::vector<dof_id_type> & di,
                   const unsigned int vn,
                   ScalarDofsFunctor scalar_dofs_functor,
                   FieldDofsFunctor field_dofs_functor,
                   int p_level = -12345) const;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the \p node.
   */
  void dof_indices (const Node * const node,
                    std::vector<dof_id_type> & di) const;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the \p node, for one variable \p vn.
   */
  void dof_indices (const Node * const node,
                    std::vector<dof_id_type> & di,
                    const unsigned int vn) const override;

  /**
   * Appends to the vector \p di the global degree of freedom indices
   * for \p elem.node_ref(n), for one variable \p vn.  On hanging
   * nodes with both vertex and non-vertex DoFs, only those indices
   * which are directly supported on \p elem are included.
   */
  void dof_indices (const Elem & elem,
                    unsigned int n,
                    std::vector<dof_id_type> & di,
                    const unsigned int vn) const;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Appends to the vector \p di the old global degree of freedom
   * indices for \p elem.node_ref(n), for one variable \p vn.  On
   * hanging nodes with both vertex and non-vertex DoFs, only those
   * indices which are directly supported on \p elem are included.
   */
  void old_dof_indices (const Elem & elem,
                        unsigned int n,
                        std::vector<dof_id_type> & di,
                        const unsigned int vn) const;

#endif // LIBMESH_ENABLE_AMR

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * corresponding to the SCALAR variable vn. If old_dofs=true,
   * the old SCALAR dof indices are returned.
   *
   * \note We do not need to pass in an element since SCALARs are
   * global variables.
   */
  void SCALAR_dof_indices (std::vector<dof_id_type> & di,
                           const unsigned int vn,
                           const bool old_dofs=false) const;

  /**
   * \returns \p true if degree of freedom index \p dof_index
   * is either a local index or in the \p send_list.
   *
   * \note This is an O(logN) operation for a send_list of size N; we
   * don't cache enough information for O(1) right now.
   */
  bool semilocal_index (dof_id_type dof_index) const;

  /**
   * \returns \p true if all degree of freedom indices in \p
   * dof_indices are either local indices or in the \p send_list.
   *
   * \note This is an O(logN) operation for a send_list of size N; we
   * don't cache enough information for O(1) right now.
   */
  bool all_semilocal_indices (const std::vector<dof_id_type> & dof_indices) const;

  /**
   * \returns \p true if degree of freedom index \p dof_index
   * is a local index.
   */
  bool local_index (dof_id_type dof_index) const
  { return (dof_index >= this->first_dof()) && (dof_index < this->end_dof()); }

  /**
   * \returns \p true iff our solutions can be locally evaluated on
   * \p obj (which should be an Elem or a Node) for variable number \p
   * var_num (for all variables, if \p var_num is invalid_uint)
   */
  template <typename DofObjectSubclass>
  bool is_evaluable(const DofObjectSubclass & obj,
                    unsigned int var_num = libMesh::invalid_uint) const;

  /**
   * Allow the implicit_neighbor_dofs flag to be set programmatically.
   * This overrides the --implicit_neighbor_dofs commandline option.
   * We can use this to set the implicit neighbor dofs option differently
   * for different systems, whereas the commandline option is the same
   * for all systems.
   */
  void set_implicit_neighbor_dofs(bool implicit_neighbor_dofs);

  /**
   * Set the _verify_dirichlet_bc_consistency flag.
   */
  void set_verify_dirichlet_bc_consistency(bool val);

  /**
   * Tells other library functions whether or not this problem
   * includes coupling between dofs in neighboring cells, as can
   * currently be specified on the command line or inferred from
   * the use of all discontinuous variables.
   */
  bool use_coupled_neighbor_dofs(const MeshBase & mesh) const;

  /**
   * Builds the local element vector \p Ue from the global vector \p Ug,
   * accounting for any constrained degrees of freedom.  For an element
   * without constrained degrees of freedom this is the trivial mapping
   * \f$ Ue[i] = Ug[dof_indices[i]] \f$
   *
   * \note The user must ensure that the element vector \p Ue is
   * properly sized when calling this method.  This is because there
   * is no \p resize() method in the \p DenseVectorBase<> class.
   */
  void extract_local_vector (const NumericVector<Number> & Ug,
                             const std::vector<dof_id_type> & dof_indices,
                             DenseVectorBase<Number> & Ue) const;

  /**
   * If T == dof_id_type, counts, if T == std::vector<dof_id_type>, fills an
   * array of, those dof indices which belong to the given variable number and
   * live on the current processor.
   */
  template <typename T, std::enable_if_t<std::is_same_v<T, dof_id_type> ||
                                         std::is_same_v<T, std::vector<dof_id_type>>, int> = 0>
  void local_variable_indices(T & idx,
                              const MeshBase & mesh,
                              unsigned int var_num) const;

  /**
   * If T == dof_id_type, counts, if T == std::vector<dof_id_type>, fills an
   * array of, those dof indices which belong to the given variable number and
   * live on the current processor.
   */
  template <typename T,
            std::enable_if_t<std::is_same_v<T, dof_id_type> ||
                                 std::is_same_v<T, std::vector<dof_id_type>>,
                             int> = 0>
  void local_variable_indices(T & idx, unsigned int var_num) const
  { this->local_variable_indices(idx, this->_mesh, var_num); }

#ifdef LIBMESH_ENABLE_CONSTRAINTS

  //--------------------------------------------------------------------
  // Constraint-specific methods
  /**
   * \returns The total number of constrained degrees of freedom
   * in the problem.
   */
  dof_id_type n_constrained_dofs() const;

  /**
   * \returns The number of constrained degrees of freedom
   * on this processor.
   */
  dof_id_type n_local_constrained_dofs() const;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * \returns The total number of constrained Nodes
   * in the mesh.
   */
  dof_id_type n_constrained_nodes() const
  { return cast_int<dof_id_type>(_node_constraints.size()); }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  /**
   * Rebuilds the raw degree of freedom and DofObject constraints,
   * based on attached DirichletBoundary objects and on non-conforming
   * interface in adapted meshes.
   *
   * A time is specified for use in building time-dependent Dirichlet
   * constraints.
   */
  void create_dof_constraints (const MeshBase &, Real time=0);

  /**
   * Gathers constraint equation dependencies from other processors
   */
  void allgather_recursive_constraints (MeshBase &);

  /**
   * Sends constraint equations to constraining processors
   */
  void scatter_constraints (MeshBase &);

  /**
   * Helper function for querying about constraint equations on other
   * processors.  If any id in \p requested_dof_ids is constrained on
   * another processor, its constraint will be added on this processor
   * as well.  If \p look_for_constrainees is true, then constraints
   * will also be returned if the id appears as a constraining value
   * not just if it appears as a constrained value.
   *
   * This function operates recursively: if the constraint for a
   * constrained dof is newly added locally, then any other dofs which
   * constrain it are queried to see if they are in turn constrained,
   * and so on.
   */
  void gather_constraints (MeshBase & mesh,
                           std::set<dof_id_type> & unexpanded_dofs,
                           bool look_for_constrainees);

  /**
   * Postprocesses any constrained degrees of freedom
   * to be constrained only in terms of unconstrained dofs, then adds
   * unconstrained dofs to the send_list and prepares that for use.
   * This should be run after both system (create_dof_constraints) and
   * user constraints have all been added.
   */
  void process_constraints (MeshBase &);

  /**
   * Throw an error if we detect any constraint loops, i.e.
   * A -> B -> C -> A
   * that is, "dof A is constrained in terms of dof B which is
   * constrained in terms of dof C which is constrained in terms of
   * dof A", since these are not supported by libMesh and give
   * erroneous results if they are present.
   *
   * \note The original "cyclic constraint" terminology was
   * unfortunate since the word cyclic is used by some software to
   * indicate an actual type of rotational/angular constraint and not
   * (as here) a cyclic graph. The former nomenclature will eventually
   * be deprecated in favor of "constraint loop".
   */
  void check_for_cyclic_constraints();
  void check_for_constraint_loops();

  /**
   * Adds a copy of the user-defined row to the constraint matrix, using
   * an inhomogeneous right-hand-side for the constraint equation.
   */
  void add_constraint_row (const dof_id_type dof_number,
                           const DofConstraintRow & constraint_row,
                           const Number constraint_rhs,
                           const bool forbid_constraint_overwrite);

  /**
   * Adds a copy of the user-defined row to the constraint matrix,
   * using an inhomogeneous right-hand-side for the adjoint constraint
   * equation.
   *
   * \p forbid_constraint_overwrite here only tests for overwriting
   * the rhs.  This method should only be used when an equivalent
   * constraint (with a potentially different rhs) already exists for
   * the primal problem.
   */
  void add_adjoint_constraint_row (const unsigned int qoi_index,
                                   const dof_id_type dof_number,
                                   const DofConstraintRow & constraint_row,
                                   const Number constraint_rhs,
                                   const bool forbid_constraint_overwrite);

  /**
   * Adds a copy of the user-defined row to the constraint matrix, using
   * a homogeneous right-hand-side for the constraint equation.
   * By default, produces an error if the DOF was already constrained.
   */
  void add_constraint_row (const dof_id_type dof_number,
                           const DofConstraintRow & constraint_row,
                           const bool forbid_constraint_overwrite = true)
  { add_constraint_row(dof_number, constraint_row, 0., forbid_constraint_overwrite); }

  /**
   * \returns An iterator pointing to the first DoF constraint row.
   */
  DofConstraints::const_iterator constraint_rows_begin() const
  { return _dof_constraints.begin(); }

  /**
   * \returns An iterator pointing just past the last DoF constraint row.
   */
  DofConstraints::const_iterator constraint_rows_end() const
  { return _dof_constraints.end(); }

  /**
   * Provide a const accessor to the DofConstraints map. This allows the user
   * to quickly search the data structure rather than just iterating over it.
   */
  const DofConstraints & get_dof_constraints() const { return _dof_constraints; }

  void stash_dof_constraints()
  {
    libmesh_assert(_stashed_dof_constraints.empty());
    _dof_constraints.swap(_stashed_dof_constraints);
  }

  void unstash_dof_constraints()
  {
    libmesh_assert(_dof_constraints.empty());
    _dof_constraints.swap(_stashed_dof_constraints);
  }

  /**
   * Similar to the stash/unstash_dof_constraints() API, but swaps
   * _dof_constraints and _stashed_dof_constraints without asserting
   * that the source or destination is empty first.
   *
   * \note There is an implicit assumption that swapping between sets
   * of Constraints does not change the sparsity pattern or expand the
   * send_list, since the only thing changed is the DofConstraints
   * themselves.  This is intended to work for swapping between
   * DofConstraints A and B, where A is used to define the send_list,
   * and B is a subset of A.
   */
  void swap_dof_constraints()
  {
    _dof_constraints.swap(_stashed_dof_constraints);
  }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * \returns An iterator pointing to the first Node constraint row.
   */
  NodeConstraints::const_iterator node_constraint_rows_begin() const
  { return _node_constraints.begin(); }

  /**
   * \returns An iterator pointing just past the last Node constraint row.
   */
  NodeConstraints::const_iterator node_constraint_rows_end() const
  { return _node_constraints.end(); }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  /**
   * \returns \p true if the degree of freedom dof is constrained,
   * \p false otherwise.
   */
  bool is_constrained_dof (const dof_id_type dof) const;

  /**
   * \returns \p true if the system has any heterogeneous constraints for
   * adjoint solution \p qoi_num, \p false otherwise.
   */
  bool has_heterogeneous_adjoint_constraints (const unsigned int qoi_num) const;

  /**
   * Backwards compatibility with misspelling.
   */
  bool has_heterogenous_adjoint_constraints (const unsigned int qoi_num) const
  {
    return this->has_heterogeneous_adjoint_constraints (qoi_num);
  }

  /**
   * \returns The heterogeneous constraint value if the degree of
   * freedom \p dof has a heterogeneous constraint for adjoint solution
   * \p qoi_num, zero otherwise.
   */
  Number has_heterogeneous_adjoint_constraint (const unsigned int qoi_num,
                                               const dof_id_type dof) const;

  /**
   * Backwards compatibility with misspelling.
   */
  Number has_heterogenous_adjoint_constraint (const unsigned int qoi_num,
                                              const dof_id_type dof) const
  {
    return this->has_heterogeneous_adjoint_constraint (qoi_num, dof);
  }

  /**
   * \returns A reference to the set of right-hand-side values in
   * primal constraint equations
   */
  DofConstraintValueMap & get_primal_constraint_values();

  /**
   * \returns \p true if the Node is constrained,
   * false otherwise.
   */
  bool is_constrained_node (const Node * node) const;

  /**
   * Prints (from processor 0) all DoF and Node constraints.  If \p
   * print_nonlocal is true, then each constraint is printed once for
   * each processor that knows about it, which may be useful for \p
   * DistributedMesh debugging.
   */
  void print_dof_constraints(std::ostream & os=libMesh::out,
                             bool print_nonlocal=false) const;

  /**
   * Gets a string reporting all DoF and Node constraints local to
   * this processor.  If \p print_nonlocal is true, then nonlocal
   * constraints which are locally known are included.
   */
  std::string get_local_constraints(bool print_nonlocal=false) const;


  /**
   * Tests the constrained degrees of freedom on the numeric vector \p v, which
   * represents a solution defined on the mesh, returning a pair whose first
   * entry is the maximum absolute error on a constrained DoF and whose second
   * entry is the maximum relative error.  Useful for debugging purposes.
   *
   * If \p v == nullptr, the system solution vector is tested.
   */
  std::pair<Real, Real> max_constraint_error(const System & system,
                                             NumericVector<Number> * v = nullptr) const;

#endif // LIBMESH_ENABLE_CONSTRAINTS

  //--------------------------------------------------------------------
  // Constraint-specific methods
  // Some of these methods are enabled (but inlined away to nothing)
  // when constraints are disabled at configure-time.  This is to
  // increase API compatibility of user code with different library
  // builds.

  /**
   * Constrains the element matrix.  This method requires the
   * element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   *
   * If \p asymmetric_constraint_rows is set to true (as it is by
   * default), constraint row equations will be reinforced in a way
   * which breaks matrix symmetry but makes inexact linear solver
   * solutions more likely to satisfy hanging node constraints.
   */
  void constrain_element_matrix (DenseMatrix<Number> & matrix,
                                 std::vector<dof_id_type> & elem_dofs,
                                 bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element matrix.  This method allows the
   * element matrix to be non-square, in which case the row_dofs
   * and col_dofs may be of different size and correspond to
   * variables approximated in different spaces.
   */
  void constrain_element_matrix (DenseMatrix<Number> & matrix,
                                 std::vector<dof_id_type> & row_dofs,
                                 std::vector<dof_id_type> & col_dofs,
                                 bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element vector.
   */
  void constrain_element_vector (DenseVector<Number> & rhs,
                                 std::vector<dof_id_type> & dofs,
                                 bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element matrix and vector.  This method requires
   * the element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   */
  void constrain_element_matrix_and_vector (DenseMatrix<Number> & matrix,
                                            DenseVector<Number> & rhs,
                                            std::vector<dof_id_type> & elem_dofs,
                                            bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element matrix and vector.  This method requires
   * the element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   *
   * The heterogeneous version of this method creates linear systems in
   * which heterogeneously constrained degrees of freedom will solve to
   * their correct offset values, as would be appropriate for finding
   * a solution to a linear problem in a single algebraic solve.  The
   * non-heterogeneous version of this method creates linear systems in
   * which even heterogeneously constrained degrees of freedom are
   * solved without offset values taken into account, as would be
   * appropriate for finding linearized updates to a solution in which
   * heterogeneous constraints are already satisfied.
   *
   * By default, the constraints for the primal solution of this
   * system are used.  If a non-negative \p qoi_index is passed in,
   * then the constraints for the corresponding adjoint solution are
   * used instead.
   */
  void heterogeneously_constrain_element_matrix_and_vector (DenseMatrix<Number> & matrix,
                                                            DenseVector<Number> & rhs,
                                                            std::vector<dof_id_type> & elem_dofs,
                                                            bool asymmetric_constraint_rows = true,
                                                            int qoi_index = -1) const;

  /*
   * Backwards compatibility with misspelling.
   */
  void heterogenously_constrain_element_matrix_and_vector (DenseMatrix<Number> & matrix,
                                                           DenseVector<Number> & rhs,
                                                           std::vector<dof_id_type> & elem_dofs,
                                                           bool asymmetric_constraint_rows = true,
                                                           int qoi_index = -1) const
  {
    return this->heterogeneously_constrain_element_matrix_and_vector
      (matrix, rhs, elem_dofs, asymmetric_constraint_rows, qoi_index);
  }

  /**
   * Constrains the element vector.  This method requires
   * the element matrix to be square and not-yet-constrained, in which
   * case the elem_dofs correspond to the global DOF indices of both
   * the rows and columns of the element matrix.
   *
   * The heterogeneous version of this method creates linear systems in
   * which heterogeneously constrained degrees of freedom will solve to
   * their correct offset values, as would be appropriate for finding
   * a solution to a linear problem in a single algebraic solve.  The
   * non-heterogeneous version of this method creates linear systems in
   * which even heterogeneously constrained degrees of freedom are
   * solved without offset values taken into account, as would be
   * appropriate for finding linearized updates to a solution in which
   * heterogeneous constraints are already satisfied.
   *
   * Note the sign difference from the nonlinear heterogeneous constraint
   * method: Solving u:=K\f has the opposite sign convention from
   * u:=u_in-J\r, and we apply heterogeneous constraints accordingly.
   *
   * By default, the constraints for the primal solution of this
   * system are used.  If a non-negative \p qoi_index is passed in,
   * then the constraints for the corresponding adjoint solution are
   * used instead.
   */
  void heterogeneously_constrain_element_vector (const DenseMatrix<Number> & matrix,
                                                 DenseVector<Number> & rhs,
                                                 std::vector<dof_id_type> & elem_dofs,
                                                 bool asymmetric_constraint_rows = true,
                                                 int qoi_index = -1) const;

  /*
   * Backwards compatibility with misspelling.
   */
  void heterogenously_constrain_element_vector (const DenseMatrix<Number> & matrix,
                                                DenseVector<Number> & rhs,
                                                std::vector<dof_id_type> & elem_dofs,
                                                bool asymmetric_constraint_rows = true,
                                                int qoi_index = -1) const
  {
    return this->heterogeneously_constrain_element_vector
      (matrix, rhs, elem_dofs, asymmetric_constraint_rows, qoi_index);
  }

  /**
   * Constrains the element Jacobian and residual.  The element
   * Jacobian is square, and the elem_dofs should correspond to the
   * global DOF indices of both the rows and columns of the element
   * matrix.
   *
   * The residual-constraining version of this method creates linear
   * systems in which heterogeneously constrained degrees of freedom
   * create non-zero residual terms when not at their correct offset
   * values, as would be appropriate for finding a solution to a
   * nonlinear problem in a quasi-Newton solve.
   *
   * Note the sign difference from the linear heterogeneous constraint
   * method: Solving u:=u_in-J\r has the opposite sign convention from
   * u:=K\f, and we apply heterogeneous constraints accordingly.
   *
   * The \p solution vector passed in should be a serialized or
   * ghosted primal solution
   */
  void heterogeneously_constrain_element_jacobian_and_residual (DenseMatrix<Number> & matrix,
                                                                DenseVector<Number> & rhs,
                                                                std::vector<dof_id_type> & elem_dofs,
                                                                NumericVector<Number> & solution_local) const;

  /**
   * Constrains the element residual.  The element Jacobian is square,
   * and the elem_dofs should correspond to the global DOF indices of
   * both the rows and columns of the element matrix.
   *
   * The residual-constraining version of this method creates linear
   * systems in which heterogeneously constrained degrees of freedom
   * create non-zero residual terms when not at their correct offset
   * values, as would be appropriate for finding a solution to a
   * nonlinear problem in a quasi-Newton solve.
   *
   * The \p solution vector passed in should be a serialized or
   * ghosted primal solution
   */
  void heterogeneously_constrain_element_residual (DenseVector<Number> & rhs,
                                                   std::vector<dof_id_type> & elem_dofs,
                                                   NumericVector<Number> & solution_local) const;


  /**
   * Constrains the element residual.  The element Jacobian is square,
   * and the elem_dofs should correspond to the global DOF indices of
   * both the rows and columns of the element matrix, and the dof
   * constraint should not include any heterogeneous terms.
   *
   * The residual-constraining version of this method creates linear
   * systems in which heterogeneously constrained degrees of freedom
   * create non-zero residual terms when not at their correct offset
   * values, as would be appropriate for finding a solution to a
   * nonlinear problem in a quasi-Newton solve.
   *
   * The \p solution vector passed in should be a serialized or
   * ghosted primal solution
   */
  void constrain_element_residual (DenseVector<Number> & rhs,
                                   std::vector<dof_id_type> & elem_dofs,
                                   NumericVector<Number> & solution_local) const;

  /**
   * Constrains a dyadic element matrix B = v w'.  This method
   * requires the element matrix to be square, in which case the
   * elem_dofs correspond to the global DOF indices of both the rows
   * and columns of the element matrix.  For this case the rows and
   * columns of the matrix necessarily correspond to variables of the
   * same approximation order.
   */
  void constrain_element_dyad_matrix (DenseVector<Number> & v,
                                      DenseVector<Number> & w,
                                      std::vector<dof_id_type> & row_dofs,
                                      bool asymmetric_constraint_rows = true) const;

  /**
   * Does not actually constrain anything, but modifies \p dofs in the
   * same way as any of the constrain functions would do, i.e. adds
   * those dofs in terms of which any of the existing dofs is
   * constrained.
   */
  void constrain_nothing (std::vector<dof_id_type> & dofs) const;

  /**
   * Constrains the numeric vector \p v, which represents a solution defined on
   * the mesh.  This may need to be used after a linear solve, if your linear
   * solver's solutions do not satisfy your DoF constraints to a tight enough
   * tolerance.
   *
   * If \p v == nullptr, the system solution vector is constrained
   *
   * If \p homogeneous == true, heterogeneous constraints are enforced
   * as if they were homogeneous.  This might be appropriate for e.g. a
   * vector representing a difference between two
   * heterogeneously-constrained solutions.
   */
  void enforce_constraints_exactly (const System & system,
                                    NumericVector<Number> * v = nullptr,
                                    bool homogeneous = false) const;

  /**
   * Heterogeneously constrains the numeric vector \p v, which
   * represents an adjoint solution defined on the mesh for quantity
   * fo interest \p q.  For homogeneous constraints, use \p
   * enforce_constraints_exactly instead
   */
  void enforce_adjoint_constraints_exactly (NumericVector<Number> & v,
                                            unsigned int q) const;

  void enforce_constraints_on_residual (const NonlinearImplicitSystem & system,
                                        NumericVector<Number> * rhs,
                                        NumericVector<Number> const * solution,
                                        bool homogeneous = true) const;

  void enforce_constraints_on_jacobian (const NonlinearImplicitSystem & system,
                                        SparseMatrix<Number> * jac) const;

#ifdef LIBMESH_ENABLE_PERIODIC

  //--------------------------------------------------------------------
  // PeriodicBoundary-specific methods

  /**
   * Adds a copy of the specified periodic boundary to the system.
   */
  void add_periodic_boundary (const PeriodicBoundaryBase & periodic_boundary);

  /**
   * Add a periodic boundary pair
   *
   * \param boundary - primary boundary
   * \param inverse_boundary - inverse boundary
   */
  void add_periodic_boundary (const PeriodicBoundaryBase & boundary, const PeriodicBoundaryBase & inverse_boundary);

  /**
   * \returns \p true if the boundary given by \p boundaryid is periodic,
   * false otherwise
   */
  bool is_periodic_boundary (const boundary_id_type boundaryid) const;

  PeriodicBoundaries * get_periodic_boundaries()
  {
    return _periodic_boundaries.get();
  }

  const PeriodicBoundaries * get_periodic_boundaries() const
  {
    return _periodic_boundaries.get();
  }

#endif // LIBMESH_ENABLE_PERIODIC


#ifdef LIBMESH_ENABLE_DIRICHLET

  //--------------------------------------------------------------------
  // DirichletBoundary-specific methods

  /**
   * Adds a copy of the specified Dirichlet boundary to the system.
   *
   * The constraints implied by DirichletBoundary objects are imposed
   * in the same order in which DirichletBoundary objects are added to
   * the DofMap. When multiple DirichletBoundary objects would impose
   * competing constraints on a given DOF, the *first*
   * DirichletBoundary to constrain the DOF "wins". This distinction
   * is important when e.g. two surfaces (sidesets) intersect. The
   * nodes on the intersection will be constrained according to
   * whichever sideset's DirichletBoundary object was added to the
   * DofMap first.
   */
  void add_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary);

  /**
   * Adds a copy of the specified Dirichlet boundary to the system,
   * corresponding to the adjoint problem defined by Quantity of
   * Interest \p q.
   */
  void add_adjoint_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary,
                                       unsigned int q);

  /**
   * Removes the specified Dirichlet boundary from the system.
   */
  void remove_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary);

  /**
   * Removes from the system the specified Dirichlet boundary for the
   * adjoint equation defined by Quantity of interest index q
   */
  void remove_adjoint_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary,
                                          unsigned int q);

  const DirichletBoundaries * get_dirichlet_boundaries() const
  {
    return _dirichlet_boundaries.get();
  }

  DirichletBoundaries * get_dirichlet_boundaries()
  {
    return _dirichlet_boundaries.get();
  }

  bool has_adjoint_dirichlet_boundaries(unsigned int q) const;

  const DirichletBoundaries *
  get_adjoint_dirichlet_boundaries(unsigned int q) const;

  DirichletBoundaries *
  get_adjoint_dirichlet_boundaries(unsigned int q);

  /**
   * Check that all the ids in dirichlet_bcids are actually present in the mesh.
   * If not, this will throw an error.
   */
  void check_dirichlet_bcid_consistency (const MeshBase & mesh,
                                         const DirichletBoundary & boundary) const;
#endif // LIBMESH_ENABLE_DIRICHLET


#ifdef LIBMESH_ENABLE_AMR

  //--------------------------------------------------------------------
  // AMR-specific methods

  /**
   * After a mesh is refined and repartitioned it is possible that the
   * \p _send_list will need to be augmented.  This is the case when an
   * element is refined and its children end up on different processors
   * than the parent.  These children will need values from the parent
   * when projecting the solution onto the refined mesh, hence the parent's
   * DOF indices need to be included in the \p _send_list.
   */
  // void augment_send_list_for_projection(const MeshBase &);

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Fills the vector di with the global degree of freedom indices
   * for the element using the \p DofMap::old_dof_object.
   * If no variable number is specified then all
   * variables are returned.
   */
  void old_dof_indices (const Elem * const elem,
                        std::vector<dof_id_type> & di,
                        const unsigned int vn = libMesh::invalid_uint) const;

#endif // LIBMESH_ENABLE_AMR

  /**
   * Constrains degrees of freedom on side \p s of element \p elem which
   * correspond to variable number \p var and to p refinement levels
   * above \p p.
   */
  void constrain_p_dofs (unsigned int var,
                         const Elem * elem,
                         unsigned int s,
                         unsigned int p);

#endif // LIBMESH_ENABLE_AMR

  /**
   * Reinitialize the underlying data structures conformal to the current mesh.
   */
  void reinit
    (MeshBase & mesh,
     const std::map<const Node *, std::set<subdomain_id_type>> &
       constraining_subdomains);

  /**
   * Free all new memory associated with the object, but restore its
   * original state, with the mesh pointer and any default ghosting.
   */
  virtual void clear () override;

  /**
   * Prints summary info about the sparsity bandwidth and constraints.
   */
  void print_info(std::ostream & os=libMesh::out) const;

  /**
   * Gets summary info about the sparsity bandwidth and constraints.
   */
  std::string get_info() const;

  /**
   * Degree of freedom coupling.  If left empty each DOF
   * couples to all others.  Can be used to reduce memory
   * requirements for sparse matrices.  DOF 0 might only
   * couple to itself, in which case \p dof_coupling(0,0)
   * should be 1 and \p dof_coupling(0,j) = 0 for j not equal
   * to 0.
   *
   * This variable is named as though it were class private,
   * but it is in the public interface.  Also there are no
   * public methods for accessing it...  This typically means
   * you should only use it if you know what you are doing.
   */
  CouplingMatrix * _dof_coupling;

  /**
   * \returns The number of the system we are responsible for.
   */
  unsigned int sys_number() const;

  /**
   * Builds a sparsity pattern for matrices using the current
   * degree-of-freedom numbering and coupling.
   *
   * By default, ignores constraint equations, for build speed; this
   * is valid for the combination of !need_full_sparsity_pattern and
   * constraints which only come from periodic boundary conditions and
   * adaptive mesh refinement, where matrix constraint adds some
   * matrix entries but removes equally many (or more) other entries.
   *
   * Can be told to calculate sparsity for the constrained matrix,
   * which may be necessary in the case of spline control node
   * constraints or sufficiently many user constraints.
   */
  std::unique_ptr<SparsityPattern::Build> build_sparsity(const MeshBase & mesh,
                                                         bool calculate_constrained = false,
                                                         bool use_condensed_system = false) const;

  /**
   * Describe whether the given variable group should be p-refined. If this API is not called with
   * \p false, the default is to p-refine
   */
  void should_p_refine(unsigned int g, bool p_refine);

  /**
   * Whether the given variable group should be p-refined
   */
  bool should_p_refine(unsigned int g) const;

  /**
   * Whether the given variable should be p-refined
   */
  bool should_p_refine_var(unsigned int var) const;

  // Prevent bad user implicit conversions
  void should_p_refine(FEFamily, bool) = delete;
  void should_p_refine(Order, bool) = delete;
  bool should_p_refine(FEFamily) const = delete;
  bool should_p_refine(Order) const = delete;

  /**
   * Add a static condensation class
   */
  void create_static_condensation(MeshBase & mesh, System & system);

  /**
   * Checks whether we have static condensation
   */
  bool has_static_condensation() const { return _sc.get(); }

  /**
   * @returns the static condensation class. This should have been already added with a call to \p
   * add_static_condensation()
   */
  StaticCondensationDofMap & get_static_condensation();

  /**
   * Calls reinit on the static condensation map if it exists
   */
  void reinit_static_condensation();

private:

  /**
   * Helper function that gets the dof indices on the current element
   * for a non-SCALAR type variable, where the variable is identified
   * by its variable group number \p vg and its offset \p vig from the
   * first variable in that group.
   *
   * In DEBUG mode, the tot_size parameter will add up the total
   * number of dof indices that should have been added to di, and v
   * will be the variable number corresponding to vg and vig.
   */
  void _dof_indices (const Elem & elem,
                     int p_level,
                     std::vector<dof_id_type> & di,
                     const unsigned int vg,
                     const unsigned int vig,
                     const Node * const * nodes,
                     unsigned int       n_nodes,
                     const unsigned int v
#ifdef DEBUG
                     ,
                     std::size_t & tot_size
#endif
                    ) const;

  /**
   * As above except a \p field_dofs_functor must be provided. This method is useful when the caller
   * wants to do more than simply fill a degree of freedom container
   * @param field_dofs_functor This functor has the interface:
   *                           void field_dofs_functor(const Elem & elem,
   *                                                   const unsigned int node_num,
   *                                                   const unsigned int var_num,
   *                                                   std::vector<dof_id_type> & di,
   *                                                   const dof_id_type field_dof)
   *                           where \p field_dof represents a field degree of freedom to act on and
   *                           is associated with \p node_num and \p var_num. If the degree of
   *                           freedom is elemental than \p node_num will be \p invalid_uint. \p di
   *                           is the degree of freedom container provided to the \p _dof_indices
   *                           method
   */
  template <typename FieldDofsFunctor>
  void _dof_indices (const Elem & elem,
                     int p_level,
                     std::vector<dof_id_type> & di,
                     const unsigned int vg,
                     const unsigned int vig,
                     const Node * const * nodes,
                     unsigned int       n_nodes,
                     const unsigned int v,
#ifdef DEBUG
                     std::size_t & tot_size,
#endif
                     FieldDofsFunctor field_dofs_functor) const;

  /**
   * Helper function that implements the element-nodal versions of
   * dof_indices and old_dof_indices
   */
  void _node_dof_indices (const Elem & elem,
                          unsigned int n,
                          const DofObject & obj,
                          std::vector<dof_id_type> & di,
                          const unsigned int vn) const;

  /**
   * Invalidates all active DofObject dofs for this system
   */
  void invalidate_dofs(MeshBase & mesh) const;

  /**
   * \returns The Node pointer with index \p i from the \p mesh.
   */
  DofObject * node_ptr(MeshBase & mesh, dof_id_type i) const;

  /**
   * \returns The Elem pointer with index \p i from the \p mesh.
   */
  DofObject * elem_ptr(MeshBase & mesh, dof_id_type i) const;

  /**
   * A member function type like \p node_ptr() or \p elem_ptr().
   */
  typedef DofObject * (DofMap::*dofobject_accessor)
    (MeshBase & mesh, dof_id_type i) const;

  /**
   * Helper function for distributing dofs in parallel
   */
  template<typename iterator_type>
  void set_nonlocal_dof_objects(iterator_type objects_begin,
                                iterator_type objects_end,
                                MeshBase & mesh,
                                dofobject_accessor objects);

  /**
   * We may have mesh constraint rows with dependent nodes in one
   * subdomain but dependency nodes in another subdomain, and we may
   * have variables whose subdomain restriction includes the dependent
   * subdomain but not the dependency.  In those cases we need to
   * place degrees of freedom on dependency nodes anyway.
   *
   * The set value for node n will include all subdomain ids of
   * elements with nodes in subdomains constrained by n.
   *
   * We use a map<set> rather than a multimap here because we expect
   * to be inserting the same subdomain multiple times and we don't
   * need duplicate values.
   */
  std::map<const Node *, std::set<subdomain_id_type>>
  calculate_constraining_subdomains();

  /**
   * Distributes the global degrees of freedom, for dofs on
   * this processor.  In this format the local
   * degrees of freedom are in a contiguous block for each
   * variable in the system.
   * Starts at index next_free_dof, and increments it to
   * the post-final index.
   *
   * Uses the provided constraining_subdomains map from
   * calculate_constraining_subdomains() to ensure allocation of all
   * DoFs on constraining nodes.
   */
  void distribute_local_dofs_var_major
    (dof_id_type & next_free_dof,
     MeshBase & mesh,
     const std::map<const Node *, std::set<subdomain_id_type>> &
       constraining_subdomains);

  /**
   * Distributes the global degrees of freedom for dofs on this
   * processor.  In this format all the degrees of freedom at a
   * node/element are in contiguous blocks.  Starts at index \p
   * next_free_dof, and increments it to the post-final index.  If \p
   * build_send_list is \p true, builds the send list.  If \p false,
   * clears and reserves the send list.
   *
   * Uses the provided constraining_subdomains map from
   * calculate_constraining_subdomains() to ensure allocation of all
   * DoFs on constraining nodes.
   *
   * \note The degrees of freedom for a given variable are not in
   * contiguous blocks, as in the case of \p distribute_local_dofs_var_major.
   */
  void distribute_local_dofs_node_major
    (dof_id_type & next_free_dof,
     MeshBase & mesh,
     const std::map<const Node *, std::set<subdomain_id_type>> &
       constraining_subdomains);

  /*
   * Helper method for the above two to count + distriubte SCALAR dofs
   */
  void distribute_scalar_dofs (dof_id_type & next_free_dof);

#ifdef DEBUG
  /*
   * Internal assertions for distribute_local_dofs_*
   */
  void assert_no_nodes_missed(MeshBase & mesh);
#endif

  /*
   * A utility method for obtaining a set of elements to ghost along
   * with merged coupling matrices.
   */
  typedef std::set<std::unique_ptr<CouplingMatrix>, Utility::CompareUnderlying> CouplingMatricesSet;
  static void
  merge_ghost_functor_outputs (GhostingFunctor::map_type & elements_to_ghost,
                               CouplingMatricesSet & temporary_coupling_matrices,
                               const GhostingFunctorIterator & gf_begin,
                               const GhostingFunctorIterator & gf_end,
                               const MeshBase::const_element_iterator & elems_begin,
                               const MeshBase::const_element_iterator & elems_end,
                               processor_id_type p);

  /**
   * Adds entries to the \p _send_list vector corresponding to DoFs
   * on elements neighboring the current processor.
   */
  void add_neighbors_to_send_list(MeshBase & mesh);

#ifdef LIBMESH_ENABLE_CONSTRAINTS

  /**
   * Build the constraint matrix C associated with the element
   * degree of freedom indices elem_dofs. The optional parameter
   * \p called_recursively should be left at the default value
   * \p false.  This is used to handle the special case of
   * an element's degrees of freedom being constrained in terms
   * of other, local degrees of freedom.  The usual case is
   * for an elements DOFs to be constrained by some other,
   * external DOFs.
   */
  void build_constraint_matrix (DenseMatrix<Number> & C,
                                std::vector<dof_id_type> & elem_dofs,
                                const bool called_recursively=false) const;

  /**
   * Build the constraint matrix C and the forcing vector H
   * associated with the element degree of freedom indices elem_dofs.
   * The optional parameter \p called_recursively should be left at
   * the default value \p false.  This is used to handle the special
   * case of an element's degrees of freedom being constrained in
   * terms of other, local degrees of freedom.  The usual case is for
   * an elements DOFs to be constrained by some other, external DOFs
   * and/or Dirichlet conditions.
   *
   * The forcing vector will depend on which solution's heterogeneous
   * constraints are being applied.  For the default \p qoi_index this
   * will be the primal solution; for \p qoi_index >= 0 the
   * corresponding adjoint solution's constraints will be used.
   */
  void build_constraint_matrix_and_vector (DenseMatrix<Number> & C,
                                           DenseVector<Number> & H,
                                           std::vector<dof_id_type> & elem_dofs,
                                           int qoi_index = -1,
                                           const bool called_recursively=false) const;

  /**
   * Finds all the DOFS associated with the element DOFs elem_dofs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dofs (std::vector<dof_id_type> & elem_dofs) const;

  /**
   * Finds all the DofObjects associated with the set in \p objs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dof_objects (std::vector<const DofObject *> & objs) const;

  /**
   * Adds entries to the \p _send_list vector corresponding to DoFs
   * which are dependencies for constraint equations on the current
   * processor.
   */
  void add_constraints_to_send_list();

  /**
   * Adds any spline constraints from the Mesh to our DoF constraints.
   * If any Dirichlet constraints exist on spline-constrained nodes,
   * l2-projects those constraints onto the spline basis.
   */
  void process_mesh_constraint_rows(const MeshBase & mesh);

#endif // LIBMESH_ENABLE_CONSTRAINTS

  /**
   * This flag indicates whether or not we do an opt-mode check for
   * the presence of constraint loops, i.e. cases where the constraint
   * graph is cyclic.
   */
  bool _error_on_constraint_loop;

  /**
   * This flag indicates whether or not we explicitly take constraint
   * equations into account when computing a sparsity pattern.
   */
  bool _constrained_sparsity_construction;

  /**
   * The finite element type for each variable.
   */
  std::vector<Variable> _variables;

  /**
   * The finite element type for each variable group.
   */
  std::vector<VariableGroup> _variable_groups;

  /**
   * The variable group number for each variable.
   */
  std::vector<unsigned int> _variable_group_numbers;

  /**
   * A map from variable number to variable group number
   */
  std::unordered_map<unsigned int, unsigned int> _var_to_vg;

  /**
   * The number of the system we manage DOFs for.
   */
  const unsigned int _sys_number;

  /**
   * The mesh that system uses.
   */
  MeshBase & _mesh;

  /**
   * Additional matrices handled by this object.  These pointers do \e
   * not handle the memory, instead, \p System, who
   * told \p DofMap about them, owns them.
   */
  std::vector<SparseMatrix<Number> * > _matrices;

  /**
   * First DOF index for SCALAR variable v, or garbage for non-SCALAR
   * variable v
   */
  std::vector<dof_id_type> _first_scalar_df;

  /**
   * A list containing all the global DOF indices that affect the
   * solution on my processor.
   */
  std::vector<dof_id_type> _send_list;

  /**
   * Function object to call to add extra entries to the sparsity pattern
   */
  SparsityPattern::AugmentSparsityPattern * _augment_sparsity_pattern;

  /**
   * A function pointer to a function to call to add extra entries to the sparsity pattern
   */
  void (*_extra_sparsity_function)(SparsityPattern::Graph &,
                                   std::vector<dof_id_type> & n_nz,
                                   std::vector<dof_id_type> & n_oz,
                                   void *);
  /**
   * A pointer associated with the extra sparsity that can optionally be passed in
   */
  void * _extra_sparsity_context;

  /**
   * Function object to call to add extra entries to the send list
   */
  AugmentSendList * _augment_send_list;

  /**
   * A function pointer to a function to call to add extra entries to the send list
   */
  void (*_extra_send_list_function)(std::vector<dof_id_type> &, void *);

  /**
   * A pointer associated with the extra send list that can optionally be passed in
   */
  void * _extra_send_list_context;

  /**
   * The default coupling GhostingFunctor, used to implement standard
   * libMesh sparsity pattern construction.
   *
   * We use a std::unique_ptr here to reduce header dependencies.
   */
  std::unique_ptr<DefaultCoupling> _default_coupling;

  /**
   * The default algebraic GhostingFunctor, used to implement standard
   * libMesh send_list construction.
   *
   * We use a std::unique_ptr here to reduce header dependencies.
   */
  std::unique_ptr<DefaultCoupling> _default_evaluating;

  /**
   * The list of all GhostingFunctor objects to be used when
   * distributing ghosted vectors.
   *
   * The library should automatically refer these functors to the
   * MeshBase, too, so any algebraically ghosted dofs will live on
   * geometrically ghosted elements.
   *
   * Keep these in a vector so any parallel computation is done in the
   * same order on all processors.
   */
  std::vector<GhostingFunctor *> _algebraic_ghosting_functors;

  /**
   * The list of all GhostingFunctor objects to be used when
   * coupling degrees of freedom in matrix sparsity patterns.
   *
   * These objects will *also* be used as algebraic ghosting functors,
   * but not vice-versa.
   *
   * The library should automatically refer these functors to the
   * MeshBase, too, so any dofs coupled to local dofs will live on
   * geometrically ghosted elements.
   */
  std::vector<GhostingFunctor *> _coupling_functors;

  /**
   * Hang on to references to any GhostingFunctor objects we were
   * passed in shared_ptr form
   */
  std::map<GhostingFunctor *, std::shared_ptr<GhostingFunctor> > _shared_functors;

  /**
   * Default false; set to true if any attached matrix requires a full
   * sparsity pattern.
   */
  bool need_full_sparsity_pattern;

  /**
   * The sparsity pattern of the global matrix.  If
   * need_full_sparsity_pattern is true, we save the entire sparse
   * graph here.  Otherwise we save just the n_nz and n_oz vectors.
   */
  std::unique_ptr<SparsityPattern::Build> _sp;

  /**
   * The total number of SCALAR dofs associated to
   * all SCALAR variables.
   */
  dof_id_type _n_SCALAR_dofs;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * First old DOF index for SCALAR variable v, or garbage for
   * non-SCALAR variable v
   */
  std::vector<dof_id_type> _first_old_scalar_df;

  /**
   * A container of variable groups that we should not p-refine
   */
  std::unordered_set<unsigned int> _dont_p_refine;
#endif

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  /**
   * Data structure containing DOF constraints.  The ith
   * entry is the constraint matrix row for DOF i.
   */
  DofConstraints _dof_constraints, _stashed_dof_constraints;

  DofConstraintValueMap      _primal_constraint_values;

  AdjointDofConstraintValues _adjoint_constraint_values;
#endif

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Data structure containing DofObject constraints.
   */
  NodeConstraints _node_constraints;
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS


#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * Data structure containing periodic boundaries.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  std::unique_ptr<PeriodicBoundaries> _periodic_boundaries;
#endif

#ifdef LIBMESH_ENABLE_DIRICHLET
  /**
   * Data structure containing Dirichlet functions.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  std::unique_ptr<DirichletBoundaries> _dirichlet_boundaries;

  /**
   * Data structure containing Dirichlet functions.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  std::vector<std::unique_ptr<DirichletBoundaries>> _adjoint_dirichlet_boundaries;
#endif

  friend class SparsityPattern::Build;

  /**
   * Bools to indicate if we override the --implicit_neighbor_dofs
   * commandline options.
   */
  bool _implicit_neighbor_dofs_initialized;
  bool _implicit_neighbor_dofs;

  /**
   * Flag which determines whether we should do some additional
   * checking of the consistency of the DirichletBoundary objects
   * added by the user. Defaults to true, but can be disabled in cases
   * where you only want to add DirichletBoundary objects "locally"
   * and can guarantee that no repartitioning will be done, since
   * repartitioning could cause processors to own new boundary sides
   * for which they no longer have the proper DirichletBoundary
   * objects stored.
   */
  bool _verify_dirichlet_bc_consistency;

  /// Static condensation class
  std::unique_ptr<StaticCondensationDofMap> _sc;
};


// ------------------------------------------------------------
// Dof Map inline member functions
inline
unsigned int DofMap::sys_number() const
{
  return _sys_number;
}



inline
const VariableGroup & DofMap::variable_group (const unsigned int g) const
{
  libmesh_assert_less (g, _variable_groups.size());

  return _variable_groups[g];
}



inline
const Variable & DofMap::variable (const unsigned int c) const
{
  libmesh_assert_less (c, _variables.size());

  return _variables[c];
}



inline
Order DofMap::variable_order (const unsigned int c) const
{
  libmesh_assert_less (c, _variables.size());

  return _variables[c].type().order;
}



inline
Order DofMap::variable_group_order (const unsigned int vg) const
{
  libmesh_assert_less (vg, _variable_groups.size());

  return _variable_groups[vg].type().order;
}



inline
const FEType & DofMap::variable_type (const unsigned int c) const
{
  libmesh_assert_less (c, _variables.size());

  return _variables[c].type();
}



inline
const FEType & DofMap::variable_group_type (const unsigned int vg) const
{
  libmesh_assert_less (vg, _variable_groups.size());

  return _variable_groups[vg].type();
}


#ifdef LIBMESH_ENABLE_CONSTRAINTS


inline
bool DofMap::is_constrained_node (const Node *
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
                                  node
#endif
                                  ) const
{
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  if (_node_constraints.count(node))
    return true;
#endif

  return false;
}


inline
bool DofMap::is_constrained_dof (const dof_id_type dof) const
{
  if (_dof_constraints.count(dof))
    return true;

  return false;
}


inline
bool DofMap::has_heterogeneous_adjoint_constraints (const unsigned int qoi_num) const
{
  AdjointDofConstraintValues::const_iterator it =
    _adjoint_constraint_values.find(qoi_num);
  if (it == _adjoint_constraint_values.end())
    return false;
  if (it->second.empty())
    return false;

  return true;
}


inline
Number DofMap::has_heterogeneous_adjoint_constraint (const unsigned int qoi_num,
                                                    const dof_id_type dof) const
{
  AdjointDofConstraintValues::const_iterator it =
    _adjoint_constraint_values.find(qoi_num);
  if (it != _adjoint_constraint_values.end())
    {
      DofConstraintValueMap::const_iterator rhsit =
        it->second.find(dof);
      if (rhsit == it->second.end())
        return 0;
      else
        return rhsit->second;
    }

  return 0;
}



inline
DofConstraintValueMap & DofMap::get_primal_constraint_values()
{
  return _primal_constraint_values;
}



#else

//--------------------------------------------------------------------
// Constraint-specific methods get inlined into nothing if
// constraints are disabled, so there's no reason for users not to
// use them.

inline void DofMap::constrain_element_matrix (DenseMatrix<Number> &,
                                              std::vector<dof_id_type> &,
                                              bool) const {}

inline void DofMap::constrain_element_matrix (DenseMatrix<Number> &,
                                              std::vector<dof_id_type> &,
                                              std::vector<dof_id_type> &,
                                              bool) const {}

inline void DofMap::constrain_element_vector (DenseVector<Number> &,
                                              std::vector<dof_id_type> &,
                                              bool) const {}

inline void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number> &,
                                                         DenseVector<Number> &,
                                                         std::vector<dof_id_type> &,
                                                         bool) const {}

inline void DofMap::heterogeneously_constrain_element_matrix_and_vector
  (DenseMatrix<Number> &, DenseVector<Number> &,
   std::vector<dof_id_type> &, bool, int) const {}

inline void DofMap::heterogeneously_constrain_element_vector
  (const DenseMatrix<Number> &, DenseVector<Number> &,
   std::vector<dof_id_type> &, bool, int) const {}

inline void DofMap::constrain_element_dyad_matrix (DenseVector<Number> &,
                                                   DenseVector<Number> &,
                                                   std::vector<dof_id_type> &,
                                                   bool) const {}

inline void DofMap::constrain_nothing (std::vector<dof_id_type> &) const {}

inline void DofMap::enforce_constraints_exactly (const System &,
                                                 NumericVector<Number> *,
                                                 bool) const {}

inline void DofMap::enforce_adjoint_constraints_exactly (NumericVector<Number> &,
                                                         unsigned int) const {}


inline void DofMap::enforce_constraints_on_residual
  (const NonlinearImplicitSystem &,
   NumericVector<Number> *,
   NumericVector<Number> const *,
   bool) const {}

inline void DofMap::enforce_constraints_on_jacobian
  (const NonlinearImplicitSystem &,
   SparseMatrix<Number> *) const {}

#endif // LIBMESH_ENABLE_CONSTRAINTS



inline
void DofMap::set_constrained_sparsity_construction(bool use_constraints)
{
  // This got only partly finished...
  if (use_constraints)
    libmesh_not_implemented();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  _constrained_sparsity_construction = use_constraints;
#endif
  libmesh_ignore(use_constraints);
}

inline
void DofMap::full_sparsity_pattern_needed()
{
  need_full_sparsity_pattern = true;
}

inline
bool DofMap::constrained_sparsity_construction()
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  return _constrained_sparsity_construction;
#else
  return true;
#endif
}

inline
void DofMap::should_p_refine(const unsigned int g, const bool p_refine)
{
#ifdef LIBMESH_ENABLE_AMR
  if (p_refine)
  {
    auto it = _dont_p_refine.find(g);
    if (it != _dont_p_refine.end())
      _dont_p_refine.erase(it);
  }
  else
    _dont_p_refine.insert(g);
#else
  libmesh_ignore(g, p_refine);
#endif
}

inline
bool DofMap::should_p_refine(const unsigned int g) const
{
#ifdef LIBMESH_ENABLE_AMR
  return !_dont_p_refine.count(g);
#else
  libmesh_ignore(g);
  return false;
#endif
}

inline
unsigned int DofMap::var_group_from_var_number(const unsigned int var_num) const
{
  libmesh_assert(var_num < n_variables());
  return libmesh_map_find(_var_to_vg, var_num);
}

inline
bool DofMap::should_p_refine_var(const unsigned int var) const
{
#ifdef LIBMESH_ENABLE_AMR
  const auto vg = this->var_group_from_var_number(var);
  return !_dont_p_refine.count(vg);
#else
  libmesh_ignore(var);
  return false;
#endif
}

template <typename FieldDofsFunctor>
void DofMap::_dof_indices (const Elem & elem,
                           int p_level,
                           std::vector<dof_id_type> & di,
                           const unsigned int vg,
                           const unsigned int vig,
                           const Node * const * nodes,
                           unsigned int       n_nodes,
                           const unsigned int v,
#ifdef DEBUG
                           std::size_t & tot_size,
#endif
                           FieldDofsFunctor field_dofs_functor) const
{
  const VariableGroup & var = this->variable_group(vg);

  if (var.active_on_subdomain(elem.subdomain_id()))
    {
      const ElemType type        = elem.type();
      const unsigned int sys_num = this->sys_number();
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
      const bool is_inf          = elem.infinite();
#endif

      const bool extra_hanging_dofs =
        FEInterface::extra_hanging_dofs(var.type());

      FEType fe_type = var.type();

      const bool add_p_level =
#ifdef LIBMESH_ENABLE_AMR
          !_dont_p_refine.count(vg);
#else
          false;
#endif

#ifdef DEBUG
      // The number of dofs per element is non-static for subdivision FE
      if (var.type().family == SUBDIVISION)
        tot_size += n_nodes;
      else
        // FIXME: Is the passed-in p_level just elem.p_level()? If so,
        // this seems redundant.
        tot_size += FEInterface::n_dofs(fe_type, add_p_level*p_level, &elem);
#endif

      // The total Order is not required when getting the function
      // pointer, it is only needed when the function is called (see
      // below).
      const FEInterface::n_dofs_at_node_ptr ndan =
        FEInterface::n_dofs_at_node_function(fe_type, &elem);

      // Get the node-based DOF numbers
      for (unsigned int n=0; n != n_nodes; n++)
        {
          const Node & node = *nodes[n];

          // Cache the intermediate lookups that are common to every
          // component
#ifdef DEBUG
          const std::pair<unsigned int, unsigned int>
            vg_and_offset = node.var_to_vg_and_offset(sys_num,v);
          libmesh_assert_equal_to (vg, vg_and_offset.first);
          libmesh_assert_equal_to (vig, vg_and_offset.second);
#endif
          const unsigned int n_comp = node.n_comp_group(sys_num,vg);

          // There is a potential problem with h refinement.  Imagine a
          // quad9 that has a linear FE on it.  Then, on the hanging side,
          // it can falsely identify a DOF at the mid-edge node. This is why
          // we go through FEInterface instead of node.n_comp() directly.
          const unsigned int nc =
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
            is_inf ?
            FEInterface::n_dofs_at_node(fe_type, add_p_level*p_level, &elem, n) :
#endif
            ndan (type, fe_type.order + add_p_level*p_level, n);

          // If this is a non-vertex on a hanging node with extra
          // degrees of freedom, we use the non-vertex dofs (which
          // come in reverse order starting from the end, to
          // simplify p refinement)
          if (extra_hanging_dofs && !elem.is_vertex(n))
            {
              const int dof_offset = n_comp - nc;

              // We should never have fewer dofs than necessary on a
              // node unless we're getting indices on a parent element,
              // and we should never need the indices on such a node
              if (dof_offset < 0)
                {
                  libmesh_assert(!elem.active());
                  di.resize(di.size() + nc, DofObject::invalid_id);
                }
              else
                for (int i=int(n_comp)-1; i>=dof_offset; i--)
                  {
                    const dof_id_type d =
                      node.dof_number(sys_num, vg, vig, i, n_comp);
                    libmesh_assert_not_equal_to (d, DofObject::invalid_id);
                    field_dofs_functor(elem, n, v, di, d);
                  }
            }
          // If this is a vertex or an element without extra hanging
          // dofs, our dofs come in forward order coming from the
          // beginning
          else
            {
              // We have a good component index only if it's being
              // used on this FE type (nc) *and* it's available on
              // this DofObject (n_comp).
              const unsigned int good_nc = std::min(n_comp, nc);
              for (unsigned int i=0; i!=good_nc; ++i)
                {
                  const dof_id_type d =
                    node.dof_number(sys_num, vg, vig, i, n_comp);
                  libmesh_assert_not_equal_to (d, DofObject::invalid_id);
                  libmesh_assert_less (d, this->n_dofs());
                  field_dofs_functor(elem, n, v, di, d);
                }

              // With fewer good component indices than we need, e.g.
              // due to subdomain expansion, the remaining expected
              // indices are marked invalid.
              if (n_comp < nc)
                for (unsigned int i=n_comp; i!=nc; ++i)
                  di.push_back(DofObject::invalid_id);
            }
        }

      // If there are any element-based DOF numbers, get them
      const unsigned int nc = FEInterface::n_dofs_per_elem(fe_type, add_p_level*p_level, &elem);

      // We should never have fewer dofs than necessary on an
      // element unless we're getting indices on a parent element
      // (and we should never need those indices) or off-domain for a
      // subdomain-restricted variable (where invalid_id is the
      // correct thing to return)
      if (nc != 0)
        {
          const unsigned int n_comp = elem.n_comp_group(sys_num,vg);
          if (elem.n_systems() > sys_num && nc <= n_comp)
            {
              for (unsigned int i=0; i<nc; i++)
                {
                  const dof_id_type d =
                    elem.dof_number(sys_num, vg, vig, i, n_comp);
                  libmesh_assert_not_equal_to (d, DofObject::invalid_id);

                  field_dofs_functor(elem, invalid_uint, v, di, d);
                }
            }
          else
            {
              libmesh_assert(!elem.active() || fe_type.family == LAGRANGE || fe_type.family == SUBDIVISION);
              di.resize(di.size() + nc, DofObject::invalid_id);
            }
        }
    }
}



template <typename ScalarDofsFunctor, typename FieldDofsFunctor>
void DofMap::dof_indices (const Elem * const elem,
                          std::vector<dof_id_type> & di,
                          const unsigned int vn,
                          ScalarDofsFunctor scalar_dofs_functor,
                          FieldDofsFunctor field_dofs_functor,
                          int p_level) const
{
  // We now allow elem==nullptr to request just SCALAR dofs
  // libmesh_assert(elem);

  LOG_SCOPE("dof_indices()", "DofMap");

  // Clear the DOF indices vector
  di.clear();

  // Use the default p refinement level?
  if (p_level == -12345)
    p_level = elem ? elem->p_level() : 0;

  const unsigned int vg = this->_variable_group_numbers[vn];
  const VariableGroup & var = this->variable_group(vg);
  const unsigned int vig = vn - var.number();

#ifdef DEBUG
  // Check that sizes match in DEBUG mode
  std::size_t tot_size = 0;
#endif

  if (elem && elem->type() == TRI3SUBDIVISION)
    {
      // Subdivision surface FE require the 1-ring around elem
      const Tri3Subdivision * sd_elem = static_cast<const Tri3Subdivision *>(elem);

      // Ghost subdivision elements have no real dofs
      if (!sd_elem->is_ghost())
        {
          // Determine the nodes contributing to element elem
          std::vector<const Node *> elem_nodes;
          MeshTools::Subdivision::find_one_ring(sd_elem, elem_nodes);

          _dof_indices(*elem, p_level, di, vg, vig, elem_nodes.data(),
                       cast_int<unsigned int>(elem_nodes.size()), vn,
#ifdef DEBUG
                       tot_size,
#endif
                       field_dofs_functor);
        }

      return;
    }

  // Get the dof numbers
  if (var.type().family == SCALAR &&
      (!elem ||
       var.active_on_subdomain(elem->subdomain_id())))
    {
#ifdef DEBUG
      tot_size += var.type().order;
#endif
      std::vector<dof_id_type> di_new;
      this->SCALAR_dof_indices(di_new,vn);
      scalar_dofs_functor(*elem, di, di_new);
    }
  else if (elem)
    _dof_indices(*elem, p_level, di, vg, vig, elem->get_nodes(),
                 elem->n_nodes(), vn,
#ifdef DEBUG
                 tot_size,
#endif
                 field_dofs_functor);

#ifdef DEBUG
  libmesh_assert_equal_to (tot_size, di.size());
#endif
}

inline
StaticCondensationDofMap & DofMap::get_static_condensation()
{
  libmesh_assert(_sc);
  return *_sc;
}

} // namespace libMesh

#endif // LIBMESH_DOF_MAP_H
