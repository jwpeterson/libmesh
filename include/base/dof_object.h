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



#ifndef LIBMESH_DOF_OBJECT_H
#define LIBMESH_DOF_OBJECT_H

// Local includes
#include "libmesh/id_types.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/reference_counted_object.h"

// C++ includes
#include <cstddef>
#include <cstring>
#include <vector>
#include <memory>

namespace libMesh
{

// Forward declarations
class DofObject;


/**
 * The \p DofObject defines an abstract base class for objects that
 * have degrees of freedom associated with them.  Examples of such
 * objects are the \p Node and \p Elem classes.  This class can
 * not be instantiated, only derived from.
 *
 * \author Benjamin S. Kirk
 * \date 2003, 2011
 */

class DofObject : public ReferenceCountedObject<DofObject>
{
#ifdef LIBMESH_IS_UNIT_TESTING
public:
#else
protected:
#endif

  /**
   * Constructor. Protected so that you can't instantiate one of these
   * except as a part of a Node or Elem.
   */
  DofObject ();

public:

  /**
   * Destructor.
   */
  ~DofObject () = default;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * This object on the last mesh.  Useful for projecting
   * solutions from one mesh to another.
   *
   * Public access to old_dof_object is now officially deprecated and will
   * be removed in future libMesh versions.  Use the \p get_old_dof_object()
   * accessor instead.
   */
#ifndef LIBMESH_ENABLE_DEPRECATED
protected:
#endif
  std::unique_ptr<DofObject> old_dof_object;

public:
  /**
   * Pointer accessor for previously public old_dof_object. If you
   * want to assert that the old_dof_object pointer is valid as well,
   * consider using the get_old_dof_object_ref() accessor instead.
   */
  DofObject * get_old_dof_object() { return old_dof_object.get(); }
  const DofObject * get_old_dof_object() const  { return old_dof_object.get(); }

  /**
   * As above, but do not use in situations where the old_dof_object
   * may be nullptr, since this function asserts that the
   * old_dof_object is valid before returning a reference to it.
   */
  DofObject & get_old_dof_object_ref()
  {
    libmesh_assert(old_dof_object);
    return *old_dof_object;
  }

  const DofObject & get_old_dof_object_ref() const
  {
    libmesh_assert(old_dof_object);
    return *old_dof_object;
  }

  /**
   * Sets the \p old_dof_object to nullptr
   */
  void clear_old_dof_object ();

  /**
   * Sets the \p old_dof_object to a copy of \p this
   */
  void set_old_dof_object ();

#endif

  /**
   * Clear the \p DofMap data structures holding degree of freedom
   * data.
   *
   * If any extra integers are associated with this \p DofObject,
   * their count and values are unchanged.
   */
  void clear_dofs ();

  /**
   * Sets all degree of freedom numbers to \p invalid_id
   */
  void invalidate_dofs (const unsigned int sys_num = libMesh::invalid_uint);

  /**
   * Sets the id to \p invalid_id
   */
  void invalidate_id ();

  /**
   * Sets the processor id to \p invalid_processor_id
   */
  void invalidate_processor_id ();

  /**
   * Invalidates all the indices for this \p DofObject
   */
  void invalidate ();

  /**
   * \returns The number of degrees of freedom associated with
   * system \p s directly stored on this object. Optionally only
   * degrees of freedom for variable number \p var are counted.  Does
   * not count degrees of freedom only indirectly associated with this
   * object, such as those stored on an element's nodes when \p
   * n_dofs() is called on the element itself.
   */
  unsigned int n_dofs (const unsigned int s,
                       const unsigned int var =
                       libMesh::invalid_uint) const;

  /**
   * \returns The \p id for this \p DofObject
   */
  dof_id_type id () const;

  /**
   * \returns The \p id for this \p DofObject as a writable reference.
   */
  dof_id_type & set_id ();

  /**
   * \returns The globally \p unique_id for this \p DofObject
   */
  unique_id_type unique_id () const;

#ifdef LIBMESH_ENABLE_DEPRECATED
  /**
   * \returns The globally \p unique_id for this \p DofObject as a
   * writable reference.  Deprecated; use the API taking an input
   * instead.
   */
  unique_id_type & set_unique_id ();
#endif // LIBMESH_ENABLE_DEPRECATED

  /**
   * Sets the \p unique_id for this \p DofObject
   */
  void set_unique_id (unique_id_type new_id);

  /**
   * Sets the \p id for this \p DofObject
   */
  void set_id (const dof_id_type dofid)
  { this->set_id() = dofid; }

  /**
   * \returns \p true if this \p DofObject has a valid \p id set,
   * \p false otherwise.
   */
  bool valid_id () const;

  /**
   * \returns \p true if this \p DofObject has a valid \p unique_id set,
   * \p false otherwise.
   */
  bool valid_unique_id () const;

  /**
   * \returns The processor that this DofObject belongs to.
   *
   * When partitioning and DoF numbering have been performed by
   * libMesh, every current DoF on this DofObject will belong to its
   * processor.
   */
  processor_id_type processor_id () const;

  /**
   * \returns The processor that this DofObject belongs to as a
   * writable reference.
   */
  processor_id_type & processor_id ();

  /**
   * Sets the \p processor_id for this \p DofObject.
   */
  void processor_id (const processor_id_type pid);

  /**
   * \returns \p true if this \p DofObject has a valid \p id set,
   * \p false otherwise.
   */
  bool valid_processor_id () const;

  /**
   * \returns The number of systems associated with this
   * \p DofObject
   */
  unsigned int n_systems() const;

  /**
   * \returns The total number of pseudo-systems associated with this
   * \p DofObject :
   * n_systems(), plus one iff \p this->has_extra_integers()
   */
  unsigned int n_pseudo_systems() const;

  /**
   * Sets the number of systems for this \p DofObject.  If this number
   * is a change, also clears all variable count and DoF indexing
   * associated with this \p DofObject.
   *
   * If any extra integers are associated with this \p DofObject,
   * their count and values are unchanged.
   */
  void set_n_systems (const unsigned int s);

  /**
   * Sets the value on this object of the extra integer associated
   * with \p index, which should have been obtained via a call to \p
   * MeshBase::add_elem_integer or \p MeshBase::add_node_integer
   */
  void set_extra_integer (const unsigned int index, const dof_id_type value);

  /**
   * Gets the value on this object of the extra integer associated
   * with \p index, which should have been obtained via a call to \p
   * MeshBase::add_elem_integer or \p MeshBase::add_node_integer
   */
  dof_id_type get_extra_integer (const unsigned int index) const;

  /**
   * Sets the value on this object of the extra datum associated
   * with \p index, which should have been obtained via a call to \p
   * MeshBase::add_elem_datum or \p MeshBase::add_node_datum using
   * the same type T.
   */
  template <typename T>
  void set_extra_datum (const unsigned int index, const T value);

  /**
   * Gets the value on this object of the extra datum associated
   * with \p index, which should have been obtained via a call to \p
   * MeshBase::add_elem_datum or \p MeshBase::add_node_datum using
   * the same type T.
   */
  template <typename T>
  T get_extra_datum (const unsigned int index) const;


  /**
   * Adds an additional system to the \p DofObject
   */
  void add_system ();

  /**
   * \returns The number of \p VariableGroup variable groups
   * associated with system \p s for this \p DofObject
   */
  unsigned int n_var_groups(const unsigned int s) const;

  /**
   * \returns The number of \p Variable variables associated
   * with \p VariableGroup \p vg in system \p s for this \p DofObject
   */
  unsigned int n_vars(const unsigned int s,
                      const unsigned int vg) const;

  /**
   * \returns The number of \p Variable variables associated
   * with system \p s for this \p DofObject
   */
  unsigned int n_vars(const unsigned int s) const;

  /**
   * Sets number of variables in each group associated with system \p s for this
   * \p DofObject. Implicit in this is also setting the number of \p VariableGroup
   * variable groups for the system.
   * Has the effect of setting the number of components
   * to 0 even when called even with (nvg == this->n_var_groups(s)).
   */
  void set_n_vars_per_group(const unsigned int s,
                            const std::vector<unsigned int> & nvpg);

  /**
   * \returns The number of components for variable \p var
   * of system \p s associated with this \p DofObject.
   * For example, the \p HIERARCHIC shape functions may
   * have multiple DoFs associated with one node.  Another
   * example is the \p MONOMIALs, where only the elements
   * hold the DoFs.  For the different spatial directions,
   * and orders, see \p FE.
   */
  unsigned int n_comp(const unsigned int s,
                      const unsigned int var) const;

  /**
   * \returns The number of components for \p VariableGroup \p vg
   * of system \p s associated with this \p DofObject.
   * For example, the \p HIERARCHIC shape functions may
   * have multiple DoFs associated with one node.  Another
   * example is the \p MONOMIALs, where only the elements
   * hold the DoFs.  For the different spatial directions,
   * and orders, see \p FE.
   */
  unsigned int n_comp_group(const unsigned int s,
                            const unsigned int vg) const;

  /**
   * Sets the number of components for \p Variable \p var
   * of system \p s associated with this \p DofObject
   */
  void set_n_comp(const unsigned int s,
                  const unsigned int var,
                  const unsigned int ncomp);

  /**
   * Sets the number of components for \p VariableGroup \p vg
   * of system \p s associated with this \p DofObject
   */
  void set_n_comp_group(const unsigned int s,
                        const unsigned int vg,
                        const unsigned int ncomp);

  /**
   * \returns The global degree of freedom number for variable \p var,
   * component \p comp for system \p s associated with this \p DofObject
   *
   * When partitioning and DoF numbering have been performed by
   * libMesh, every current DoF on this DofObject will belong to its
   * processor.
   */
  dof_id_type dof_number(const unsigned int s,
                         const unsigned int var,
                         const unsigned int comp) const;

  /**
   * \returns The global degree of freedom number for variable group
   * \p vg, variable index \p vig within the group, component \p comp
   * out of \p n_comp, for system \p s on this \p DofObject
   *
   * Even users who need to call dof_number from user code probably
   * don't want to call this overload.
   */
  dof_id_type dof_number(const unsigned int s,
                         const unsigned int vg,
                         const unsigned int vig,
                         const unsigned int comp,
                         const unsigned int n_comp) const;

  /**
   * \returns A pair consisting of the variable group number and the
   * offset index from the start of that group for variable \p var on
   * system \p s associated with this \p DofObject
   */
  std::pair<unsigned int, unsigned int>
  var_to_vg_and_offset(const unsigned int s,
                       const unsigned int var) const;

  /**
   * Sets the global degree of freedom number for variable \p var,
   * component \p comp for system \p s associated with this \p DofObject
   */
  void set_dof_number(const unsigned int s,
                      const unsigned int var,
                      const unsigned int comp,
                      const dof_id_type dn);

  /**
   * \returns \p true if any system has variables which have been assigned,
   * \p false otherwise.
   */
  bool has_dofs(const unsigned int s=libMesh::invalid_uint) const;

  /**
   * \p VariableGroup DoF indices are indexed as
   * id = base + var_in_vg*ncomp + comp
   * This method allows for direct access to the base.
   */
  void set_vg_dof_base(const unsigned int s,
                       const unsigned int vg,
                       const dof_id_type db);

  /**
   * \p VariableGroup DoF indices are indexed as
   * id = base + var_in_vg*ncomp + comp
   * This method allows for direct access to the base.
   */
  dof_id_type vg_dof_base(const unsigned int s,
                          const unsigned int vg) const;

  /**
   * Assigns a set of extra integers to this \p DofObject.  There will
   * now be \p n_integers associated; this *replaces*, not augments,
   * any previous count.
   *
   * Any newly-added values will initially be DofObject::invalid_id
   *
   * If non-integer data is in the set, each datum of type T should be
   * counted sizeof(T)/sizeof(dof_id_type) times in \p n_integers.
   */
  void add_extra_integers (const unsigned int n_integers);

  /**
   * Assigns a set of extra integers to this \p DofObject.  There will
   * now be \p n_integers associated; this *replaces*, not augments,
   * any previous count.
   *
   * Any newly-added values will be copied from \p default_values.
   *
   * If non-integer data is in the set, each datum of type T should be
   * counted sizeof(T)/sizeof(dof_id_type) times in \p n_integers, and
   * its data should be expressed in \p default_values as per memcpy.
   */
  void add_extra_integers (const unsigned int n_integers,
                           const std::vector<dof_id_type> & default_values);

  /**
   * Returns how many extra integers are associated to the \p DofObject
   *
   * If non-integer data has been associated, each datum of type T
   * counts for sizeof(T)/sizeof(dof_id_type) times in the return
   * value.
   */
  unsigned int n_extra_integers () const;

  /**
   * Returns whether extra integers are associated to the \p DofObject
   */
  bool has_extra_integers () const;

  /**
   * An invalid \p id to distinguish an uninitialized \p DofObject
   */
  static const dof_id_type invalid_id = static_cast<dof_id_type>(-1);

  /**
   * An invalid \p unique_id to distinguish an uninitialized \p DofObject
   */
  static const unique_id_type invalid_unique_id = static_cast<unique_id_type>(-1);

  /**
   * An invalid \p processor_id to distinguish DoFs that have
   * not been assigned to a processor.
   */
  static const processor_id_type invalid_processor_id = static_cast<processor_id_type>(-1);

  /**
   * If we pack our indices into an buffer for communications, how
   * many ints do we need?
   */
  unsigned int packed_indexing_size() const;

  /**
   * If we have indices packed into an buffer for communications, how
   * much of that buffer applies to this dof object?
   */
  static unsigned int unpackable_indexing_size
  (std::vector<largest_id_type>::const_iterator begin);

  /**
   * A method for creating our index buffer from packed data -
   * basically with our current implementation we investigate the size
   * term and then copy.
   */
  void unpack_indexing(std::vector<largest_id_type>::const_iterator begin);

  /**
   * A method for creating packed data from our index buffer -
   * basically a copy with prepended size with our current
   * implementation.
   */
  void pack_indexing(std::back_insert_iterator<std::vector<largest_id_type>> target) const;

  /**
   * Print our buffer for debugging.
   */
  void debug_buffer () const;

  /**
   * Print out info for debugging.
   */
  void print_dof_info() const;

  // Deep copy (or almost-copy) of DofObjects is solely for a couple
  // tricky internal uses.
private:

  /**
   * "Copy"-constructor.  Does not copy old_dof_object, but leaves it
   * null in the new object.
   */
  DofObject (const DofObject &);

  /**
   * Convenient factory function that calls either the (deep) copy
   * constructor or the default constructor depending on the input
   * arg. Like the copy constructor, this function is also private. We
   * can't use std::make_unique to construct a DofObject since the
   * copy constructor is private, but we can at least encapsulate the
   * code which calls "new" directly.
   */
  std::unique_ptr<DofObject>
  construct(const DofObject * other = nullptr);

  /**
   * Deep-copying assignment operator
   */
  DofObject & operator= (const DofObject & dof_obj);

  /**
   * Utility function - for variable \p var in system \p s, figure out what
   * variable group it lives in.
   */
  unsigned int var_to_vg (const unsigned int s,
                          const unsigned int var) const;

  /**
   * Utility function - for variable \p var in system \p s, figure out what
   * variable group it lives in.
   */
  unsigned int system_var_to_vg_var (const unsigned int s,
                                     const unsigned int vg,
                                     const unsigned int var) const;

  /**
   * A globally unique id, guaranteed not to change as the mesh is repartitioned or adapted
   */
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  unique_id_type _unique_id;
#endif

  /**
   * The \p id of the \p DofObject
   */
  dof_id_type _id;

  /**
   * The \p processor_id of the \p DofObject.
   * Degrees of freedom are wholly owned by processors,
   * however they may be duplicated on other processors.
   *
   * This is stored as an unsigned short int since we cannot
   * expect to be solving on 65000+ processors any time soon,
   * can we??
   */
  processor_id_type _processor_id;

  /**
   * DoF index information.  This is packed into a contiguous buffer of the following format:
   *
   * \verbatim
   * [hdr end_0 end_1 ... end_{nps-2} (ncv_0 idx_0 ncv_1 idx_1 ... ncv_nv idx_nv)_0
   *                                  (ncv_0 idx_0 ncv_1 idx_1 ... ncv_nv idx_nv)_1
   *                                   ...
   *                                  (ncv_0 idx_0 ncv_1 idx_1 ... ncv_nv idx_nv)_{nps-2} ]
   * \endverbatim
   *
   * 'hdr' determines whether this \p DofObject \p has_extra_integers()
   * associated with it; iff so then it is negative.
   *
   * The total number of "pseudo systems" is nps := abs(hdr).
   *
   * The total number of true systems is
   * \verbatim
   * ns = hdr,            hdr >= 0
   *    = abs(hdr) - 1,   otherwise.
   * \endverbatim
   *
   * 'end_s' is the index past the end of the variable group (or
   * integer) storage for (pseudo) system \p s.
   *
   * \note We specifically do not store the end for the last (pseudo)
   * system - this always _idx_buf.size().
   *
   * As a first example, consider the case of 4 systems, with 3, 0, 1,
   * 2 variable groups, respectively.  The _idx_buf then looks like:
   *
   * \verbatim
   * [4 10 10 12 () (ncv_0 idx_0 ncv_1 idx_1 ncv_2 idx_2) () (ncv_0 idx_0) (ncv_0 idx_0 ncv_1 idx_1)]
   * [0  1  2  3         4     5     6     7     8     9         10    11      12    13    14    15]
   * \endverbatim
   *
   * The ending index for each (pseudo) system is then given by:
   *
   * \verbatim
   * end_s = _idx_buf.size(),                        s == (nps-1),
   *       = _idx_buf[s+1] + has_extra_integers(),   otherwise.
   * \endverbatim
   *
   * The starting indices are not specifically stored, but rather inferred as follows:
   *
   * start_s = abs(_idx_buf[s])
   *
   * Now, the defining characteristic of the \p VariableGroup is that it supports
   * an arbitrary number of variables of the same type.  At the \p DofObject level, what
   * that means is that each \p Variable in the \p VariableGroup will have the same number
   * of nonzero components, and they can all be indexed from the same base number.  We use this
   * information in the ncv_# and idx_# entries as follows:
   *
   * ncv_# = n_vars*ncv_magic + n_comp      for variable group #
   * idx_# = base_offset                    for variable group #
   *
   * the DoF index for a particular component c of variable v within that group is then given by
   *
   * idx_var = idx_# + n_comp*v + c
   *
   * \note There is a subtlety here - "variable v within that group" usually means nothing to the
   * user. This class is either indexed with variable group numbers, or variable numbers counted
   * *within the system*. So for a system with 2 variable groups, 4 and 8 variables each,
   * the 5th variable in the system is the 1st variable in 2nd variable group.
   * (Now of course 0-base everything...  but you get the idea.)
   *
   * When hdr is *negative* when cast to a signed type, then we
   * interpret that to mean there exists one pseudo-system following
   * the true systems, one for which the _idx_buf data stores the
   * values associated with add_extra_integer entries, not ncv and idx
   * data associated with system variables.  We still return only the
   * number of true systems for n_systems(), but we report
   * has_extra_integers() as true iff hdr is negative, and abs(hdr)
   * will reflect the total number of pseudo-systems, n_systems()+1.
   *
   * E.g. if we had added two extra integers to the example case
   * above, the _idx_buf then looks like:
   *
   * \verbatim
   * [-5 11 11 13 17 () (ncv_0 idx_0 ncv_1 idx_1 ncv_2 idx_2) () (ncv_0 idx_0) (ncv_0 idx_0 ncv_1 idx_1) (xtra1 xtra2)]
   * [0   1  2  3  4         5     6     7     8     9    10         11    12      13    14    15    16      17    18]
   * \endverbatim
   */
  typedef dof_id_type index_t;
  typedef std::vector<index_t> index_buffer_t;
  index_buffer_t _idx_buf;

  /**
   * Above we introduced the chimera ncv, which is a hybrid of the form
   * ncv = ncv_magic*nv + nc
   * where nv are the number of identical variables of a given type,
   * and nc is the number of components for this set of variables.
   *
   * It is hoped that by setting this to a power of two, an optimizing compiler
   * will recognize later that  #/ncv_magic is simply a bitshift
   */
  static const index_t ncv_magic = 256; // = 2^8, in case we want to manually bitshift
  static const index_t ncv_magic_exp = 8; // Let's manually bitshift

  /**
   * The starting index for system \p s.
   */
  unsigned int start_idx(const unsigned int s) const;

  /**
   * The ending index for system \p s.
   */
  unsigned int end_idx(const unsigned int s) const;

  /**
   * The starting index for an extra_integers pseudosystem
   */
  unsigned int start_idx_ints() const;

  /**
   * The ending index for an extra_integers pseudosystem
   */
  unsigned int end_idx_ints() const;

  // methods only available for unit testing
#ifdef LIBMESH_IS_UNIT_TESTING
public:
  void set_buffer (const std::vector<dof_id_type> & buf)
  { _idx_buf = buf; }
#endif
};



//------------------------------------------------------
// Inline functions
inline
DofObject::DofObject () :
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _unique_id (invalid_unique_id),
#endif
  _id (invalid_id),
  _processor_id (invalid_processor_id)
{
  this->invalidate();
}



inline
std::unique_ptr<DofObject>
DofObject::construct(const DofObject * other)
{
  return other
    ? std::unique_ptr<DofObject>(new DofObject(*other))
    : std::unique_ptr<DofObject>(new DofObject());
}



inline
void DofObject::invalidate_dofs (const unsigned int sys_num)
{
  const unsigned int n_sys = this->n_systems();
  // If the user does not specify the system number...
  if (sys_num >= n_sys)
    {
      for (auto s : make_range(n_sys))
        for (auto vg : make_range(this->n_var_groups(s)))
          if (this->n_comp_group(s,vg))
            this->set_vg_dof_base(s,vg,invalid_id);
    }
  // ...otherwise invalidate the dofs for all systems
  else
    for (auto vg : make_range(this->n_var_groups(sys_num)))
      if (this->n_comp_group(sys_num,vg))
        this->set_vg_dof_base(sys_num,vg,invalid_id);
}



inline
void DofObject::invalidate_id ()
{
  this->set_id (invalid_id);
}



inline
void DofObject::invalidate_processor_id ()
{
  this->processor_id (invalid_processor_id);
}



inline
void DofObject::invalidate ()
{
  this->invalidate_dofs ();
  this->invalidate_id ();
  this->invalidate_processor_id ();
}



inline
void DofObject::clear_dofs ()
{
  this->set_n_systems(0);
}



inline
unsigned int DofObject::n_dofs (const unsigned int s,
                                const unsigned int var) const
{
  libmesh_assert_less (s, this->n_systems());

  unsigned int num = 0;

  // Count all variables
  if (var == libMesh::invalid_uint)
    for (auto v : make_range(this->n_vars(s)))
      num += this->n_comp(s,v);

  // Only count specified variable
  else
    num = this->n_comp(s,var);

  return num;
}



inline
dof_id_type DofObject::id () const
{
  return _id;
}



inline
dof_id_type & DofObject::set_id ()
{
  return _id;
}



inline
unique_id_type DofObject::unique_id () const
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  libmesh_assert (this->valid_unique_id());
  return _unique_id;
#else
  return invalid_unique_id;
#endif
}



#ifdef LIBMESH_ENABLE_DEPRECATED
inline
unique_id_type & DofObject::set_unique_id ()
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  libmesh_deprecated();
  return _unique_id;
#else
  libmesh_not_implemented();
#endif
}
#endif // LIBMESH_ENABLE_DEPRECATED



inline
void DofObject::set_unique_id (unique_id_type new_id)
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _unique_id = new_id;
#else
  libmesh_ignore(new_id);
  libmesh_not_implemented();
#endif
}



inline
bool DofObject::valid_id () const
{
  return (DofObject::invalid_id != _id);
}



inline
bool DofObject::valid_unique_id () const
{
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  return (DofObject::invalid_unique_id != _unique_id);
#else
  return false;
#endif
}



inline
processor_id_type DofObject::processor_id () const
{
  return _processor_id;
}



inline
processor_id_type & DofObject::processor_id ()
{
  return _processor_id;
}



inline
void DofObject::processor_id (const processor_id_type pid)
{
  this->processor_id() = pid;
}



inline
bool DofObject::valid_processor_id () const
{
  return (DofObject::invalid_processor_id != _processor_id);
}



inline
unsigned int DofObject::n_systems () const
{
  const int hdr = _idx_buf.empty() ?
    0 : cast_int<int>(dof_id_signed_type(_idx_buf[0]));
  return hdr >= 0 ? hdr : (-hdr-1);
}



inline
unsigned int DofObject::n_pseudo_systems () const
{
  const int hdr = _idx_buf.empty() ?
    0 : cast_int<int>(dof_id_signed_type(_idx_buf[0]));
  return std::abs(hdr);
}



inline
unsigned int DofObject::n_var_groups(const unsigned int s) const
{
  libmesh_assert_less (s, this->n_systems());

  return (this->end_idx(s) - this->start_idx(s)) / 2;
}



inline
unsigned int DofObject::n_vars(const unsigned int s,
                               const unsigned int vg) const
{
  libmesh_assert_less (s,  this->n_systems());
  libmesh_assert_less (vg, this->n_var_groups(s));

  const unsigned int start_idx_sys = this->start_idx(s);

  libmesh_assert_less ((start_idx_sys + 2*vg), _idx_buf.size());

  return (cast_int<unsigned int>
          (_idx_buf[start_idx_sys + 2*vg]) >> ncv_magic_exp);
}



inline
unsigned int DofObject::n_vars(const unsigned int s) const
{
  libmesh_assert_less (s, this->n_systems());

  const unsigned int nvg = this->n_var_groups(s);

  unsigned int val=0;

  for (unsigned int vg=0; vg<nvg; vg++)
    val += this->n_vars(s,vg);

  return val;
}




inline
unsigned int DofObject::n_comp(const unsigned int s,
                               const unsigned int var) const
{
  libmesh_assert_less (s,   this->n_systems());
  libmesh_assert_less (var, this->n_vars(s));

  return this->n_comp_group(s,this->var_to_vg(s,var));
}




inline
unsigned int DofObject::n_comp_group(const unsigned int s,
                                     const unsigned int vg) const
{
  libmesh_assert_less (s,  this->n_systems());
  libmesh_assert_less (vg, this->n_var_groups(s));

  const unsigned int
    start_idx_sys = this->start_idx(s);

  libmesh_assert_less ((start_idx_sys + 2*vg), _idx_buf.size());

  return (_idx_buf[start_idx_sys + 2*vg] % ncv_magic);
}



inline
dof_id_type DofObject::dof_number(const unsigned int s,
                                  const unsigned int var,
                                  const unsigned int comp) const
{
  libmesh_assert_less (s,    this->n_systems());
  libmesh_assert_less (var,  this->n_vars(s));
  libmesh_assert_less (comp, this->n_comp(s,var));

  const std::pair<unsigned int, unsigned int>
    vg_vig = this->var_to_vg_and_offset(s,var);

  const unsigned int
    n_comp = this->n_comp_group(s,vg_vig.first);

  return this->dof_number(s, vg_vig.first, vg_vig.second,
                          comp, n_comp);
}



inline
dof_id_type DofObject::dof_number(const unsigned int s,
                                  const unsigned int vg,
                                  const unsigned int vig,
                                  const unsigned int comp,
                                  const unsigned int n_comp) const
{
  libmesh_assert_less (s,   this->n_systems());
  libmesh_assert_less (vg,  this->n_var_groups(s));
  libmesh_assert_less (vig, this->n_vars(s,vg));

  const unsigned int
    start_idx_sys = this->start_idx(s);

  libmesh_assert_less ((start_idx_sys + 2*vg + 1), _idx_buf.size());

  const dof_id_type
    base_idx = _idx_buf[start_idx_sys + 2*vg + 1];

  // if the first component is invalid, they
  // are all invalid
  if (base_idx == invalid_id)
    return invalid_id;

  // otherwise the index is the first component
  // index augmented by the component number
  else
    return cast_int<dof_id_type>(base_idx + vig*n_comp + comp);
}



inline
void
DofObject::set_extra_integer(const unsigned int index,
                             const dof_id_type value)
{
  libmesh_assert_less(index, this->n_extra_integers());
  libmesh_assert_less(this->n_pseudo_systems(), _idx_buf.size());

  const unsigned int start_idx_i = this->start_idx_ints();

  libmesh_assert_less(start_idx_i+index, _idx_buf.size());
  _idx_buf[start_idx_i+index] = value;
}



inline
dof_id_type
DofObject::get_extra_integer (const unsigned int index) const
{
  libmesh_assert_less(index, this->n_extra_integers());
  libmesh_assert_less(this->n_systems(), _idx_buf.size());

  const unsigned int start_idx_i = this->start_idx_ints();

  libmesh_assert_less(start_idx_i+index, _idx_buf.size());
  return _idx_buf[start_idx_i+index];
}



// If we're using a type T that's a class with no trivial
// copy-assignment, -Wclass-memaccess will scream about doing it with
// memcpy, even if (as with boost::multiprecision::float128) this is a
// false positive.
#include "libmesh/ignore_warnings.h"



template <typename T>
inline
void
DofObject::set_extra_datum(const unsigned int index,
                           const T value)
{
#ifndef NDEBUG
  const unsigned int n_more_integers = (sizeof(T)-1)/sizeof(dof_id_type);
#endif
  libmesh_assert_less(index+n_more_integers, this->n_extra_integers());
  libmesh_assert_less(this->n_pseudo_systems(), _idx_buf.size());

  const unsigned int start_idx_i = this->start_idx_ints();

  libmesh_assert_less(start_idx_i+index+n_more_integers, _idx_buf.size());
  std::memcpy(&_idx_buf[start_idx_i+index], &value, sizeof(T));
}



template <typename T>
inline
T
DofObject::get_extra_datum (const unsigned int index) const
{
#ifndef NDEBUG
  const unsigned int n_more_integers = (sizeof(T)-1)/sizeof(dof_id_type);
#endif
  libmesh_assert_less(index+n_more_integers, this->n_extra_integers());
  libmesh_assert_less(this->n_systems(), _idx_buf.size());

  const unsigned int start_idx_i = this->start_idx_ints();

  libmesh_assert_less(start_idx_i+index+n_more_integers, _idx_buf.size());
  T returnval;
  std::memcpy(&returnval, &_idx_buf[start_idx_i+index], sizeof(T));
  return returnval;
}



#include "libmesh/restore_warnings.h"



inline
unsigned int
DofObject::n_extra_integers () const
{
  if (_idx_buf.empty())
    return 0;

  const int hdr = dof_id_signed_type(_idx_buf[0]);
  if (hdr >= 0)
    return 0;

  const unsigned int start_idx_i = this->start_idx_ints();

  return _idx_buf.size() - start_idx_i;
}



inline
bool
DofObject::has_extra_integers () const
{
  if (_idx_buf.empty())
    return 0;

  return (dof_id_signed_type(_idx_buf[0]) < 0);
}



inline
std::pair<unsigned int, unsigned int>
DofObject::var_to_vg_and_offset(const unsigned int s,
                                const unsigned int var) const
{
  std::pair<unsigned int, unsigned int> returnval(0,0);

  unsigned int & vg = returnval.first;
  unsigned int & offset = returnval.second;

  unsigned int vg_start = 0;
  for (; ; vg++)
    {
      libmesh_assert_less(vg, this->n_var_groups(s));

      const unsigned int vg_end = vg_start + this->n_vars(s,vg);
      if (var < vg_end)
        {
          offset = var - vg_start;
          return returnval;
        }
      vg_start = vg_end;
    }
}



inline
bool DofObject::has_dofs (const unsigned int sys) const
{
  if (sys == libMesh::invalid_uint)
    {
      for (auto s : make_range(this->n_systems()))
        if (this->n_vars(s))
          return true;
    }

  else
    {
      libmesh_assert_less (sys, this->n_systems());

      if (this->n_vars(sys))
        return true;
    }

  return false;
}



inline
unsigned int DofObject::start_idx (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_systems());
  libmesh_assert_less (s, _idx_buf.size());

  return cast_int<unsigned int>(std::abs(dof_id_signed_type(_idx_buf[s])));
}



inline
unsigned int DofObject::end_idx (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_systems());
  libmesh_assert_less (s, _idx_buf.size());

  return ((s+1) == this->n_pseudo_systems()) ?
    cast_int<unsigned int>(_idx_buf.size()) :
    cast_int<unsigned int>(_idx_buf[s+1]);
}



inline
unsigned int DofObject::start_idx_ints () const
{
  libmesh_assert (this->has_extra_integers());

  unsigned int n_sys = this->n_systems();

  libmesh_assert_less(this->n_systems(), _idx_buf.size());
  return n_sys ? cast_int<unsigned int>(_idx_buf[this->n_systems()]) :
                 (n_sys+1);
}



inline
unsigned int DofObject::end_idx_ints () const
{
  libmesh_assert (this->has_extra_integers());

  return cast_int<unsigned int>(_idx_buf.size());
}



inline
void DofObject::set_vg_dof_base(const unsigned int s,
                                const unsigned int vg,
                                const dof_id_type db)
{
  libmesh_assert_less (s,  this->n_systems());
  libmesh_assert_less (vg, this->n_var_groups(s));

  const unsigned int
    start_idx_sys = this->start_idx(s);

  libmesh_assert_less ((start_idx_sys + 2*vg + 1), _idx_buf.size());

  _idx_buf[start_idx_sys + 2*vg + 1] = db;

  libmesh_assert_equal_to (this->vg_dof_base(s,vg), db);
}



inline
dof_id_type DofObject::vg_dof_base(const unsigned int s,
                                   const unsigned int vg) const
{
  libmesh_assert_less (s,  this->n_systems());
  libmesh_assert_less (vg, this->n_var_groups(s));

  const unsigned int
    start_idx_sys = this->start_idx(s);

  libmesh_assert_less ((start_idx_sys + 2*vg + 1), _idx_buf.size());

  // #ifdef DEBUG
  //   std::cout << " [ ";
  //   for (auto i : _idx_buf)
  //     std::cout << i << " ";
  //   std::cout << "]\n";
  // #endif

  return _idx_buf[start_idx_sys + 2*vg + 1];
}



inline
unsigned int DofObject::var_to_vg (const unsigned int s,
                                   const unsigned int var) const
{
  const unsigned int
    nvg = this->n_var_groups(s);

  for (unsigned int vg=0, vg_end=0; vg<nvg; vg++)
    {
      vg_end += this->n_vars(s,vg);
      if (var < vg_end) return vg;
    }

  libmesh_error_msg("Error: could not map variable " << var << " to variable group.");
}



inline
unsigned int DofObject::system_var_to_vg_var (const unsigned int s,
                                              const unsigned int vg,
                                              const unsigned int var) const
{
  unsigned int accumulated_sum=0;

  for (unsigned int vgc=0; vgc<vg; vgc++)
    accumulated_sum += this->n_vars(s,vgc);

  libmesh_assert_less_equal (accumulated_sum, var);

  return (var - accumulated_sum);
}

/**
 * Comparison object to use with DofObject pointers.  This sorts by id(),
 * so when we iterate over a set of DofObjects we visit the objects in
 * order of increasing ID.
 */
struct CompareDofObjectsByID
{
  bool operator()(const DofObject * a,
                  const DofObject * b) const
  {
    libmesh_assert (a);
    libmesh_assert (b);

    return a->id() < b->id();
  }
};

struct CompareDofObjectsByPIDAndThenID
{
  bool operator()(const DofObject * a,
                  const DofObject * b) const
  {
    libmesh_assert (a);
    libmesh_assert (b);

    if (a->processor_id() < b->processor_id())
      return true;
    if (b->processor_id() < a->processor_id())
      return false;

    return a->id() < b->id();
  }
};

} // namespace libMesh


#endif // #ifndef LIBMESH_DOF_OBJECT_H
