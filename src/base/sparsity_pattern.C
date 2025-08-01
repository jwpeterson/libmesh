// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local includes
#include "libmesh/sparsity_pattern.h"

// libMesh includes
#include "libmesh/coupling_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/ghosting_functor.h"
#include "libmesh/hashword.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/utility.h"
#include "libmesh/static_condensation_dof_map.h"

// TIMPI includes
#include "timpi/communicator.h"


namespace libMesh
{
namespace SparsityPattern
{

//-------------------------------------------------------
// we need to implement these constructors here so that
// a full DofMap definition is available.
Build::Build (const DofMap & dof_map_in,
              const CouplingMatrix * dof_coupling_in,
              const std::vector<GhostingFunctor *> & coupling_functors_in,
              const bool implicit_neighbor_dofs_in,
              const bool need_full_sparsity_pattern_in,
              const bool calculate_constrained_in,
              const StaticCondensationDofMap * const sc_in) :
  ParallelObject(dof_map_in),
  dof_map(dof_map_in),
  dof_coupling(dof_coupling_in),
  coupling_functors(coupling_functors_in),
  implicit_neighbor_dofs(implicit_neighbor_dofs_in),
  need_full_sparsity_pattern(need_full_sparsity_pattern_in),
  calculate_constrained(calculate_constrained_in),
  sc(sc_in),
  sparsity_pattern(),
  nonlocal_pattern(),
  n_nz(),
  n_oz()
{}



Build::Build (Build & other, Threads::split) :
  ParallelObject(other),
  dof_map(other.dof_map),
  dof_coupling(other.dof_coupling),
  coupling_functors(other.coupling_functors),
  implicit_neighbor_dofs(other.implicit_neighbor_dofs),
  need_full_sparsity_pattern(other.need_full_sparsity_pattern),
  calculate_constrained(other.calculate_constrained),
  sc(other.sc),
  hashed_dof_sets(other.hashed_dof_sets),
  sparsity_pattern(),
  nonlocal_pattern(),
  n_nz(),
  n_oz()
{}



#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)

void _dummy_function(void) {}

#endif



void Build::sorted_connected_dofs(const Elem * elem,
                                  std::vector<dof_id_type> & dofs_vi,
                                  unsigned int vi)
{
  if (this->sc)
  {
    // We build a sparsity pattern that will match the size of the condensed system. This is so that
    // we have the data necessary to init the reduced system matrix
    dofs_vi.clear();

    auto total_and_uncondensed_from_scalar_dofs_functor =
        [&dofs_vi](const Elem & /*elem*/,
                   std::vector<dof_id_type> & dof_indices,
                   const std::vector<dof_id_type> & scalar_dof_indices)
    {
      dof_indices.insert(dof_indices.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
      dofs_vi.insert(dofs_vi.end(), scalar_dof_indices.begin(), scalar_dof_indices.end());
    };

    auto total_and_uncondensed_from_field_dofs_functor =
        [&dofs_vi, this](const Elem & functor_elem,
                        const unsigned int node_num,
                        const unsigned int var_num,
                        std::vector<dof_id_type> & dof_indices,
                        const dof_id_type field_dof)
    {
      dof_indices.push_back(field_dof);
      if (this->sc->uncondensed_vars().count(var_num) ||
          (node_num != invalid_uint && !functor_elem.is_internal(node_num)))
        dofs_vi.push_back(field_dof);
    };

    dof_map.dof_indices(elem,
                        dummy_vec,
                        vi,
                        total_and_uncondensed_from_scalar_dofs_functor,
                        total_and_uncondensed_from_field_dofs_functor,
                        elem->p_level());
  }
  else
    dof_map.dof_indices (elem, dofs_vi, vi);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  dof_map.find_connected_dofs (dofs_vi);
#endif
  // We can be more efficient if we sort the element DOFs into
  // increasing order
  std::sort(dofs_vi.begin(), dofs_vi.end());

  // Handle cases where duplicate nodes are intentionally assigned to
  // a single element.
  dofs_vi.erase(std::unique(dofs_vi.begin(), dofs_vi.end()), dofs_vi.end());
}



void Build::handle_vi_vj(const std::vector<dof_id_type> & element_dofs_i,
                         const std::vector<dof_id_type> & element_dofs_j)
{
  const unsigned int n_dofs_on_element_i =
    cast_int<unsigned int>(element_dofs_i.size());

  const processor_id_type proc_id     = dof_map.processor_id();
  const dof_id_type first_dof_on_proc = dof_map.first_dof(proc_id);
  const dof_id_type end_dof_on_proc   = dof_map.end_dof(proc_id);

  std::vector<dof_id_type>
    dofs_to_add;

  const unsigned int n_dofs_on_element_j =
    cast_int<unsigned int>(element_dofs_j.size());

  // It only makes sense to compute hashes and see if we can skip
  // doing work when there are a "large" amount of DOFs for a given
  // element. The cutoff for "large" is somewhat arbitrarily chosen
  // based on a test case with a spider node that resulted in O(10^3)
  // entries in element_dofs_i for O(10^3) elements. Making this
  // number larger will disable the hashing optimization in more
  // cases.
  bool dofs_seen = false;
  if (n_dofs_on_element_j > 0 && n_dofs_on_element_i > 256)
    {
      auto hash_i = Utility::hashword(element_dofs_i);
      auto hash_j = Utility::hashword(element_dofs_j);
      auto final_hash = Utility::hashword2(hash_i, hash_j);
      auto result = hashed_dof_sets.insert(final_hash);
      // if insert failed, we have already seen these dofs
      dofs_seen = !result.second;
    }

  // there might be 0 dofs for the other variable on the same element
  // (when subdomain variables do not overlap) and that's when we do
  // not do anything
  if (n_dofs_on_element_j > 0 && !dofs_seen)
    {
      for (unsigned int i=0; i<n_dofs_on_element_i; i++)
        {
          const dof_id_type ig = element_dofs_i[i];

          SparsityPattern::Row * row;

          // We save non-local row components for now so we can
          // communicate them to other processors later.

          if ((ig >= first_dof_on_proc) &&
              (ig <  end_dof_on_proc))
            {
              // This is what I mean
              // libmesh_assert_greater_equal ((ig - first_dof_on_proc), 0);
              // but do the test like this because ig and
              // first_dof_on_proc are unsigned ints
              libmesh_assert_greater_equal (ig, first_dof_on_proc);
              libmesh_assert_less (ig, (sparsity_pattern.size() +
                                        first_dof_on_proc));

              row = &sparsity_pattern[ig - first_dof_on_proc];
            }
          else
            {
              row = &nonlocal_pattern[ig];
            }

          // If the row is empty we will add *all*
          // the element j DOFs, so just do that.
          if (row->empty())
            {
              row->insert(row->end(),
                          element_dofs_j.begin(),
                          element_dofs_j.end());
            }
          else
            {
              // Build a list of the DOF indices not found in the
              // sparsity pattern
              dofs_to_add.clear();

              // Cache iterators.  Low will move forward, subsequent
              // searches will be on smaller ranges
              SparsityPattern::Row::iterator
                low  = std::lower_bound
                (row->begin(), row->end(), element_dofs_j.front()),
                high = std::upper_bound
                (low,          row->end(), element_dofs_j.back());

              for (unsigned int j=0; j<n_dofs_on_element_j; j++)
                {
                  const dof_id_type jg = element_dofs_j[j];

                  // See if jg is in the sorted range
                  std::pair<SparsityPattern::Row::iterator,
                            SparsityPattern::Row::iterator>
                    pos = std::equal_range (low, high, jg);

                  // Must add jg if it wasn't found
                  if (pos.first == pos.second)
                    dofs_to_add.push_back(jg);

                  // pos.first is now a valid lower bound for any
                  // remaining element j DOFs. (That's why we sorted them.)
                  // Use it for the next search
                  low = pos.first;
                }

              // Add to the sparsity pattern
              if (!dofs_to_add.empty())
                {
                  const std::size_t old_size = row->size();

                  row->insert (row->end(),
                               dofs_to_add.begin(),
                               dofs_to_add.end());

                  SparsityPattern::sort_row
                    (row->begin(), row->begin()+old_size,
                     row->end());
                }
            }
        } // End dofs-of-var-i loop
    } // End if-dofs-of-var-j
}



void Build::operator()(const ConstElemRange & range)
{
  // Compute the sparsity structure of the global matrix.  This can be
  // fed into a PetscMatrixBase to allocate exactly the number of nonzeros
  // necessary to store the matrix.  This algorithm should be linear
  // in the (# of elements)*(# nodes per element)
  sparsity_pattern.resize(dof_map.n_local_dofs());

  // Handle dof coupling specified by library and user coupling functors
  {
    const unsigned int n_var = dof_map.n_variables();

    std::vector<std::vector<dof_id_type> > element_dofs_i(n_var);

    std::vector<const Elem *> coupled_neighbors;
    for (const auto & elem : range)
      {
        // Make some fake element iterators defining a range
        // pointing to only this element.
        Elem * const * elempp = const_cast<Elem * const *>(&elem);
        Elem * const * elemend = elempp+1;

        const MeshBase::const_element_iterator fake_elem_it =
          MeshBase::const_element_iterator(elempp,
                                           elemend,
                                           Predicates::NotNull<Elem * const *>());

        const MeshBase::const_element_iterator fake_elem_end =
          MeshBase::const_element_iterator(elemend,
                                           elemend,
                                           Predicates::NotNull<Elem * const *>());

        GhostingFunctor::map_type elements_to_couple;
        DofMap::CouplingMatricesSet temporary_coupling_matrices;

        dof_map.merge_ghost_functor_outputs(elements_to_couple,
                                            temporary_coupling_matrices,
                                            dof_map.coupling_functors_begin(),
                                            dof_map.coupling_functors_end(),
                                            fake_elem_it,
                                            fake_elem_end,
                                            DofObject::invalid_processor_id);
        for (unsigned int vi=0; vi<n_var; vi++)
          this->sorted_connected_dofs(elem, element_dofs_i[vi], vi);

        for (unsigned int vi=0; vi<n_var; vi++)
          for (const auto & [partner, ghost_coupling] : elements_to_couple)
            {
              // Loop over coupling matrix row variables if we have a
              // coupling matrix, or all variables if not.
              if (ghost_coupling)
                {
                  libmesh_assert_equal_to (ghost_coupling->size(), n_var);
                  ConstCouplingRow ccr(vi, *ghost_coupling);

                  for (const auto & idx : ccr)
                    {
                      if (partner == elem)
                        this->handle_vi_vj(element_dofs_i[vi], element_dofs_i[idx]);
                      else
                        {
                          std::vector<dof_id_type> partner_dofs;
                          this->sorted_connected_dofs(partner, partner_dofs, idx);
                          this->handle_vi_vj(element_dofs_i[vi], partner_dofs);
                        }
                    }
                }
              else
                {
                  for (unsigned int vj = 0; vj != n_var; ++vj)
                    {
                      if (partner == elem)
                        this->handle_vi_vj(element_dofs_i[vi], element_dofs_i[vj]);
                      else
                        {
                          std::vector<dof_id_type> partner_dofs;
                          this->sorted_connected_dofs(partner, partner_dofs, vj);
                          this->handle_vi_vj(element_dofs_i[vi], partner_dofs);
                        }
                    }
                }
            } // End ghosted element loop
      } // End range element loop
  } // End ghosting functor section
}



void Build::join (const SparsityPattern::Build & other)
{
  libmesh_assert_equal_to (sparsity_pattern.size(), other.sparsity_pattern.size());

  for (dof_id_type r=0; r<dof_map.n_local_dofs(); r++)
    {
      // increment the number of on and off-processor nonzeros in this row
      // (note this will be an upper bound unless we need the full sparsity pattern)
      SparsityPattern::Row       & my_row    = sparsity_pattern[r];
      const SparsityPattern::Row & their_row = other.sparsity_pattern[r];

      // simple copy if I have no dofs
      if (my_row.empty())
        my_row = their_row;

      // otherwise add their DOFs to mine, resort, and re-unique the row
      else if (!their_row.empty()) // do nothing for the trivial case where
        {                          // their row is empty
          my_row.insert (my_row.end(),
                         their_row.begin(),
                         their_row.end());

          // We cannot use SparsityPattern::sort_row() here because it expects
          // the [begin,middle) [middle,end) to be non-overlapping.  This is not
          // necessarily the case here, so use std::sort()
          std::sort (my_row.begin(), my_row.end());

          my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
        }
    }

  // Move nonlocal row information to ourselves; the other thread
  // won't need it in the map after that.
  for (const auto & p : other.nonlocal_pattern)
    {
#ifndef NDEBUG
      const dof_id_type dof_id = p.first;

      processor_id_type dbg_proc_id = 0;
      while (dof_id >= dof_map.end_dof(dbg_proc_id))
        dbg_proc_id++;
      libmesh_assert (dbg_proc_id != this->processor_id());
#endif

      const SparsityPattern::Row & their_row = p.second;

      // We should have no empty values in a map
      libmesh_assert (!their_row.empty());

      if (auto my_it = nonlocal_pattern.find(p.first);
          my_it == nonlocal_pattern.end())
        {
          //          nonlocal_pattern[it->first].swap(their_row);
          nonlocal_pattern[p.first] = their_row;
        }
      else
        {
          SparsityPattern::Row & my_row = my_it->second;

          my_row.insert (my_row.end(),
                         their_row.begin(),
                         their_row.end());

          // We cannot use SparsityPattern::sort_row() here because it expects
          // the [begin,middle) [middle,end) to be non-overlapping.  This is not
          // necessarily the case here, so use std::sort()
          std::sort (my_row.begin(), my_row.end());

          my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
        }
    }

  // Combine the other thread's hashed_dof_sets with ours.
  hashed_dof_sets.insert(other.hashed_dof_sets.begin(),
                         other.hashed_dof_sets.end());
}



void Build::parallel_sync ()
{
  parallel_object_only();
  libmesh_assert(this->comm().verify(need_full_sparsity_pattern));

  const auto n_dofs_on_proc  = dof_map.n_local_dofs();
  const auto local_first_dof = dof_map.first_dof();

  // The data to send
  std::map<processor_id_type, std::vector<dof_id_type>> ids_to_send;
  std::map<processor_id_type, std::vector<Row>> rows_to_send;

  // Loop over the nonlocal rows and transform them into the new datastructure
  NonlocalGraph::iterator it = nonlocal_pattern.begin();
  while (it != nonlocal_pattern.end())
  {
    const auto dof_id = it->first;
    auto & row = it->second;

    processor_id_type proc_id = 0;
    while (dof_id >= dof_map.end_dof(proc_id))
      proc_id++;

    ids_to_send[proc_id].push_back(dof_id);

    // Note this invalidates the data in nonlocal_pattern
    rows_to_send[proc_id].push_back(std::move(row));

    // Might as well remove it since it's invalidated anyway
    it = nonlocal_pattern.erase(it);
  }

  std::map<processor_id_type, std::vector<dof_id_type>> received_ids_map;

  auto ids_action_functor =
    [& received_ids_map]
    (processor_id_type pid,
     const std::vector<dof_id_type> & received_ids)
    {
      received_ids_map.emplace(pid, received_ids);
    };

  Parallel::push_parallel_vector_data(this->comm(), ids_to_send,
                                      ids_action_functor);

  auto rows_action_functor =
    [this,
     & received_ids_map,
     local_first_dof]
    (processor_id_type pid,
     const std::vector<Row> & received_rows)
    {
      const std::vector<dof_id_type> & received_ids = libmesh_map_find(received_ids_map, pid);

      std::size_t n_rows = received_rows.size();
      libmesh_assert_equal_to(n_rows, received_ids.size());

      for (auto i : IntRange<std::size_t>(0, n_rows))
        {
          const auto r = received_ids[i];
          libmesh_assert(dof_map.local_index(r));

          const auto my_r = r - local_first_dof;

          auto & their_row = received_rows[i];

          auto & my_row = sparsity_pattern[my_r];

          // They wouldn't have sent an empty row
          libmesh_assert(!their_row.empty());

          // We can end up with an empty row on a dof that touches our
          // inactive elements but not our active ones
          if (my_row.empty())
            {
              my_row.assign (their_row.begin(), their_row.end());
            }
          else
            {
              my_row.insert (my_row.end(),
                             their_row.begin(),
                             their_row.end());

              // We cannot use SparsityPattern::sort_row() here because it expects
              // the [begin,middle) [middle,end) to be non-overlapping.  This is not
              // necessarily the case here, so use std::sort()
              std::sort (my_row.begin(), my_row.end());

              my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
            }

        }
    };

  Parallel::push_parallel_vector_data(this->comm(), rows_to_send,
                                      rows_action_functor);

  // We should have sent everything at this point.
  libmesh_assert (nonlocal_pattern.empty());

  // assert these are empty because std::vector::resize will only append the specified element value
  // if the new size is greater than the current size. Elements whose indices are less than the
  // current size are untouched
  libmesh_assert(n_nz.empty());
  libmesh_assert(n_oz.empty());
  n_nz.resize (n_dofs_on_proc, 0);
  n_oz.resize (n_dofs_on_proc, 0);

  const dof_id_type first_dof_on_proc = dof_map.first_dof();
  const dof_id_type end_dof_on_proc   = dof_map.end_dof();

  for (dof_id_type i=0; i<n_dofs_on_proc; i++)
    {
      // Get the row of the sparsity pattern
      SparsityPattern::Row & row = sparsity_pattern[i];

      for (const auto & df : row)
        if ((df < first_dof_on_proc) || (df >= end_dof_on_proc))
          n_oz[i]++;
        else
          n_nz[i]++;

      libmesh_assert(n_nz[i] <= n_dofs_on_proc);

      // If we're not building a full sparsity pattern, then we want
      // to avoid overcounting these entries as much as possible.
      if (!need_full_sparsity_pattern)
        row.clear();
    }
}


void Build::apply_extra_sparsity_object(SparsityPattern::AugmentSparsityPattern & asp)
{
  asp.augment_sparsity_pattern (sparsity_pattern, n_nz, n_oz);
}


std::size_t Build::n_nonzeros() const
{
  // At some point I'll remember that "C++17" compilers don't always
  // come with complete C++17 standard libraries.
  // std::size_t total_nonzeros = std::reduce(n_nz.begin(), n_nz.end(), std::size_t(0));
  // total_nonzeros += std::reduce(n_oz.begin(), n_oz.end(), std::size_t(0));

  std::size_t total_nonzeros = 0;
  for (auto nnzi : n_nz)
    total_nonzeros += nnzi;
  for (auto nozi : n_oz)
    total_nonzeros += nozi;

  this->comm().sum(total_nonzeros);
  return total_nonzeros;
}


} // namespace SparsityPattern
} // namespace libMesh
