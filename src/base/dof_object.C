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



// Local includes
#include "libmesh/dof_object.h"

namespace libMesh
{

// ------------------------------------------------------------
// DofObject class static member -now initialized in header
const dof_id_type       DofObject::invalid_id;
const unique_id_type    DofObject::invalid_unique_id;
const processor_id_type DofObject::invalid_processor_id;



// Copy Constructor
DofObject::DofObject (const DofObject & dof_obj) :
  ReferenceCountedObject<DofObject>(),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _unique_id     (dof_obj._unique_id),
#endif
  _id            (dof_obj._id),
  _processor_id  (dof_obj._processor_id),
  _idx_buf       (dof_obj._idx_buf)
{
  // DO NOT copy old_dof_object, because this isn't a *real* copy
  // constructor, it's a "copy almost everything" constructor that
  // is intended to be used solely for internal construction of
  // old_dof_object, never for a true deep copy where the newly
  // created object really matches the source object.

  // Check that everything worked
#ifdef DEBUG

  libmesh_assert_equal_to (this->n_systems(), dof_obj.n_systems());

  for (auto s : make_range(this->n_systems()))
    {
      libmesh_assert_equal_to (this->n_vars(s),       dof_obj.n_vars(s));
      libmesh_assert_equal_to (this->n_var_groups(s), dof_obj.n_var_groups(s));

      for (auto vg : make_range(this->n_var_groups(s)))
        libmesh_assert_equal_to (this->n_vars(s,vg), dof_obj.n_vars(s,vg));

      for (auto v : make_range(this->n_vars(s)))
        {
          libmesh_assert_equal_to (this->n_comp(s,v), dof_obj.n_comp(s,v));

          for (auto c : make_range(this->n_comp(s,v)))
            libmesh_assert_equal_to (this->dof_number(s,v,c), dof_obj.dof_number(s,v,c));
        }
    }

#endif
}


// Deep-copying assignment operator
DofObject & DofObject::operator= (const DofObject & dof_obj)
{
  if (&dof_obj == this)
    return *this;

#ifdef LIBMESH_ENABLE_AMR
  this->clear_old_dof_object();
  this->old_dof_object = this->construct(dof_obj.get_old_dof_object());
#endif

  _id           = dof_obj._id;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _unique_id    = dof_obj._unique_id;
#endif
  _processor_id = dof_obj._processor_id;
  _idx_buf      = dof_obj._idx_buf;


  // Check that everything worked
#ifdef DEBUG

  libmesh_assert_equal_to (this->n_systems(), dof_obj.n_systems());

  for (auto s : make_range(this->n_systems()))
    {
      libmesh_assert_equal_to (this->n_vars(s),       dof_obj.n_vars(s));
      libmesh_assert_equal_to (this->n_var_groups(s), dof_obj.n_var_groups(s));

      for (auto vg : make_range(this->n_var_groups(s)))
        libmesh_assert_equal_to (this->n_vars(s,vg), dof_obj.n_vars(s,vg));

      for (auto v : make_range(this->n_vars(s)))
        {
          libmesh_assert_equal_to (this->n_comp(s,v), dof_obj.n_comp(s,v));

          for (auto c : make_range(this->n_comp(s,v)))
            libmesh_assert_equal_to (this->dof_number(s,v,c), dof_obj.dof_number(s,v,c));
        }
    }

#endif

  return *this;
}





#ifdef LIBMESH_ENABLE_AMR

void  DofObject::clear_old_dof_object ()
{
  this->old_dof_object.reset(nullptr);
}



void DofObject::set_old_dof_object ()
{
  this->clear_old_dof_object();

  libmesh_assert (!this->old_dof_object);

  // Make a new DofObject, assign a copy of \p this.
  // Make sure the copy ctor for DofObject works!!
  this->old_dof_object = this->construct(this);
}

#endif



void DofObject::set_n_systems (const unsigned int ns)
{
  const unsigned int old_ns = this->n_systems();

  // Check for trivial return
  if (ns == old_ns)
    return;

  const unsigned int nei = this->n_extra_integers();
  const dof_id_type header_size = ns + bool(nei);
  const dof_id_type hdr = nei ?
    static_cast<dof_id_type>(-static_cast<std::ptrdiff_t>(header_size))
    : header_size;
  index_buffer_t new_buf(header_size + nei, hdr);
  if (nei)
    {
      const unsigned int start_idx_ints = old_ns ?
        cast_int<unsigned int>(_idx_buf[old_ns]) :
        1;
      libmesh_assert_less(start_idx_ints, _idx_buf.size());
      std::copy(_idx_buf.begin()+start_idx_ints,
                _idx_buf.end(),
                new_buf.begin()+header_size);
      if (ns)
        std::fill(new_buf.begin()+1, new_buf.begin()+ns+1, ns+1);
    }

  // vector swap trick to force deallocation when shrinking
  new_buf.swap(_idx_buf);

#ifdef DEBUG
  libmesh_assert_equal_to(nei, this->n_extra_integers());

  // check that all systems now exist and that they have 0 size
  libmesh_assert_equal_to (ns, this->n_systems());
  for (auto s : make_range(this->n_systems()))
    {
      libmesh_assert_equal_to (this->n_vars(s),       0);
      libmesh_assert_equal_to (this->n_var_groups(s), 0);
    }
#endif
}



void DofObject::add_system()
{
  // quick return?
  if (this->n_systems() == 0)
    {
      this->set_n_systems(1);
      return;
    }

  // cache this value before we screw it up!
  const unsigned int ns_orig = this->n_systems();

  DofObject::index_buffer_t::iterator it = _idx_buf.begin() + ns_orig;

  // Create the entry for the new system indicating 0 variables.
  //
  // increment the number of systems and the offsets for each of
  // the systems including the new one we just added.
  if (this->has_extra_integers())
    {
      // this inserts the extra_integers' start position as the start
      // position for the new system.  We'll increment all those
      // counts in one sweep next, to account for header expansion.
      _idx_buf.insert(it, *it);

      _idx_buf[0]--;
      for (unsigned int i=1; i<ns_orig+2; i++)
        {
          libmesh_assert_less(i, _idx_buf.size());
          _idx_buf[i]++;
        }
    }
  else
    {
      // this inserts the current vector size at the position for the
      // new system
      _idx_buf.insert(it, cast_int<dof_id_type>(_idx_buf.size()));

      for (unsigned int i=0; i<ns_orig+1; i++)
        {
          libmesh_assert_less(i, _idx_buf.size());
          _idx_buf[i]++;
        }
    }

  libmesh_assert_equal_to (this->n_systems(), (ns_orig+1));
  libmesh_assert_equal_to (this->n_vars(ns_orig), 0);
  libmesh_assert_equal_to (this->n_var_groups(ns_orig), 0);
}



void DofObject::set_n_vars_per_group(const unsigned int s,
                                     const std::vector<unsigned int> & nvpg)
{
  const unsigned int n_sys = this->n_systems();

  libmesh_assert_less (s, n_sys);

  // number of variable groups for this system - inferred
  const unsigned int nvg = cast_int<unsigned int>(nvpg.size());

  // BSK - note that for compatibility with the previous implementation
  // calling this method when (nvars == this->n_vars()) requires that
  // we invalidate the DOF indices and set the number of components to 0.
  // Note this was a bit of a surprise to me - there was no quick return in
  // the old method, which caused removal and readdition of the DOF indices
  // even in the case of (nvars == this->n_vars()), resulting in n_comp(s,v)
  // implicitly becoming 0 regardless of any previous value.
  // quick return?
  if (nvg == this->n_var_groups(s))
    {
      for (unsigned int vg=0; vg<nvg; vg++)
        {
          this->set_n_comp_group(s,vg,0);
          libmesh_assert_equal_to (this->n_vars(s,vg), nvpg[vg]);
        }
      return;
    }

  const bool hei = this->has_extra_integers();

  // since there is ample opportunity to screw up other systems, let us
  // cache their current sizes and later assert that they are unchanged.
#ifdef DEBUG
  const unsigned int nei = this->n_extra_integers();

  DofObject::index_buffer_t old_system_sizes, old_extra_integers;
  old_system_sizes.reserve(n_sys);
  old_extra_integers.reserve(nei);

  for (unsigned int s_ctr=0; s_ctr<n_sys; s_ctr++)
    old_system_sizes.push_back(this->n_var_groups(s_ctr));

  for (unsigned int ei=0; ei != nei; ++ei)
    old_extra_integers.push_back(this->get_extra_integer(ei));
#endif

  // remove current indices if we have some
  if (this->n_var_groups(s) != 0)
    {
      const unsigned int old_nvg_s = this->n_var_groups(s);

      DofObject::index_buffer_t::iterator
        it  = _idx_buf.begin(),
        end = _idx_buf.begin();

      std::advance(it,  this->start_idx(s));
      std::advance(end, this->end_idx(s));
      _idx_buf.erase(it,end);

      for (unsigned int ctr=(s+1); ctr<n_sys; ctr++)
        _idx_buf[ctr] -= 2*old_nvg_s;

      if (hei)
        _idx_buf[n_sys] -= 2*old_nvg_s;
    }

  // better not have any now!
  libmesh_assert_equal_to (this->n_var_groups(s), 0);

  // Make sure we didn't screw up any of our sizes!
#ifdef DEBUG
  for (auto s_ctr : make_range(this->n_systems()))
    if (s_ctr != s)
      libmesh_assert_equal_to (this->n_var_groups(s_ctr), old_system_sizes[s_ctr]);

  libmesh_assert_equal_to (nei, this->n_extra_integers());

  for (unsigned int ei=0; ei != nei; ++ei)
    libmesh_assert_equal_to(old_extra_integers[ei], this->get_extra_integer(ei));
#endif

  // OK, if the user requested 0 that is what we have
  if (nvg == 0)
    return;

  {
    // array to hold new indices
    DofObject::index_buffer_t var_idxs(2*nvg);
    for (unsigned int vg=0; vg<nvg; vg++)
      {
        var_idxs[2*vg    ] = ncv_magic*nvpg[vg] + 0;
        var_idxs[2*vg + 1] = invalid_id - 1;
      }

    DofObject::index_buffer_t::iterator it = _idx_buf.begin();
    std::advance(it, this->end_idx(s));
    _idx_buf.insert(it, var_idxs.begin(), var_idxs.end());

    for (unsigned int ctr=(s+1); ctr<n_sys; ctr++)
      _idx_buf[ctr] += 2*nvg;

    if (hei)
      _idx_buf[n_sys] += 2*nvg;

    // resize _idx_buf to fit so no memory is wasted.
    DofObject::index_buffer_t(_idx_buf).swap(_idx_buf);
  }

  libmesh_assert_equal_to (nvg, this->n_var_groups(s));

#ifdef DEBUG

  libmesh_assert_equal_to (this->n_var_groups(s), nvpg.size());

  for (auto vg : make_range(this->n_var_groups(s)))
    {
      libmesh_assert_equal_to (this->n_vars(s,vg), nvpg[vg]);
      libmesh_assert_equal_to (this->n_comp_group(s,vg), 0);
    }

  for (auto v : make_range(this->n_vars(s)))
    libmesh_assert_equal_to (this->n_comp(s,v), 0);

  // again, all other system sizes should be unchanged!
  for (auto s_ctr : make_range(this->n_systems()))
    if (s_ctr != s)
      libmesh_assert_equal_to (this->n_var_groups(s_ctr), old_system_sizes[s_ctr]);

  // Extra integers count and values should also be unchanged!
  libmesh_assert_equal_to (nei, this->n_extra_integers());

  for (unsigned int ei=0; ei != nei; ++ei)
    libmesh_assert_equal_to(old_extra_integers[ei], this->get_extra_integer(ei));
#endif
}



void DofObject::set_n_comp(const unsigned int s,
                           const unsigned int var,
                           const unsigned int ncomp)
{
  libmesh_assert_less (s,   this->n_systems());
  libmesh_assert_less (var, this->n_vars(s));

  this->set_n_comp_group(s, this->var_to_vg(s,var), ncomp);
}



void DofObject::set_n_comp_group(const unsigned int s,
                                 const unsigned int vg,
                                 const unsigned int ncomp)
{
  libmesh_assert_less (s,  this->n_systems());
  libmesh_assert_less (vg, this->n_var_groups(s));

  // Check for trivial return
  if (ncomp == this->n_comp_group(s,vg)) return;

#ifndef NDEBUG
  if (ncomp >= ncv_magic)
    {
      const index_t ncvm = ncv_magic;
      libmesh_error_msg("ERROR: ncomp must be less than DofObject::ncv_magic!\n" \
                        << "ncomp = "                                   \
                        << ncomp                                \
                        << ", ncv_magic = "                     \
                        << ncvm                                 \
                        << "\nrecompile and try again!");
    }
#endif

  const unsigned int
    start_idx_sys = this->start_idx(s),
    n_vars_group  = this->n_vars(s,vg),
    base_offset   = start_idx_sys + 2*vg;

  libmesh_assert_less ((base_offset + 1), _idx_buf.size());

  // if (ncomp)
  //   libMesh::out << "s,vg,ncomp="
  //       << s  << ","
  //       << vg << ","
  //       << ncomp << '\n';

  // set the number of components, maintaining the number
  // of variables in the group
  _idx_buf[base_offset] = ncv_magic*n_vars_group + ncomp;

  // We use (invalid_id - 1) to signify no
  // components for this object
  _idx_buf[base_offset + 1] = (ncomp == 0) ? invalid_id - 1 : invalid_id;

  // this->debug_buffer();
  // libMesh::out << "s,vg = " << s << "," << vg << '\n'
  //     << "base_offset=" << base_offset << '\n'
  //     << "this->n_comp(s,vg)=" << this->n_comp(s,vg) << '\n'
  //     << "this->n_comp_group(s,vg)=" << this->n_comp_group(s,vg) << '\n'
  //     << "this->n_vars(s,vg)=" << this->n_vars(s,vg) << '\n'
  //     << "this->n_var_groups(s)=" << this->n_var_groups(s) << '\n';

  libmesh_assert_equal_to (ncomp, this->n_comp_group(s,vg));
}



void DofObject::set_dof_number(const unsigned int s,
                               const unsigned int var,
                               const unsigned int comp,
                               const dof_id_type dn)
{
  libmesh_assert_less (s,    this->n_systems());
  libmesh_assert_less (var,  this->n_vars(s));
  libmesh_assert_less (comp, this->n_comp(s,var));

  const unsigned int
    vg            = this->var_to_vg(s,var),
#ifndef NDEBUG
    ncg           = this->n_comp_group(s,vg),
#endif
    vig           = this->system_var_to_vg_var(s,vg,var),
    start_idx_sys = this->start_idx(s);

  libmesh_assert_less ((start_idx_sys + 2*vg + 1), _idx_buf.size());

  dof_id_type & base_idx = _idx_buf[start_idx_sys + 2*vg + 1];

  // We intend to change all dof numbers together or not at all
  if (comp || vig)
    libmesh_assert ((dn == invalid_id && base_idx == invalid_id) ||
                    (dn == base_idx + vig*ncg + comp));

  // only explicitly store the base index for vig==0, comp==0
  else
    base_idx = dn;

  libmesh_assert_equal_to (this->dof_number(s, var, comp), dn);
}



void
DofObject::add_extra_integers (const unsigned int n_integers)
{
  if (_idx_buf.empty())
    {
      if (n_integers)
        {
          _idx_buf.resize(n_integers+1, DofObject::invalid_id);
          _idx_buf[0] = dof_id_type(-1);
        }
      return;
    }
  else
    {
      const int hdr = dof_id_signed_type(_idx_buf[0]);

      // We already have some extra integers, but may need more or
      // less now.
      if (hdr < 0)
        {
          const unsigned int old_n_integers = this->n_extra_integers();
          if (n_integers != old_n_integers)
            {
              // Make or remove space as needed by count change
              _idx_buf.resize(_idx_buf.size()+n_integers-old_n_integers, DofObject::invalid_id);

              // The start index for the extra integers is unchanged.
            }
        }
      else if (n_integers)
      // We had no extra integers, but need to add some
        {
          // Mark the DofObject as holding extra integers
          _idx_buf[0] = dof_id_type(-hdr-1);

          // Insert the integer start position
          DofObject::index_buffer_t::iterator it = _idx_buf.begin() + hdr;
          _idx_buf.insert(it, _idx_buf.size()+1);

          // Increment the previous system start positions to account
          // for the new header entry creating an offset
          for (int i=1; i<hdr; i++)
            _idx_buf[i]++;

          // Append space for extra integers
          _idx_buf.resize(_idx_buf.size()+n_integers, DofObject::invalid_id);
        }
    }
}



void
DofObject::add_extra_integers (const unsigned int n_integers,
                               const std::vector<dof_id_type> & default_values)
{
  libmesh_assert_equal_to(n_integers, default_values.size());

  const unsigned int n_old_integers = this->n_extra_integers();
  this->add_extra_integers(n_integers);
  if (n_integers > n_old_integers)
    {
      const unsigned int n_more_integers = n_integers - n_old_integers;
      std::copy(default_values.begin()+n_old_integers,
                default_values.end(),
                _idx_buf.end()-n_more_integers);
    }
}



unsigned int DofObject::packed_indexing_size() const
{
  return
    cast_int<unsigned int> (
#ifdef LIBMESH_ENABLE_AMR
                            ((old_dof_object == nullptr) ? 0 : old_dof_object->packed_indexing_size()) + 2 +
#else
                            1 +
#endif
                            _idx_buf.size());
}



unsigned int
DofObject::unpackable_indexing_size(std::vector<largest_id_type>::const_iterator begin)
{
#ifdef LIBMESH_ENABLE_AMR
  const bool has_old_dof_object = cast_int<bool>(*begin++);

  static const int dof_header_size = 2;
#else
  static const bool has_old_dof_object = false;
  static const int dof_header_size = 1;
#endif

  const largest_id_type this_indexing_size = *begin++;

  return cast_int<unsigned int>
    (dof_header_size + this_indexing_size +
     (has_old_dof_object ?
      unpackable_indexing_size(begin+this_indexing_size) : 0));
}


void DofObject::unpack_indexing(std::vector<largest_id_type>::const_iterator begin)
{
  _idx_buf.clear();

#ifdef LIBMESH_ENABLE_AMR
  this->clear_old_dof_object();
  const bool has_old_dof_object = cast_int<bool>(*begin++);
#endif

  const largest_id_type size = *begin++;
  _idx_buf.reserve(size);
  std::copy(begin, begin+size, back_inserter(_idx_buf));

  // Check as best we can for internal consistency now
  libmesh_assert(_idx_buf.empty() ||
                 (std::abs(dof_id_signed_type(_idx_buf[0])) <=
                  cast_int<dof_id_signed_type>(_idx_buf.size())));
#ifdef DEBUG
  if (!_idx_buf.empty())
    {
      const int hdr = cast_int<int>(dof_id_signed_type(_idx_buf[0]));
      const unsigned int ns = hdr >= 0 ? hdr : (-hdr-1);
      for (unsigned int i=1; i < ns; ++i)
        {
          if (hdr > 0 || i > 1)
            libmesh_assert_greater_equal (_idx_buf[i], _idx_buf[i-1]);
          else
            libmesh_assert_greater_equal (_idx_buf[i], ns);
          libmesh_assert_equal_to ((_idx_buf[i] - _idx_buf[i-1])%2, 0);
          libmesh_assert_less_equal (_idx_buf[i], _idx_buf.size());
        }
      if (hdr < 0 && ns > 0)
        libmesh_assert_less_equal(_idx_buf[ns], _idx_buf.size());
    }
#endif

#ifdef LIBMESH_ENABLE_AMR
  if (has_old_dof_object)
    {
      this->old_dof_object = this->construct();
      this->old_dof_object->unpack_indexing(begin+size);
    }
#endif
}


void
DofObject::pack_indexing(std::back_insert_iterator<std::vector<largest_id_type>> target) const
{
#ifdef LIBMESH_ENABLE_AMR
  // We might need to pack old_dof_object too
  *target++ = (old_dof_object == nullptr) ? 0 : 1;
#endif

  *target++ = _idx_buf.size();
  std::copy(_idx_buf.begin(), _idx_buf.end(), target);

#ifdef LIBMESH_ENABLE_AMR
  if (old_dof_object)
    old_dof_object->pack_indexing(target);
#endif
}



void DofObject::debug_buffer () const
{
  libMesh::out << " [ ";
  for (const auto & idx : _idx_buf)
    libMesh::out << idx << " ";
  libMesh::out << "]\n";
}



void DofObject::print_dof_info() const
{
  libMesh::out << this->id() << " [ ";

  for (auto s : make_range(this->n_systems()))
    {
      libMesh::out << "s:" << s << " ";
      for (auto var : make_range(this->n_vars(s)))
        {
          libMesh::out << "v:" << var << " ";
          for (auto comp : make_range(this->n_comp(s,var)))
            {
              libMesh::out << "c:" << comp << " dof:" << this->dof_number(s,var,comp) << " ";
            }
        }
    }

  libMesh::out << "]\n";
}



} // namespace libMesh
