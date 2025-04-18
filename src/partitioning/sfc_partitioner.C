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
#include "libmesh/libmesh_config.h"
#include "libmesh/elem.h"
#include "libmesh/enum_partitioner_type.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/sfc_partitioner.h"

#ifdef LIBMESH_HAVE_SFCURVES
namespace Sfc {
extern "C" {
#     include "sfcurves.h"
}
}
#else
#  include "libmesh/linear_partitioner.h"
#endif

namespace libMesh
{


PartitionerType SFCPartitioner::type() const
{
  return SFC_PARTITIONER;
}


void SFCPartitioner::partition_range(MeshBase & mesh,
                                     MeshBase::element_iterator beg,
                                     MeshBase::element_iterator end,
                                     unsigned int n)
{
  // Check for easy returns
  if (beg == end)
    return;

  if (n == 1)
    {
      this->single_partition_range (beg, end);
      return;
    }

  libmesh_assert_greater (n, 0);

  // What to do if the sfcurves library IS NOT present
#ifndef LIBMESH_HAVE_SFCURVES

  libmesh_do_once(
    libMesh::out << "ERROR: The library has been built without"    << std::endl
                 << "Space Filling Curve support.  Using a linear" << std::endl
                 << "partitioner instead!" << std::endl;);

  LinearPartitioner lp;
  lp.partition_range (mesh, beg, end, n);

  // What to do if the sfcurves library IS present
#else

  LOG_SCOPE("partition_range()", "SFCPartitioner");

  // We don't yet support distributed meshes with this Partitioner
  if (!mesh.is_serial())
    libmesh_not_implemented();

  const dof_id_type n_range_elem = std::distance(beg, end);
  const dof_id_type max_elem_id = mesh.max_elem_id();

  // The forward_map maps the range's element ids into a contiguous
  // block of indices.
  std::vector<dof_id_type> forward_map (max_elem_id, DofObject::invalid_id);

  // the reverse_map maps the contiguous ids back
  // to active elements
  std::vector<Elem *> reverse_map (n_range_elem, nullptr);

  std::vector<double> x      (n_range_elem);
  std::vector<double> y      (n_range_elem);
  std::vector<double> z      (n_range_elem);
  std::vector<int>    table  (n_range_elem);

  // Map the range's element ids into a contiguous range.
  dof_id_type el_num = 0;

  for (auto & elem : as_range(beg, end))
    {
      libmesh_assert_less (elem->id(), forward_map.size());
      libmesh_assert_less (el_num, reverse_map.size());

      forward_map[elem->id()] = el_num;
      reverse_map[el_num] = elem;
      el_num++;
    }
  libmesh_assert_equal_to (el_num, n_range_elem);

  // Get the vertex average for each range element.
  for (const auto & elem : as_range(beg, end))
    {
      libmesh_assert_less (elem->id(), forward_map.size());

      const Point p = elem->vertex_average();

      x[forward_map[elem->id()]] = double(p(0));
      y[forward_map[elem->id()]] = double(p(1));
      z[forward_map[elem->id()]] = double(p(2));
    }

  // We need an integer reference to pass to the Sfc interface.
  int size = static_cast<int>(n_range_elem);

  // build the space-filling curve
  if (_sfc_type == "Hilbert")
    Sfc::hilbert (x.data(),
                  y.data(),
                  z.data(),
                  &size,
                  table.data());

  else if (_sfc_type == "Morton")
    Sfc::morton  (x.data(),
                  y.data(),
                  z.data(),
                  &size,
                  table.data());

  else
    {
      libMesh::out << "ERROR: Unknown type: " << _sfc_type << std::endl
                   << " Valid types are"                   << std::endl
                   << "  \"Hilbert\""                      << std::endl
                   << "  \"Morton\""                       << std::endl
                   << " "                                  << std::endl
                   << "Partitioning with a Hilbert curve." << std::endl;

      Sfc::hilbert (x.data(),
                    y.data(),
                    z.data(),
                    &size,
                    table.data());
    }


  // Assign the partitioning to the range elements
  {
    // {
    //   std::ofstream out ("sfc.dat");
    //   out << "variables=x,y,z" << std::endl;
    //   out << "zone f=point" << std::endl;
    //   for (unsigned int i=0; i<n_range_elem; i++)
    //     out << x[i] << " " << y[i] << " " << z[i] << std::endl;
    // }

    const dof_id_type blksize = (n_range_elem + n - 1) / n;

    for (dof_id_type i=0; i<n_range_elem; i++)
      {
        libmesh_assert_less (table[i] - 1, reverse_map.size());

        Elem * elem = reverse_map[table[i] - 1];

        elem->processor_id() = cast_int<processor_id_type>(i/blksize);
      }
  }

#endif
}



void SFCPartitioner::_do_partition (MeshBase & mesh,
                                    const unsigned int n)
{
  this->partition_range(mesh,
                        mesh.active_elements_begin(),
                        mesh.active_elements_end(),
                        n);
}

} // namespace libMesh
