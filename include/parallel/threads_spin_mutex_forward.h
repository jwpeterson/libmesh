// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_THREADS_SPIN_MUTEX_FORWARD_H
#define LIBMESH_THREADS_SPIN_MUTEX_FORWARD_H

#include "libmesh/libmesh_config.h"

// Lightweight forward declaration of Threads::spin_mutex for headers that
// need the type without including the full threads.h.
//
// In the TBB backend spin_mutex is a typedef for tbb::spin_mutex, which is
// itself introduced into namespace tbb via an inline-namespace using-alias
// (not a direct class definition).  A plain "class spin_mutex" forward
// declaration in namespace tbb would therefore conflict with that alias, so
// we must include tbb/spin_mutex.h to make the type complete and then
// re-declare the typedef.  For the pthreads and serial backends spin_mutex
// is a proper class and a forward declaration suffices.

#ifdef LIBMESH_HAVE_TBB_API

#  include "libmesh/ignore_warnings.h"
#  include "tbb/spin_mutex.h"
#  include "libmesh/restore_warnings.h"

namespace libMesh
{
namespace Threads
{
typedef tbb::spin_mutex spin_mutex;
} // namespace Threads
} // namespace libMesh

#else

namespace libMesh
{
namespace Threads
{
class spin_mutex;
} // namespace Threads
} // namespace libMesh

#endif // LIBMESH_HAVE_TBB_API

#endif // LIBMESH_THREADS_SPIN_MUTEX_FORWARD_H
