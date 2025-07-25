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



#ifndef LIBMESH_WRAPPED_FUNCTION_H
#define LIBMESH_WRAPPED_FUNCTION_H

// Local Includes
#include "libmesh/dense_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/function_base.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"
#include "libmesh/system.h"
#include "libmesh/fe_interface.h"

// C++ includes
#include <cstddef>
#include <memory>

namespace libMesh
{

/**
 * \brief Wrap a libMesh-style function pointer into a FunctionBase object.
 *
 * This class provides a wrapper with which to evaluate a
 * (libMesh-style) function pointer in a FunctionBase-compatible
 * interface. All overridden virtual functions are documented in
 * function_base.h.
 *
 * \author Roy Stogner
 * \date 2012
 * \note To wrap an ordinary function pointer use the AnalyticFunction class.
 */
template <typename Output=Number>
class WrappedFunction : public FunctionBase<Output>
{
public:

  /**
   * Constructor to wrap scalar-valued function pointers.
   */
  WrappedFunction (const System & sys,
                   Output fptr(const Point & p,
                               const Parameters & parameters,
                               const std::string & sys_name,
                               const std::string & unknown_name) = nullptr,
                   const Parameters * parameters = nullptr,
                   unsigned int varnum=0)
    : _sys(sys),
      _fptr(fptr),
      _parameters(parameters),
      _varnum(varnum)
  {
    this->_initialized = true;
    if (!parameters)
      _parameters = &sys.get_equation_systems().parameters;
  }

  /**
   * The move/copy ctor and destructor are defaulted for this class.
   */
  WrappedFunction (WrappedFunction &&) = default;
  WrappedFunction (const WrappedFunction &) = default;
  virtual ~WrappedFunction () = default;

  /**
   * This class contains a const reference so it can't be assigned.
   */
  WrappedFunction & operator= (const WrappedFunction &) = delete;
  WrappedFunction & operator= (WrappedFunction &&) = delete;

  virtual std::unique_ptr<FunctionBase<Output>> clone () const override;

  virtual Output operator() (const Point & p,
                             const Real time = 0.) override;

  virtual void operator() (const Point & p,
                           const Real time,
                           DenseVector<Output> & output) override;

  virtual Output component (unsigned int i,
                            const Point & p,
                            Real time=0.) override;

protected:

  const System & _sys;

  Output (*_fptr)(const Point & p,
                  const Parameters & parameters,
                  const std::string & sys_name,
                  const std::string & unknown_name);

  const Parameters * _parameters;

  unsigned int _varnum;
};


// ------------------------------------------------------------
// WrappedFunction inline methods


template <typename Output>
inline
Output WrappedFunction<Output>::operator() (const Point & p,
                                            const Real /*time*/)
{
  libmesh_assert(_fptr);
  libmesh_assert(_parameters);
  return _fptr(p,
               *_parameters,
               _sys.name(),
               _sys.variable_name(_varnum));
}


template <typename Output>
inline
std::unique_ptr<FunctionBase<Output>>
WrappedFunction<Output>::clone () const
{
  return std::make_unique<WrappedFunction<Output>>
    (_sys, _fptr, _parameters, _varnum);
}


template <typename Output>
inline
void WrappedFunction<Output>::operator() (const Point & p,
                                          const Real /*time*/,
                                          DenseVector<Output> & output)
{
  libmesh_assert(_fptr);
  libmesh_assert(_parameters);

  // We fill each entry of output with a single scalar component of
  // the data in our System
  libmesh_assert_equal_to (output.size(), _sys.n_components());

  // Loop over variables, then over each component in
  // vector-valued variables, evaluating each.
  const unsigned int n_vars = _sys.n_vars();
  for (unsigned int v = 0; v != n_vars; ++v)
    {
      const auto n_components =
        _sys.variable(v).n_components(_sys.get_mesh());
      if (n_components == 1)
        output(_sys.variable_scalar_number(v,0)) =
          _fptr(p, *_parameters, _sys.name(), _sys.variable_name(v));
      else
        {
          // Right now our only non-scalar variable type is the
          // SCALAR variables.  The irony is priceless.
          libmesh_assert_equal_to (_sys.variable(v).type().family, SCALAR);

          // We pass the point (j,0,0) to an old-style fptr function
          // pointer to distinguish the different scalars within the
          // SCALAR variable.
          for (unsigned int j=0; j != n_components; ++j)
            output(_sys.variable_scalar_number(v,j)) =
              _fptr(Point(j,0,0), *_parameters,
                    _sys.name(), _sys.variable_name(v));
        }
    }
}


template <typename Output>
inline
Output WrappedFunction<Output>::component (unsigned int i,
                                           const Point & p,
                                           Real /*time*/)
{
  libmesh_assert(_fptr);
  libmesh_assert(_parameters);

  // Loop over variables, then over each component in
  // FEFamily SCALAR variables.
  unsigned int vc = 0;
  const unsigned int n_vars = _sys.n_vars();
  for (unsigned int v = 0; v != n_vars; ++v)
    {
      const auto & var_fe_type = _sys.variable_type(v);
      const auto n_components = _sys.variable(v).n_components(_sys.get_mesh());
      if (i >= vc + n_components)
      {
        vc += n_components;
        continue;
      }

      if (n_components > 1 && var_fe_type.family != SCALAR)
        libmesh_error_msg(
            "WrappedFunction::component cannot currently evaluate vector finite element families");

      if (n_components == 1)
        return _fptr(p, *_parameters, _sys.name(), _sys.variable_name(v));
      else
        {
          libmesh_assert_equal_to (_sys.variable(i).type().family, SCALAR);

          // We pass the point (j,0,0) to an old-style fptr function
          // pointer to distinguish the different scalars within the
          // SCALAR variable.
          for (unsigned int j=0; j != n_components; ++j)
            if (i == vc + j)
              return _fptr(Point(j,0,0), *_parameters,
                           _sys.name(), _sys.variable_name(v));
        }
    }

  libmesh_error_msg("Component index " << i << " not found in system " << _sys.name());
  return Output();
}



} // namespace libMesh

#endif // LIBMESH_WRAPPED_FUNCTION_H
