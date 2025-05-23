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
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/exact_solution.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/function_base.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_function.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/quadrature.h"
#include "libmesh/wrapped_function.h"
#include "libmesh/wrapped_functor.h"
#include "libmesh/fe_interface.h"
#include "libmesh/raw_accessor.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/utility.h"

// C++ Includes
#include <memory>


namespace libMesh
{

ExactSolution::ExactSolution(const EquationSystems & es) :
  _equation_systems(es),
  _equation_systems_fine(nullptr),
  _extra_order(1)
{
  // Initialize the _errors data structure which holds all
  // the eventual values of the error.
  for (auto sys : make_range(_equation_systems.n_systems()))
    {
      // Reference to the system
      const System & system = _equation_systems.get_system(sys);

      // The name of the system
      const std::string & sys_name = system.name();

      // The SystemErrorMap to be inserted
      ExactSolution::SystemErrorMap sem;

      for (auto var : make_range(system.n_vars()))
        {
          // The name of this variable
          const std::string & var_name = system.variable_name(var);
          sem[var_name] = std::vector<Real>(5, 0.);
        }

      _errors[sys_name] = sem;
    }
}


ExactSolution::ExactSolution(ExactSolution &&) = default;
ExactSolution::~ExactSolution() = default;


void ExactSolution::
set_excluded_subdomains(const std::set<subdomain_id_type> & excluded)
{
  _excluded_subdomains = excluded;
}

void ExactSolution::attach_reference_solution (const EquationSystems * es_fine)
{
  libmesh_assert(es_fine);
  _equation_systems_fine = es_fine;

  // If we're using a fine grid solution, we're not using exact values
  _exact_values.clear();
  _exact_derivs.clear();
  _exact_hessians.clear();
}


void ExactSolution::attach_exact_value (ValueFunctionPointer fptr)
{
  libmesh_assert(fptr);

  // Clear out any previous _exact_values entries, then add a new
  // entry for each system.
  _exact_values.clear();

  for (auto sys : make_range(_equation_systems.n_systems()))
    {
      const System & system = _equation_systems.get_system(sys);
      _exact_values.emplace_back(std::make_unique<WrappedFunctor<Number>>(system, fptr, &_equation_systems.parameters));
    }

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = nullptr;
}


void ExactSolution::attach_exact_values (const std::vector<FunctionBase<Number> *> & f)
{
  // Automatically delete any previous _exact_values entries, then add a new
  // entry for each system.
  _exact_values.clear();

  for (auto ptr : f)
    _exact_values.emplace_back(ptr ? std::make_unique<WrappedFunctor<Number>>(*ptr) : nullptr);
}


void ExactSolution::attach_exact_values (const std::vector<FEMFunctionBase<Number> *> & f)
{
  // Automatically delete any previous _exact_values entries, then add a new
  // entry for each system.
  _exact_values.clear();

  for (auto ptr : f)
    _exact_values.emplace_back(ptr ? ptr->clone() : nullptr);
}


void ExactSolution::attach_exact_value (unsigned int sys_num,
                                        FunctionBase<Number> * f)
{
  if (_exact_values.size() <= sys_num)
    _exact_values.resize(sys_num+1);

  if (f)
    _exact_values[sys_num] = std::make_unique<WrappedFunctor<Number>>(*f);
}


void ExactSolution::attach_exact_value (unsigned int sys_num,
                                        FEMFunctionBase<Number> * f)
{
  if (_exact_values.size() <= sys_num)
    _exact_values.resize(sys_num+1);

  if (f)
    _exact_values[sys_num] = f->clone();
}


void ExactSolution::attach_exact_deriv (GradientFunctionPointer gptr)
{
  libmesh_assert(gptr);

  // Clear out any previous _exact_derivs entries, then add a new
  // entry for each system.
  _exact_derivs.clear();

  for (auto sys : make_range(_equation_systems.n_systems()))
    {
      const System & system = _equation_systems.get_system(sys);
      _exact_derivs.emplace_back(std::make_unique<WrappedFunctor<Gradient>>(system, gptr, &_equation_systems.parameters));
    }

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = nullptr;
}


void ExactSolution::attach_exact_derivs (const std::vector<FunctionBase<Gradient> *> & g)
{
  // Automatically delete any previous _exact_derivs entries, then add a new
  // entry for each system.
  _exact_derivs.clear();

  for (auto ptr : g)
    _exact_derivs.emplace_back(ptr ? std::make_unique<WrappedFunctor<Gradient>>(*ptr) : nullptr);
}


void ExactSolution::attach_exact_derivs (const std::vector<FEMFunctionBase<Gradient> *> & g)
{
  // Automatically delete any previous _exact_derivs entries, then add a new
  // entry for each system.
  _exact_derivs.clear();

  for (auto ptr : g)
    _exact_derivs.emplace_back(ptr ? ptr->clone() : nullptr);
}


void ExactSolution::attach_exact_deriv (unsigned int sys_num,
                                        FunctionBase<Gradient> * g)
{
  if (_exact_derivs.size() <= sys_num)
    _exact_derivs.resize(sys_num+1);

  if (g)
    _exact_derivs[sys_num] = std::make_unique<WrappedFunctor<Gradient>>(*g);
}


void ExactSolution::attach_exact_deriv (unsigned int sys_num,
                                        FEMFunctionBase<Gradient> * g)
{
  if (_exact_derivs.size() <= sys_num)
    _exact_derivs.resize(sys_num+1);

  if (g)
    _exact_derivs[sys_num] = g->clone();
}


void ExactSolution::attach_exact_hessian (HessianFunctionPointer hptr)
{
  libmesh_assert(hptr);

  // Clear out any previous _exact_hessians entries, then add a new
  // entry for each system.
  _exact_hessians.clear();

  for (auto sys : make_range(_equation_systems.n_systems()))
    {
      const System & system = _equation_systems.get_system(sys);
      _exact_hessians.emplace_back(std::make_unique<WrappedFunctor<Tensor>>(system, hptr, &_equation_systems.parameters));
    }

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = nullptr;
}


void ExactSolution::attach_exact_hessians (std::vector<FunctionBase<Tensor> *> h)
{
  // Automatically delete any previous _exact_hessians entries, then add a new
  // entry for each system.
  _exact_hessians.clear();

  for (auto ptr : h)
    _exact_hessians.emplace_back(ptr ? std::make_unique<WrappedFunctor<Tensor>>(*ptr) : nullptr);
}


void ExactSolution::attach_exact_hessians (std::vector<FEMFunctionBase<Tensor> *> h)
{
  // Automatically delete any previous _exact_hessians entries, then add a new
  // entry for each system.
  _exact_hessians.clear();

  for (auto ptr : h)
    _exact_hessians.emplace_back(ptr ? ptr->clone() : nullptr);
}


void ExactSolution::attach_exact_hessian (unsigned int sys_num,
                                          FunctionBase<Tensor> * h)
{
  if (_exact_hessians.size() <= sys_num)
    _exact_hessians.resize(sys_num+1);

  if (h)
    _exact_hessians[sys_num] = std::make_unique<WrappedFunctor<Tensor>>(*h);
}


void ExactSolution::attach_exact_hessian (unsigned int sys_num,
                                          FEMFunctionBase<Tensor> * h)
{
  if (_exact_hessians.size() <= sys_num)
    _exact_hessians.resize(sys_num+1);

  if (h)
    _exact_hessians[sys_num] = h->clone();
}


std::vector<Real> & ExactSolution::_check_inputs(std::string_view sys_name,
                                                 std::string_view unknown_name)
{
  // Return a reference to the proper error entry, or throw an error
  // if it doesn't exist.
  auto & system_error_map = libmesh_map_find(_errors, sys_name);
  return libmesh_map_find(system_error_map, unknown_name);
}



void ExactSolution::compute_error(std::string_view sys_name,
                                  std::string_view unknown_name)
{
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real> & error_vals = this->_check_inputs(sys_name,
                                                       unknown_name);

  libmesh_assert( _equation_systems.has_system(sys_name) );
  const System & sys = _equation_systems.get_system<System>( sys_name );

  libmesh_assert( sys.has_variable( unknown_name ) );
  switch( FEInterface::field_type(sys.variable_type( unknown_name )) )
    {
    case TYPE_SCALAR:
      {
        this->_compute_error<Real>(sys_name,
                                   unknown_name,
                                   error_vals);
        break;
      }
    case TYPE_VECTOR:
      {
        this->_compute_error<RealGradient>(sys_name,
                                           unknown_name,
                                           error_vals);
        break;
      }
    default:
      libmesh_error_msg("Invalid variable type!");
    }
}





Real ExactSolution::error_norm(std::string_view sys_name,
                               std::string_view unknown_name,
                               const FEMNormType & norm)
{
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real> & error_vals = this->_check_inputs(sys_name,
                                                       unknown_name);

  libmesh_assert(_equation_systems.has_system(sys_name));
  libmesh_assert(_equation_systems.get_system(sys_name).has_variable( unknown_name ));
  const FEType & fe_type = _equation_systems.get_system(sys_name).variable_type(unknown_name);

  switch (norm)
    {
    case L2:
      return std::sqrt(error_vals[0]);
    case H1:
      return std::sqrt(error_vals[0] + error_vals[1]);
    case H2:
      return std::sqrt(error_vals[0] + error_vals[1] + error_vals[2]);
    case HCURL:
      {
        libmesh_error_msg_if(FEInterface::field_type(fe_type) == TYPE_SCALAR,
                             "Cannot compute HCurl error norm of scalar-valued variables!");

        return std::sqrt(error_vals[0] + error_vals[5]);
      }
    case HDIV:
      {
        libmesh_error_msg_if(FEInterface::field_type(fe_type) == TYPE_SCALAR,
                             "Cannot compute HDiv error norm of scalar-valued variables!");

        return std::sqrt(error_vals[0] + error_vals[6]);
      }
    case H1_SEMINORM:
      return std::sqrt(error_vals[1]);
    case H2_SEMINORM:
      return std::sqrt(error_vals[2]);
    case HCURL_SEMINORM:
      {
        libmesh_error_msg_if(FEInterface::field_type(fe_type) == TYPE_SCALAR,
                             "Cannot compute HCurl error seminorm of scalar-valued variables!");

        return std::sqrt(error_vals[5]);
      }
    case HDIV_SEMINORM:
      {
        libmesh_error_msg_if(FEInterface::field_type(fe_type) == TYPE_SCALAR,
                             "Cannot compute HDiv error seminorm of scalar-valued variables!");

        return std::sqrt(error_vals[6]);
      }
    case L1:
      return error_vals[3];
    case L_INF:
      return error_vals[4];

    default:
      libmesh_error_msg("Currently only Sobolev norms/seminorms are supported!");
    }
}







Real ExactSolution::l2_error(std::string_view sys_name,
                             std::string_view unknown_name)
{

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real> & error_vals = this->_check_inputs(sys_name,
                                                       unknown_name);

  // Return the square root of the first component of the
  // computed error.
  return std::sqrt(error_vals[0]);
}







Real ExactSolution::l1_error(std::string_view sys_name,
                             std::string_view unknown_name)
{

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real> & error_vals = this->_check_inputs(sys_name,
                                                       unknown_name);

  // Return the square root of the first component of the
  // computed error.
  return error_vals[3];
}







Real ExactSolution::l_inf_error(std::string_view sys_name,
                                std::string_view unknown_name)
{

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real> & error_vals = this->_check_inputs(sys_name,
                                                       unknown_name);

  // Return the square root of the first component of the
  // computed error.
  return error_vals[4];
}







Real ExactSolution::h1_error(std::string_view sys_name,
                             std::string_view unknown_name)
{
  // If the user has supplied no exact derivative function, we
  // just integrate the H1 norm of the solution; i.e. its
  // difference from an "exact solution" of zero.

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real> & error_vals = this->_check_inputs(sys_name,
                                                       unknown_name);

  // Return the square root of the sum of the computed errors.
  return std::sqrt(error_vals[0] + error_vals[1]);
}


Real ExactSolution::hcurl_error(std::string_view sys_name,
                                std::string_view unknown_name)
{
  return this->error_norm(sys_name,unknown_name,HCURL);
}


Real ExactSolution::hdiv_error(std::string_view sys_name,
                               std::string_view unknown_name)
{
  return this->error_norm(sys_name,unknown_name,HDIV);
}



Real ExactSolution::h2_error(std::string_view sys_name,
                             std::string_view unknown_name)
{
  // If the user has supplied no exact derivative functions, we
  // just integrate the H2 norm of the solution; i.e. its
  // difference from an "exact solution" of zero.

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real> & error_vals = this->_check_inputs(sys_name,
                                                       unknown_name);

  // Return the square root of the sum of the computed errors.
  return std::sqrt(error_vals[0] + error_vals[1] + error_vals[2]);
}







template<typename OutputShape>
void ExactSolution::_compute_error(std::string_view sys_name,
                                   std::string_view unknown_name,
                                   std::vector<Real> & error_vals)
{
  // Make sure we aren't "overconfigured"
  libmesh_assert (!(_exact_values.size() && _equation_systems_fine));

  // We need a communicator.
  const Parallel::Communicator & communicator(_equation_systems.comm());

  // This function must be run on all processors at once
  libmesh_parallel_only(communicator);

  // Get a reference to the system whose error is being computed.
  // If we have a fine grid, however, we'll integrate on that instead
  // for more accuracy.
  const System & computed_system = _equation_systems_fine ?
    _equation_systems_fine->get_system(sys_name) :
    _equation_systems.get_system (sys_name);

  FEMContext context(computed_system);

  const MeshBase & mesh = computed_system.get_mesh();

  const Real time = _equation_systems.get_system(sys_name).time;

  const unsigned int sys_num = computed_system.number();
  const unsigned int var = computed_system.variable_number(unknown_name);
  unsigned int var_component = 0;
  for (const auto var_num : make_range(var))
  {
    const auto & var_fe_type = computed_system.variable_type(var_num);
    const auto var_vec_dim = FEInterface::n_vec_dim(mesh, var_fe_type);
    var_component += var_vec_dim;
  }

  // Prepare a global solution, a serialized mesh, and a MeshFunction
  // of the coarse system, if we need them for fine-system integration
  std::unique_ptr<MeshFunction> coarse_values;
  std::unique_ptr<NumericVector<Number>> comparison_soln =
    NumericVector<Number>::build(_equation_systems.comm());
  MeshSerializer
    serial(const_cast<MeshBase&>(_equation_systems.get_mesh()),
           _equation_systems_fine);
  if (_equation_systems_fine)
    {
      const System & comparison_system
        = _equation_systems.get_system(sys_name);

      std::vector<Number> global_soln;
      comparison_system.update_global_solution(global_soln);
      comparison_soln->init(comparison_system.solution->size(), true, SERIAL);
      (*comparison_soln) = global_soln;

      coarse_values = std::make_unique<MeshFunction>
        (_equation_systems,
         *comparison_soln,
         comparison_system.get_dof_map(),
         comparison_system.variable_number(unknown_name));
      coarse_values->init();
    }

  // Grab which element dimensions are present in the mesh
  const std::set<unsigned char> & elem_dims = mesh.elem_dimensions();

  // Initialize any functors we're going to use
  for (auto & ev : _exact_values)
    if (ev)
      {
        ev->init();
        ev->init_context(context);
      }

  for (auto & ed : _exact_derivs)
    if (ed)
      {
        ed->init();
        ed->init_context(context);
      }

  for (auto & eh : _exact_hessians)
    if (eh)
      {
        eh->init();
        eh->init_context(context);
      }

  // If we have *no* functors we intend to use (because we're using a
  // fine system, or because our exact solution is zero and we're just
  // computing norms in an outdated way) then let our FE objects know
  // we don't actually need anything from them, so they don't think
  // we've just invoked them in a deprecated "compute everything"
  // fashion.
  if (_exact_values.empty() && _exact_derivs.empty() &&
      _exact_hessians.empty())
    for (auto dim : elem_dims)
      for (auto v : make_range(computed_system.n_vars()))
        {
          FEAbstract * fe;
          context.get_element_fe(v, fe, dim);
          fe->get_nothing();
        }

  // Get a reference to the dofmap and mesh for that system
  const DofMap & computed_dof_map = computed_system.get_dof_map();

  // Zero the error before summation
  // 0 - sum of square of function error (L2)
  // 1 - sum of square of gradient error (H1 semi)
  // 2 - sum of square of Hessian error (H2 semi)
  // 3 - sum of sqrt(square of function error) (L1)
  // 4 - max of sqrt(square of function error) (Linfty)
  // 5 - sum of square of curl error (HCurl semi)
  // 6 - sum of square of div error (HDiv semi)
  error_vals = std::vector<Real>(7, 0.);

  // Construct Quadrature rule based on default quadrature order
  const FEType & fe_type  = computed_dof_map.variable_type(var);
  const auto field_type = FEInterface::field_type(fe_type);

  unsigned int n_vec_dim = FEInterface::n_vec_dim( mesh, fe_type );

  // FIXME: MeshFunction needs to be updated to support vector-valued
  //        elements before we can use a reference solution.
  if ((n_vec_dim > 1) && _equation_systems_fine)
    {
      libMesh::err << "Error calculation using reference solution not yet\n"
                   << "supported for vector-valued elements."
                   << std::endl;
      libmesh_not_implemented();
    }


  // Allow space for dims 0-3, even if we don't use them all
  std::vector<std::unique_ptr<FEGenericBase<OutputShape>>> fe_ptrs(4);
  std::vector<std::unique_ptr<QBase>> q_rules(4);

  // Prepare finite elements for each dimension present in the mesh
  for (const auto dim : elem_dims)
    {
      // Build a quadrature rule.
      q_rules[dim] = fe_type.default_quadrature_rule (dim, _extra_order);

      // Disallow rules with negative weights.  That will use more
      // quadrature points, but we're going to be taking square roots
      // of element integral results here!
      q_rules[dim]->allow_rules_with_negative_weights = false;

      // Construct finite element object
      fe_ptrs[dim] = FEGenericBase<OutputShape>::build(dim, fe_type);

      // Attach quadrature rule to FE object
      fe_ptrs[dim]->attach_quadrature_rule (q_rules[dim].get());
    }

  // The global degree of freedom indices associated
  // with the local degrees of freedom.
  std::vector<dof_id_type> dof_indices;


  //
  // Begin the loop over the elements
  //
  // TODO: this ought to be threaded (and using subordinate
  // MeshFunction objects in each thread rather than a single
  // master)
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Skip this element if it is in a subdomain excluded by the user.
      const subdomain_id_type elem_subid = elem->subdomain_id();
      if (_excluded_subdomains.count(elem_subid))
        continue;

      // The spatial dimension of the current Elem. FEs and other data
      // are indexed on dim.
      const unsigned int dim = elem->dim();

      // If the variable is not active on this subdomain, don't bother
      if (!computed_system.variable(var).active_on_subdomain(elem_subid))
        continue;

      /* If the variable is active, then we're going to restrict the
         MeshFunction evaluations to the current element subdomain.
         This is for cases such as mixed dimension meshes where we want
         to restrict the calculation to one particular domain. */
      std::set<subdomain_id_type> subdomain_id;
      subdomain_id.insert(elem_subid);

      FEGenericBase<OutputShape> * fe = fe_ptrs[dim].get();
      QBase * qrule = q_rules[dim].get();
      libmesh_assert(fe);
      libmesh_assert(qrule);

      // The Jacobian*weight at the quadrature points.
      const std::vector<Real> & JxW = fe->get_JxW();

      // The value of the shape functions at the quadrature points
      // i.e. phi(i) = phi_values[i][qp]
      const std::vector<std::vector<OutputShape>> &  phi_values = fe->get_phi();

      // The value of the shape function gradients at the quadrature points
      const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient>> &
        dphi_values = fe->get_dphi();

      // The value of the shape function curls at the quadrature points
      // Only computed for vector-valued elements
      const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputShape>> * curl_values = nullptr;

      // The value of the shape function divergences at the quadrature points
      // Only computed for vector-valued elements
      const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence>> * div_values = nullptr;

      if (field_type == TYPE_VECTOR)
        {
          curl_values = &fe->get_curl_phi();
          div_values = &fe->get_div_phi();
        }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      // The value of the shape function second derivatives at the quadrature points
      // Not computed for vector-valued elements
      const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor>> *
        d2phi_values = nullptr;

      if (field_type != TYPE_VECTOR)
        d2phi_values = &fe->get_d2phi();
#endif

      // The XYZ locations (in physical space) of the quadrature points
      const std::vector<Point> & q_point = fe->get_xyz();

      // reinitialize the element-specific data
      // for the current element
      fe->reinit (elem);

      context.pre_fe_reinit(computed_system, elem);
      context.elem_fe_reinit();

      // Get the local to global degree of freedom maps
      computed_dof_map.dof_indices    (elem, dof_indices, var);

      // The number of quadrature points
      const unsigned int n_qp = qrule->n_points();

      // The number of shape functions
      const unsigned int n_sf =
        cast_int<unsigned int>(dof_indices.size());

      //
      // Begin the loop over the Quadrature points.
      //
      for (unsigned int qp=0; qp<n_qp; qp++)
        {
          // Real u_h = 0.;
          // RealGradient grad_u_h;

          typename FEGenericBase<OutputShape>::OutputNumber u_h(0.);

          typename FEGenericBase<OutputShape>::OutputNumberGradient grad_u_h;
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          typename FEGenericBase<OutputShape>::OutputNumberTensor grad2_u_h;
#endif
          typename FEGenericBase<OutputShape>::OutputNumber curl_u_h(0.0);
          typename FEGenericBase<OutputShape>::OutputNumberDivergence div_u_h = 0.0;

          // Compute solution values at the current
          // quadrature point.  This requires a sum
          // over all the shape functions evaluated
          // at the quadrature point.
          for (unsigned int i=0; i<n_sf; i++)
            {
              // Values from current solution.
              u_h      += phi_values[i][qp]*computed_system.current_solution  (dof_indices[i]);
              grad_u_h += dphi_values[i][qp]*computed_system.current_solution (dof_indices[i]);
              if (field_type == TYPE_VECTOR)
                {
                  curl_u_h += (*curl_values)[i][qp]*computed_system.current_solution (dof_indices[i]);
                  div_u_h += (*div_values)[i][qp]*computed_system.current_solution (dof_indices[i]);
                }
              else
                {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                  grad2_u_h += (*d2phi_values)[i][qp]*computed_system.current_solution (dof_indices[i]);
#endif
                }
            }

          // Compute the value of the error at this quadrature point
          typename FEGenericBase<OutputShape>::OutputNumber exact_val(0);
          RawAccessor<typename FEGenericBase<OutputShape>::OutputNumber> exact_val_accessor( exact_val, n_vec_dim );
          if (_exact_values.size() > sys_num && _exact_values[sys_num])
            {
              for (unsigned int c = 0; c < n_vec_dim; c++)
                exact_val_accessor(c) =
                  _exact_values[sys_num]->
                  component(context, var_component+c, q_point[qp], time);
            }
          else if (_equation_systems_fine)
            {
              // FIXME: Needs to be updated for vector-valued elements
              DenseVector<Number> output(1);
              (*coarse_values)(q_point[qp],time,output,&subdomain_id);
              exact_val = output(0);
            }
          const typename FEGenericBase<OutputShape>::OutputNumber val_error = u_h - exact_val;

          // Add the squares of the error to each contribution
          Real error_sq = TensorTools::norm_sq(val_error);
          error_vals[0] += JxW[qp]*error_sq;

          Real norm = sqrt(error_sq);
          error_vals[3] += JxW[qp]*norm;

          if (error_vals[4]<norm) { error_vals[4] = norm; }

          // Compute the value of the error in the gradient at this
          // quadrature point
          typename FEGenericBase<OutputShape>::OutputNumberGradient exact_grad;
          RawAccessor<typename FEGenericBase<OutputShape>::OutputNumberGradient> exact_grad_accessor( exact_grad, LIBMESH_DIM );
          if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
            {
              for (unsigned int c = 0; c < n_vec_dim; c++)
                for (unsigned int d = 0; d < LIBMESH_DIM; d++)
                  exact_grad_accessor(d + c*LIBMESH_DIM) =
                    _exact_derivs[sys_num]->
                    component(context, var_component+c, q_point[qp], time)(d);
            }
          else if (_equation_systems_fine)
            {
              // FIXME: Needs to be updated for vector-valued elements
              std::vector<Gradient> output(1);
              coarse_values->gradient(q_point[qp],time,output,&subdomain_id);
              exact_grad = output[0];
            }

          const typename FEGenericBase<OutputShape>::OutputNumberGradient grad_error = grad_u_h - exact_grad;

          error_vals[1] += JxW[qp]*grad_error.norm_sq();


          if (field_type == TYPE_VECTOR)
            {
              // Compute the value of the error in the curl at this
              // quadrature point
              typename FEGenericBase<OutputShape>::OutputNumber exact_curl(0.0);
              if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
                {
                  exact_curl = TensorTools::curl_from_grad( exact_grad );
                }
              else if (_equation_systems_fine)
                {
                  // FIXME: Need to implement curl for MeshFunction and support reference
                  //        solution for vector-valued elements
                }

              const typename FEGenericBase<OutputShape>::OutputNumber curl_error = curl_u_h - exact_curl;

              error_vals[5] += JxW[qp]*TensorTools::norm_sq(curl_error);

              // Compute the value of the error in the divergence at this
              // quadrature point
              typename FEGenericBase<OutputShape>::OutputNumberDivergence exact_div = 0.0;
              if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
                {
                  exact_div = TensorTools::div_from_grad( exact_grad );
                }
              else if (_equation_systems_fine)
                {
                  // FIXME: Need to implement div for MeshFunction and support reference
                  //        solution for vector-valued elements
                }

              const typename FEGenericBase<OutputShape>::OutputNumberDivergence div_error = div_u_h - exact_div;

              error_vals[6] += JxW[qp]*TensorTools::norm_sq(div_error);
            }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          // Compute the value of the error in the hessian at this
          // quadrature point
          typename FEGenericBase<OutputShape>::OutputNumberTensor exact_hess;
          RawAccessor<typename FEGenericBase<OutputShape>::OutputNumberTensor> exact_hess_accessor( exact_hess, dim );
          if (_exact_hessians.size() > sys_num && _exact_hessians[sys_num])
            {
              //FIXME: This needs to be implemented to support rank 3 tensors
              //       which can't happen until type_n_tensor is fully implemented
              //       and a RawAccessor<TypeNTensor> is fully implemented
              if (field_type == TYPE_VECTOR)
                libmesh_not_implemented();

              for (unsigned int c = 0; c < n_vec_dim; c++)
                for (unsigned int d = 0; d < dim; d++)
                  for (unsigned int e =0; e < dim; e++)
                    exact_hess_accessor(d + e*dim + c*dim*dim) =
                      _exact_hessians[sys_num]->
                      component(context, var_component+c, q_point[qp], time)(d,e);

              // FIXME: operator- is not currently implemented for TypeNTensor
              const typename FEGenericBase<OutputShape>::OutputNumberTensor grad2_error = grad2_u_h - exact_hess;
              error_vals[2] += JxW[qp]*grad2_error.norm_sq();
            }
          else if (_equation_systems_fine)
            {
              // FIXME: Needs to be updated for vector-valued elements
              std::vector<Tensor> output(1);
              coarse_values->hessian(q_point[qp],time,output,&subdomain_id);
              exact_hess = output[0];

              // FIXME: operator- is not currently implemented for TypeNTensor
              const typename FEGenericBase<OutputShape>::OutputNumberTensor grad2_error = grad2_u_h - exact_hess;
              error_vals[2] += JxW[qp]*grad2_error.norm_sq();
            }
#endif

        } // end qp loop
    } // end element loop

  // Add up the error values on all processors, except for the L-infty
  // norm, for which the maximum is computed.
  Real l_infty_norm = error_vals[4];
  communicator.max(l_infty_norm);
  communicator.sum(error_vals);
  error_vals[4] = l_infty_norm;
}

// Explicit instantiations of templated member functions
template LIBMESH_EXPORT void ExactSolution::_compute_error<Real>(std::string_view, std::string_view, std::vector<Real> &);
template LIBMESH_EXPORT void ExactSolution::_compute_error<RealGradient>(std::string_view, std::string_view, std::vector<Real> &);

} // namespace libMesh
