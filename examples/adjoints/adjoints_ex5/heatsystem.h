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



#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"
#include "libmesh/parameter_pointer.h"
#include "libmesh/parameter_vector.h"

#include <memory>

using namespace libMesh;

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class HeatSystem : public FEMSystem
{
public:
  // Constructor
  HeatSystem(EquationSystems & es,
             const std::string & name_in,
             const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
      _k(1.0),
      _fe_family("LAGRANGE"),
      _fe_order(1),
      _analytic_jacobians(true),
      dp(1.e-6)
  { this->init_qois(1); }

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }
  Real & k() { return _k; }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // A function to compute and accumulate residuals
  void perturb_accumulate_residuals(ParameterVector & parameters);

  // Sensitivity Calculation
  Number & compute_final_sensitivity()
  {
    final_sensitivity = 0.0;

    // Use the trapezoidal rule to compute the sensitivity integral
    for(unsigned int i = 0; i < R_plus_dp.size()-1; i++)
    {
      Number left_contribution = -(R_plus_dp[i] - R_minus_dp[i])/(2.*dp);
      Number right_contribution = -(R_plus_dp[i+1] - R_minus_dp[i+1])/(2.*dp);

      final_sensitivity += ( (left_contribution + right_contribution)/2. )*deltat_vector[i];
    }

    return final_sensitivity;
  }

  void set_tf(Real val)
  {
    tf = val;
  }

  ParameterVector & get_parameter_vector()
  {
    if (!parameter_vector.size())
      for (std::size_t i = 0; i != parameters.size(); ++i)
        parameter_vector.push_back(std::make_unique<ParameterPointer<Number>>(&parameters[i]));

    return parameter_vector;
  }

  Number & get_QoI_value(unsigned int QoI_index)
  {
    return computed_QoI[QoI_index];
  }

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  // Constraint parts
  // virtual bool side_constraint (bool request_jacobian,
  //                               DiffContext & context);

  // RHS for adjoint problem
  virtual void element_qoi_derivative (DiffContext & context,
                                       const QoISet & /* qois */);

  //virtual void element_qoi (DiffContext & context, const QoISet & qois);

  // Parameters associated with the system
  std::vector<Number> parameters;

  // The ParameterVector object that will contain pointers to
  // the system parameters
  ParameterVector parameter_vector;

  // The parameters to solve for
  Real _k;

  // The final time parameter
  Real tf;

  // Variables to hold the computed QoIs
  Number computed_QoI[1];

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Index for T variable
  unsigned int T_var;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;

  // Variables to hold the perturbed residuals
  std::vector<Number> R_plus_dp;
  std::vector<Number> R_minus_dp;

  // A vector to hold the possibly adaptive deltats
  std::vector<Number> deltat_vector;

  // Perturbation parameter
  Number dp;

  // The final computed sensitivity
  Number final_sensitivity;
};
