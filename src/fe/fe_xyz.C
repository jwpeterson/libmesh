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
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"

namespace libMesh
{

// ------------------------------------------------------------
// XYZ-specific implementations

// Anonymous namespace for local helper functions
namespace {

void xyz_nodal_soln(const Elem * elem,
                    const Order order,
                    const std::vector<Number> & elem_soln,
                    std::vector<Number> & nodal_soln,
                    const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  nodal_soln.resize(n_nodes);

  const Order totalorder = order + add_p_level*elem->p_level();

  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
        libmesh_assert_equal_to (elem_soln.size(), 1);

        std::fill(nodal_soln.begin(), nodal_soln.end(), elem_soln[0]);

        return;
      }


      // For other orders do interpolation at the nodes
      // explicitly.
    default:
      {
        // FEType object to be passed to various FEInterface functions below.
        FEType fe_type(order, XYZ);

        const unsigned int n_sf =
          FEInterface::n_shape_functions(fe_type, elem);
        libmesh_assert_equal_to (elem_soln.size(), n_sf);

        // Zero before summation
        std::fill(nodal_soln.begin(), nodal_soln.end(), 0);

        for (unsigned int n=0; n<n_nodes; n++)
          // u_i = Sum (alpha_i phi_i)
          for (unsigned int i=0; i<n_sf; i++)
            nodal_soln[n] += elem_soln[i] *
              FEInterface::shape(fe_type, elem, i, elem->point(n));

        return;
      } // default
    } // switch
} // xyz_nodal_soln()


} // anonymous namespace







template <unsigned int Dim>
void FEXYZ<Dim>::init_shape_functions(const std::vector<Point> & qp,
                                      const Elem * elem)
{
  libmesh_assert(elem);

  // FIXME: Is this redundant here? Who's calling init_shape_functions
  // from code that hasn't already done a determine_calculations()?
  this->determine_calculations();

  // Start logging the shape function initialization
  LOG_SCOPE("init_shape_functions()", "FE");

  // The number of quadrature points.
  const std::size_t n_qp = qp.size();

  // Number of shape functions in the finite element approximation
  // space.
  const unsigned int n_approx_shape_functions =
    this->n_dofs(elem,
                 this->get_order());

  // resize the vectors to hold current data
  // Phi are the shape functions used for the FE approximation
  {
    // (note: GCC 3.4.0 requires the use of this-> here)
    if (this->calculate_phi)
      this->phi.resize     (n_approx_shape_functions);
    if (this->calculate_dphi)
      {
        this->dphi.resize    (n_approx_shape_functions);
        this->dphidx.resize  (n_approx_shape_functions);
        this->dphidy.resize  (n_approx_shape_functions);
        this->dphidz.resize  (n_approx_shape_functions);
      }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    if (this->calculate_d2phi)
      {
        this->d2phi.resize     (n_approx_shape_functions);
        this->d2phidx2.resize  (n_approx_shape_functions);
        this->d2phidxdy.resize (n_approx_shape_functions);
        this->d2phidxdz.resize (n_approx_shape_functions);
        this->d2phidy2.resize  (n_approx_shape_functions);
        this->d2phidydz.resize (n_approx_shape_functions);
        this->d2phidz2.resize  (n_approx_shape_functions);
        this->d2phidxi2.resize (n_approx_shape_functions);
      }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

    for (unsigned int i=0; i<n_approx_shape_functions; i++)
      {
        if (this->calculate_phi)
          this->phi[i].resize           (n_qp);
        if (this->calculate_dphi)
          {
            this->dphi[i].resize        (n_qp);
            this->dphidx[i].resize      (n_qp);
            this->dphidy[i].resize      (n_qp);
            this->dphidz[i].resize      (n_qp);
          }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          {
            this->d2phi[i].resize       (n_qp);
            this->d2phidx2[i].resize    (n_qp);
            this->d2phidxdy[i].resize   (n_qp);
            this->d2phidxdz[i].resize   (n_qp);
            this->d2phidy2[i].resize    (n_qp);
            this->d2phidydz[i].resize   (n_qp);
            this->d2phidz2[i].resize    (n_qp);
          }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      }
  }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  //------------------------------------------------------------
  // Initialize the data fields, which should only be used for infinite
  // elements, to some sensible values, so that using a FE with the
  // variational formulation of an InfFE, correct element matrices are
  // returned

  {
    this->weight.resize  (n_qp);
    this->dweight.resize (n_qp);
    this->dphase.resize  (n_qp);

    for (unsigned int p=0; p<n_qp; p++)
      {
        this->weight[p] = 1.;
        this->dweight[p].zero();
        this->dphase[p].zero();
      }

  }
#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (this->calculate_dual)
    this->init_dual_shape_functions(n_approx_shape_functions, n_qp);
}




template <unsigned int Dim>
void FEXYZ<Dim>::compute_shape_functions (const Elem * elem,
                                          const std::vector<Point> &)
{
  libmesh_assert(elem);

  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
  LOG_SCOPE("compute_shape_functions()", "FE");

  const std::vector<Point> & xyz_qp = this->get_xyz();

  // Compute the value of the derivative shape function i at quadrature point p
  switch (this->dim)
    {

    case 1:
      {
        if (this->calculate_phi)
          for (auto i : index_range(this->phi))
            for (auto p : index_range(this->phi[i]))
              this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz_qp[p]);

        if (this->calculate_dphi)
          for (auto i : index_range(this->dphi))
            for (auto p : index_range(this->dphi[i]))
              {
                this->dphi[i][p](0) =
                  this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->dphi[i][p](1) = this->dphidy[i][p] = 0.;
                this->dphi[i][p](2) = this->dphidz[i][p] = 0.;
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          for (auto i : index_range(this->d2phi))
            for (auto p : index_range(this->d2phi[i]))
              {
                this->d2phi[i][p](0,0) =
                  this->d2phidx2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

#if LIBMESH_DIM>1
                this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
                  this->d2phi[i][p](1,0) = 0.;
                this->d2phi[i][p](1,1) = this->d2phidy2[i][p] = 0.;
#if LIBMESH_DIM>2
                this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
                  this->d2phi[i][p](2,0) = 0.;
                this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
                  this->d2phi[i][p](2,1) = 0.;
                this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = 0.;
#endif
#endif
              }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        // All done
        break;
      }

    case 2:
      {
        if (this->calculate_phi)
          for (auto i : index_range(this->phi))
            for (auto p : index_range(this->phi[i]))
              this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz_qp[p]);

        if (this->calculate_dphi)
          for (auto i : index_range(this->dphi))
            for (auto p : index_range(this->dphi[i]))
              {
                this->dphi[i][p](0) =
                  this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->dphi[i][p](1) =
                  this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);

#if LIBMESH_DIM == 3
                this->dphi[i][p](2) = // can only assign to the Z component if LIBMESH_DIM==3
#endif
                  this->dphidz[i][p] = 0.;
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          for (auto i : index_range(this->d2phi))
            for (auto p : index_range(this->d2phi[i]))
              {
                this->d2phi[i][p](0,0) =
                  this->d2phidx2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
                  this->d2phi[i][p](1,0) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
                this->d2phi[i][p](1,1) =
                  this->d2phidy2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
#if LIBMESH_DIM>2
                this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
                  this->d2phi[i][p](2,0) = 0.;
                this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
                  this->d2phi[i][p](2,1) = 0.;
                this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = 0.;
#endif
              }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        // All done
        break;
      }

    case 3:
      {
        if (this->calculate_phi)
          for (auto i : index_range(this->phi))
            for (auto p : index_range(this->phi[i]))
              this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz_qp[p]);

        if (this->calculate_dphi)
          for (auto i : index_range(this->dphi))
            for (auto p : index_range(this->dphi[i]))
              {
                this->dphi[i][p](0) =
                  this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->dphi[i][p](1) =
                  this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);

                this->dphi[i][p](2) =
                  this->dphidz[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          for (auto i : index_range(this->d2phi))
            for (auto p : index_range(this->d2phi[i]))
              {
                this->d2phi[i][p](0,0) =
                  this->d2phidx2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
                  this->d2phi[i][p](1,0) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
                this->d2phi[i][p](1,1) =
                  this->d2phidy2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
                this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
                  this->d2phi[i][p](2,0) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 3, xyz_qp[p]);
                this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
                  this->d2phi[i][p](2,1) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 4, xyz_qp[p]);
                this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 5, xyz_qp[p]);
              }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        // All done
        break;
      }

    default:
      libmesh_error_msg("ERROR: Invalid dimension " << this->dim);
    }
}


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(XYZ, xyz_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(XYZ)


// Full specialization of n_dofs() function for every dimension
template <> unsigned int FE<0,XYZ>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<1,XYZ>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<2,XYZ>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<3,XYZ>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }

template <> unsigned int FE<0,XYZ>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<1,XYZ>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<2,XYZ>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }
template <> unsigned int FE<3,XYZ>::n_dofs(const Elem * e, const Order o) { return monomial_n_dofs(e, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
// XYZ FEMs have no dofs at nodes, only element dofs.
template <> unsigned int FE<0,XYZ>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,XYZ>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,XYZ>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,XYZ>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,XYZ>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,XYZ>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,XYZ>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,XYZ>::n_dofs_at_node(const Elem &, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <> unsigned int FE<0,XYZ>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<1,XYZ>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<2,XYZ>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
template <> unsigned int FE<3,XYZ>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }

template <> unsigned int FE<0,XYZ>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<1,XYZ>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<2,XYZ>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }
template <> unsigned int FE<3,XYZ>::n_dofs_per_elem(const Elem & e, const Order o) { return monomial_n_dofs(&e, o); }

// Full specialization of get_continuity() function for every dimension.
template <> FEContinuity FE<0,XYZ>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,XYZ>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,XYZ>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,XYZ>::get_continuity() const { return DISCONTINUOUS; }

// Full specialization of is_hierarchic() function for every dimension.
// The XYZ shape functions are hierarchic!
template <> bool FE<0,XYZ>::is_hierarchic() const { return true; }
template <> bool FE<1,XYZ>::is_hierarchic() const { return true; }
template <> bool FE<2,XYZ>::is_hierarchic() const { return true; }
template <> bool FE<3,XYZ>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR

// Full specialization of compute_constraints() function for 2D and
// 3D only.  There are no constraints for discontinuous elements, so
// we do nothing.
template <> void FE<2,XYZ>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem *) {}
template <> void FE<3,XYZ>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem *) {}

#endif // #ifdef LIBMESH_ENABLE_AMR

// Full specialization of shapes_need_reinit() function for every dimension.
template <> bool FE<0,XYZ>::shapes_need_reinit() const { return true; }
template <> bool FE<1,XYZ>::shapes_need_reinit() const { return true; }
template <> bool FE<2,XYZ>::shapes_need_reinit() const { return true; }
template <> bool FE<3,XYZ>::shapes_need_reinit() const { return true; }


// Explicit instantiations for non-static FEXYZ member functions.
// These non-static member functions map more naturally to explicit
// instantiations than the functions above:
//
// 1.)  Since they are member functions, they rely on
// private/protected member data, and therefore do not work well
// with the "anonymous function call" model we've used above for
// the specializations.
//
// 2.) There is (IMHO) less chance of the linker calling the
// wrong version of one of these member functions, since there is
// only one FEXYZ.
template LIBMESH_EXPORT void  FEXYZ<0>::init_shape_functions(const std::vector<Point> &, const Elem *);
template LIBMESH_EXPORT void  FEXYZ<1>::init_shape_functions(const std::vector<Point> &, const Elem *);
template LIBMESH_EXPORT void  FEXYZ<2>::init_shape_functions(const std::vector<Point> &, const Elem *);
template LIBMESH_EXPORT void  FEXYZ<3>::init_shape_functions(const std::vector<Point> &, const Elem *);

template LIBMESH_EXPORT void  FEXYZ<0>::compute_shape_functions(const Elem *,const std::vector<Point> &);
template LIBMESH_EXPORT void  FEXYZ<1>::compute_shape_functions(const Elem *,const std::vector<Point> &);
template LIBMESH_EXPORT void  FEXYZ<2>::compute_shape_functions(const Elem *,const std::vector<Point> &);
template LIBMESH_EXPORT void  FEXYZ<3>::compute_shape_functions(const Elem *,const std::vector<Point> &);

} // namespace libMesh
