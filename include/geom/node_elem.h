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



#ifndef LIBMESH_NODE_ELEM_H
#define LIBMESH_NODE_ELEM_H

// Local includes
#include "libmesh/elem.h"

namespace libMesh
{

/**
 * The \p NodeElem is a point element, generally used as
 * a side of a 1D element.
 *
 * \author Roy H. Stogner
 * \date 2006
 * \brief A zero-dimensional geometric entity implementing the Elem interface.
 */
class NodeElem : public Elem
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  NodeElem (Elem * p=nullptr) :
    Elem(num_nodes, num_sides, p, _elemlinks_data,
         _nodelinks_data)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 0)
      this->set_interior_parent(nullptr);
  }

  NodeElem (NodeElem &&) = delete;
  NodeElem (const NodeElem &) = delete;
  NodeElem & operator= (const NodeElem &) = delete;
  NodeElem & operator= (NodeElem &&) = delete;
  virtual ~NodeElem() = default;

  /**
   * Geometric constants for NodeElem;
   */
  static const int num_nodes = 1;
  static const int num_sides = 0;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int libmesh_dbg_var(i)) const override
  {
    libmesh_assert_equal_to (i, 0);
    return Point(0,0,0);
  }

  /**
   * \returns 0, the dimensionality of the object.
   */
  virtual unsigned short dim () const override { return 0; }

  /**
   * \returns 1.
   */
  virtual unsigned int n_nodes() const override { return 1; }

  /**
   * \returns 0.
   */
  virtual unsigned int n_sides() const override { return 0; }

  /**
   * \returns 1.  Every NodeElem is a vertex
   */
  virtual unsigned int n_vertices() const override { return 1; }

  /**
   * \returns 0.
   */
  virtual unsigned int n_edges() const override { return 0; }

  /**
   * \returns 0.
   */
  virtual unsigned int n_faces() const override { return 0; }

  /**
   * \returns 1.
   */
  virtual unsigned int n_children() const override { return 1; }

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * This should never be important for NodeElems.
   */
  virtual dof_id_type key (const unsigned int) const override
  { libmesh_error_msg("Calling NodeElem::key(side) does not make sense."); return 0; }

  /**
   * \returns An id associated with the \p s side of this element.
   * This should never be important for NodeElems.
   */
  virtual dof_id_type low_order_key (const unsigned int) const override
  { libmesh_error_msg("Calling NodeElem::low_order_key(side) does not make sense."); return 0; }

  /**
   * NodeElems don't have sides, so they can't have nodes on sides.
   */
  virtual unsigned int local_side_node(unsigned int /*side*/,
                                       unsigned int /*side_node*/) const override
  { libmesh_error_msg("Calling NodeElem::local_side_node() does not make sense."); return 0; }

  /**
   * NodeElems don't have edges, so they can't have nodes on edges.
   */
  virtual unsigned int local_edge_node(unsigned int /*edge*/,
                                       unsigned int /*edge_node*/) const override
  { libmesh_error_msg("Calling NodeElem::local_edge_node() does not make sense."); return 0; }

  /**
   * The \p Elem::side_ptr() member makes no sense for nodes.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual void side_ptr (std::unique_ptr<Elem> &, const unsigned int) override
  { libmesh_not_implemented(); }

  /**
   * The \p Elem::build_side_ptr() member makes no sense for nodes.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual void build_side_ptr (std::unique_ptr<Elem> &, const unsigned int) override
  { libmesh_not_implemented(); }

  // Avoid hiding deprecated version with different signature
  using Elem::build_side_ptr;

  /**
   * The \p Elem::build_edge_ptr() member makes no sense for nodes.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  /**
   * The \p Elem::build_edge_ptr() member makes no sense for nodes.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> &, const unsigned int) override
  { libmesh_not_implemented(); }

  /**
   * \returns 1.
   */
  virtual unsigned int n_sub_elem() const override { return 1; }

  /**
   * \returns \p true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int libmesh_dbg_var(n)) const override
  { libmesh_assert_not_equal_to (n, invalid_uint); return true; }

  /**
   * NodeElem objects don't have faces or sides.
   */
  virtual bool is_edge(const unsigned int) const override { return false; }
  virtual bool is_face(const unsigned int) const override { return false; }

  virtual bool is_child_on_side(const unsigned int,
                                const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_side(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int) const override
  {
    libmesh_not_implemented();
    return {0};
  }

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int) const override
  {
    libmesh_not_implemented();
    return {0};
  }

  virtual std::vector<unsigned int> edges_adjacent_to_node(const unsigned int) const override
  {
    return {};
  }

  virtual std::vector<unsigned int> sides_on_edge(const unsigned int) const override
  {
    libmesh_not_implemented();
    return {0};
  }

  virtual bool is_node_on_edge(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_edge_on_side(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  /**
   * \returns The "circumcenter of mass" (area-weighted average of
   * triangulation circumcenters) of the element.
   *
   * Trivial in 0D.
   */
  virtual Point quasicircumcenter () const override
  { return this->point(0); }

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override { return true; }

  /**
   * \returns \p true because it doesn't really make sense for a NodeElem.
   */
  virtual bool has_invertible_map(Real /*tol*/) const override { return true; }

  /**
   * \returns \p true if the Lagrange shape functions on this element
   * are linear.
   */
  virtual bool is_linear () const override { return true; }

  /**
   * \returns \p NODEELEM.
   */
  virtual ElemType type() const override { return NODEELEM; }

  /**
   * \returns \p CONSTANT.
   */
  virtual Order default_order() const override;

  /**
   * \returns \p MAXIMUM.
   */
  virtual Order supported_nodal_order() const override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * \returns \p false.
   */
  virtual bool infinite () const override { return false; }

#endif

  /**
   * No way to reorient a single node.
   */
  virtual unsigned int n_permutations() const override final { return 0; }

  virtual void permute(unsigned int) override final { libmesh_error(); }

  virtual void flip(BoundaryInfo *) override final { return; /* no-op */ }
  virtual bool is_flipped() const override final { return false; }

  virtual ElemType side_type (const unsigned int) const override
  {
    libmesh_not_implemented();
    return INVALID_ELEM;
  }

  /**
   * \returns \p true if the point p is within a distance of "tol" from
   * the point representing this element, false otherwise.
   *
   * The tolerance "tol" is normally treated as a relative tolerance in
   * contains_point() checks, but in this case there is no relevant length
   * to use in determing a relative tolerance, so "tol" is treated as an
   * absolute tolerance. The NodeElem contains_point() and close_to_point()
   * implementations are identical, whereas they differ for other element
   * types.
   */
  virtual bool contains_point(const Point & p, Real tol) const override;

  /**
   * \returns this->contains_point(p, tol)
   */
  virtual bool close_to_point(const Point & p, Real tol) const override;

  virtual bool on_reference_element(const Point & p,
                                    const Real eps = TOLERANCE) const override final;

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[1+(LIBMESH_DIM>0)];

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[1];


#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual Real embedding_matrix (const unsigned int i,
                                 const unsigned int j,
                                 const unsigned int k) const override
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const Real _embedding_matrix[1][1][1];

  /**
   * Matrix that allows children to inherit boundary conditions.
   */
  unsigned int side_children_matrix (const unsigned int,
                                     const unsigned int) const
  { libmesh_not_implemented(); return 0; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh

#endif // LIBMESH_NODE_ELEM_H
