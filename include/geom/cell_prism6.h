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



#ifndef LIBMESH_CELL_PRISM6_H
#define LIBMESH_CELL_PRISM6_H

// Local includes
#include "libmesh/cell_prism.h"

namespace libMesh
{

/**
 * The \p Prism6 is an element in 3D composed of 6 nodes.  It is
 * numbered like this:
 *
 * \verbatim
 *   PRISM6:
 *           5
 *           o
 *          /:\
 *         / : \             zeta
 *        /  o  \             ^   eta (into page)
 *     3 o-------o 4          | /
 *       | . 2 . |            |/
 *       |.     .|            o---> xi
 *       o-------o
 *       0       1
 * \endverbatim
 *
 * (xi, eta, zeta): { 0  <= xi   <= 1
 *                  { 0  <= eta  <= 1
 *                  { -1 <= zeta <= 1
 *                  { xi + eta   <= 1
 * are the reference element coordinates associated with the given
 * numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D prismatic element with 6 nodes.
 */
class Prism6 final : public Prism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Prism6 (Elem * p=nullptr) :
    Prism(num_nodes, p, _nodelinks_data)
  {}

  Prism6 (Prism6 &&) = delete;
  Prism6 (const Prism6 &) = delete;
  Prism6 & operator= (const Prism6 &) = delete;
  Prism6 & operator= (Prism6 &&) = delete;
  virtual ~Prism6() = default;

  /**
   * \returns \p PRISM6.
   */
  virtual ElemType type () const override { return PRISM6; }

  /**
   * \returns 1.
   */
  virtual unsigned int n_sub_elem() const override { return 1; }

  /**
   * \returns \p true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int e) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns FIRST.
   */
  virtual Order default_order() const override;

  /**
   * Builds a \p QUAD4 or \p TRI3 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p QUAD4 or \p TRI3 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  // Avoid hiding deprecated version with different signature
  using Elem::build_side_ptr;

  /**
   * Builds a \p EDGE2 or \p INFEDGE2 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p EDGE2 built coincident with edges 0 to 2, or \p
   * INFEDGE2 built coincident with edges 3 to 5.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * Geometric constants for Prism6.
   */
  static const int num_nodes = 6;
  static const int nodes_per_side = 4;
  static const int nodes_per_edge = 2;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * This maps the child elements with the associated side of the parent element
   */
  static const unsigned int side_elems_map[num_sides][nodes_per_side];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[num_edges][nodes_per_edge];

  /**
   * An Optimized numerical quadrature approach for computing the
   * centroid of the Prism6.
   */
  virtual Point true_centroid () const override;

  /**
   * Specialized function for computing the element volume.
   */
  virtual Real volume () const override;

  /**
   * Builds a bounding box out of the nodal positions
   */
  virtual BoundingBox loose_bounding_box () const override;

  virtual void permute(unsigned int perm_num) override final;

  virtual void flip(BoundaryInfo *) override final;

  ElemType side_type (const unsigned int s) const override final;

protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[num_nodes];



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
  static const Real _embedding_matrix[num_children][num_nodes][num_nodes];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh

#endif // LIBMESH_CELL_PRISM6_H
