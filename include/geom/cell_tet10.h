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



#ifndef LIBMESH_CELL_TET10_H
#define LIBMESH_CELL_TET10_H

// Local includes
#include "libmesh/cell_tet.h"

namespace libMesh
{

/**
 * The \p Tet10 is an element in 3D composed of 10 nodes.  It is
 * numbered like this:
 *
 * \verbatim
 *               3
 *   TET10:      o
 *              /|\
 *             / | \
 *         7  /  |  \9
 *           o   |   o              zeta
 *          /    |8   \              ^
 *         /     o     \             |
 *        /    6 |      \            |
 *     0 o.....o.|.......o 2         o---> eta
 *        \      |      /             \
 *         \     |     /               \
 *          \    |    /                 xi (out of page)
 *         4 o   |   o 5
 *            \  |  /
 *             \ | /
 *              \|/
 *               o
 *               1
 * \endverbatim
 *
 * (xi, eta, zeta): { 0 <= xi   <= 1
 *                  { 0 <= eta  <= 1
 *                  { 0 <= zeta <= 1
 *                  { xi + eta + zeta <= 1
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D tetrahedral element with 10 nodes.
 */
class Tet10 final : public Tet
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tet10 (Elem * p=nullptr) :
    Tet(num_nodes, p, _nodelinks_data)
  {}

  Tet10 (Tet10 &&) = delete;
  Tet10 (const Tet10 &) = delete;
  Tet10 & operator= (const Tet10 &) = delete;
  Tet10 & operator= (Tet10 &&) = delete;
  virtual ~Tet10() = default;

  /**
   * \returns \p TET10.
   */
  virtual ElemType type () const override { return TET10; }

  /**
   * \returns 10.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 8.
   */
  virtual unsigned int n_sub_elem() const override { return 8; }

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
   * \returns \p true if the specified child is on the
   * specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override;

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const override;

  /**
   * \returns \p Tet10::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Tet10::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * Builds a \p TRI6 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) override;

  /**
   * Rebuilds a TRI6 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  // Avoid hiding deprecated version with different signature
  using Elem::build_side_ptr;

  /**
   * Builds a \p EDGE3 built coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p EDGE3 coincident with edge i.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for all \p n.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override
  { return 2; }

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 4 \le n < 10 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const override;

  /**
   * \returns The child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  See
   * elem.h for further details.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * Geometric constants for Tet10.
   */
  static const int num_nodes = 10;
  static const int nodes_per_side = 6;
  static const int nodes_per_edge = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[num_edges][nodes_per_edge];

  /**
   * A specialization for computing the volume of a Tet10.
   */
  virtual Real volume () const override;

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
                                 const unsigned int k) const override;

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const Real _embedding_matrix[num_children][num_nodes][num_nodes];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

private:

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[6][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[10];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[10];
};

} // namespace libMesh


#endif // LIBMESH_CELL_TET10_H
