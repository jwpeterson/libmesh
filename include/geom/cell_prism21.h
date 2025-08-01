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



#ifndef LIBMESH_CELL_PRISM21_H
#define LIBMESH_CELL_PRISM21_H

// Local includes
#include "libmesh/cell_prism.h"

namespace libMesh
{

/**
 * The \p Prism21 is an element in 3D composed of 21 nodes.
 * It is numbered like this:
 *
 * \verbatim
 *   PRISM21:
 *               5
 *               o
 *              /:\
 *             / : \
 *            /  :  \
 *           /   :   \
 *       14 o    :    o 13
 *         /     :     \
 *        /      :      \
 *       /       o 11    \
 *    3 /        :        \4
 *     o---------o---------o
 *     |         :12       |
 *     |         :         |
 *     |    o    :    o    |            zeta
 *     |   17    o   16    |             ^   eta (into page)
 *     |        .2.        |             | /
 *     |       .   .       |             |/
 *   9 o      .  o  .      o 10          o---> xi
 *     |     .  15   .     |
 *     |  8 o         o 7  |
 *     |   .           .   |
 *     |  .             .  |
 *     | .               . |
 *     |.                 .|
 *     o---------o---------o
 *     0         6         1
 * \endverbatim
 *
 * And it also includes two face nodes:
 * Node 18, centroid on side 0, arithmetic mean of 0/1/2 or 6/7/8
 * Node 19, centroid on side 4, arithmetic mean of 3/4/5 or 12/13/14
 *
 * And one midpoint node:
 * Node 20, on the interior, arithmetic mean of 9/10/11 or 15/16/17
 *
 * (xi, eta, zeta): { 0  <= xi   <= 1
 *                  { 0  <= eta  <= 1
 *                  { -1 <= zeta <= 1
 *                  { xi + eta   <= 1
 * are the reference element coordinates associated with the given
 * numbering.
 *
 * \author Roy H. Stogner
 * \date 2022
 * \brief A 3D prismatic element with 21 nodes.
 */
class Prism21 final : public Prism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Prism21 (Elem * p=nullptr) :
    Prism(num_nodes, p, _nodelinks_data)
  {}

  Prism21 (Prism21 &&) = delete;
  Prism21 (const Prism21 &) = delete;
  Prism21 & operator= (const Prism21 &) = delete;
  Prism21 & operator= (Prism21 &&) = delete;
  virtual ~Prism21() = default;

  /**
   * \returns \p PRISM21.
   */
  virtual ElemType type () const override { return PRISM21; }

  /**
   * \returns 21.
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
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.
   *
   * We reimplement this method here for the \p Prism21 since we can
   * use the center node of each quad face to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p Prism21::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Prism21::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * Builds a \p QUAD9 or \p TRI6 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p QUAD9 or \p TRI6 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  // Avoid hiding deprecated version with different signature
  using Elem::build_side_ptr;

  /**
   * Builds a \p EDGE3 or \p INFEDGE2 built coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  /**
   * Rebuilds a \p EDGE3 built coincident with edges 0 to 2, or \p
   * INFEDGE2 built coincident with edges 3 to 5.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for all edge nodes, 4 for quad face nodes, 3 for tri
   * face nodes, 6 for the interior node.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override;

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 6 \le n < 21 \f$.
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
   * Geometric constants for Prism21.
   */
  static const int num_nodes = 21;
  static const int nodes_per_side = 9;
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

  virtual void permute(unsigned int perm_num) override final;

  virtual void flip(BoundaryInfo *) override final;

#ifdef LIBMESH_ENABLE_AMR
  virtual
  const std::vector<std::pair<unsigned char, unsigned char>> &
  parent_bracketing_nodes(unsigned int c,
                          unsigned int n) const override
  { return _parent_bracketing_nodes[c][n]; }
#endif

  unsigned int center_node_on_side(const unsigned short side) const override final;

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

  /**
   * Pairs of nodes that bracket child nodes when doing mesh
   * refinement.
   */
  static const std::vector<std::pair<unsigned char, unsigned char>>
    _parent_bracketing_nodes[num_children][num_nodes];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.  This matrix
   * handles only the second-order nodes that are unique
   * to \p Prism20.  Earlier second-order nodes are identical
   * with \p Prism15, and are therefore handled through a
   * matrix contained in \p cell_prism.C; the final second-order node
   * is just a simple identity map.
   */
  static const unsigned short int _remaining_second_order_adjacent_vertices[5][4];
};

} // namespace libMesh

#endif // LIBMESH_CELL_PRISM21_H
