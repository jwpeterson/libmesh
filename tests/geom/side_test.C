#include <libmesh/elem.h>

#include <libmesh/cell_hex20.h>
#include <libmesh/cell_hex27.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/cell_inf_hex16.h>
#include <libmesh/cell_inf_hex18.h>
#include <libmesh/cell_inf_hex8.h>
#include <libmesh/cell_inf_prism12.h>
#include <libmesh/cell_inf_prism6.h>
#include <libmesh/cell_prism15.h>
#include <libmesh/cell_prism18.h>
#include <libmesh/cell_prism20.h>
#include <libmesh/cell_prism21.h>
#include <libmesh/cell_prism6.h>
#include <libmesh/cell_pyramid13.h>
#include <libmesh/cell_pyramid14.h>
#include <libmesh/cell_pyramid5.h>
#include <libmesh/cell_tet10.h>
#include <libmesh/cell_tet14.h>
#include <libmesh/cell_tet4.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/edge_edge3.h>
#include <libmesh/edge_edge4.h>
#include <libmesh/edge_inf_edge2.h>
#include <libmesh/face_inf_quad4.h>
#include <libmesh/face_inf_quad6.h>
#include <libmesh/face_quad4.h>
#include <libmesh/face_quad8.h>
#include <libmesh/face_quad9.h>
#include <libmesh/face_tri3.h>
#include <libmesh/face_tri6.h>
#include <libmesh/face_tri7.h>

#include <vector>

#include "libmesh_cppunit.h"

#define SIDETEST                                \
  CPPUNIT_TEST( testIsNodeOnSide );             \
  CPPUNIT_TEST( testNodesOnSide );              \
  CPPUNIT_TEST( testSidePtr );                  \
  CPPUNIT_TEST( testSidePtrFill );              \
  CPPUNIT_TEST( testBuildSidePtr );             \
  CPPUNIT_TEST( testBuildSidePtrFill );         \

using namespace libMesh;


template <typename ElemClass, ElemType side_type,
          unsigned short indexbegin, unsigned short indexend>
class SideTest : public CppUnit::TestCase {

private:
  ElemClass elem;
  std::vector<std::unique_ptr<Node>> nodes;

protected:
  std::string libmesh_suite_name;

public:
  void setUp() {
    elem.set_id() = 0;
#ifdef LIBMESH_ENABLE_AMR
    // Do tests with an Elem having a non-default p_level to ensure
    // that sides which are built have a matching p_level. p-refinement
    // is only available if LIBMESH_ENABLE_AMR is defined.
    elem.set_p_level(1);
#endif
    Point dummy;
    for (auto i : elem.node_index_range())
      {
        nodes.push_back(std::make_unique<Node>(dummy, /*id=*/i));
        elem.set_node(i, nodes[i].get());
      }
  }

  void tearDown() {}

  void testIsNodeOnSide()
  {
    LOG_UNIT_TEST;

    for (auto s : make_range(indexbegin, indexend))
      {
        std::unique_ptr<Elem> side = elem.build_side_ptr(s);
        for (auto n : elem.node_index_range())
          {
            const Node * node = elem.node_ptr(n);
            bool found_node = false;
            for (auto sn : side->node_index_range())
              if (node == side->node_ptr(sn))
                {
                  found_node = true;
                  break;
                }

            if (elem.is_node_on_side(n, s))
              {
                CPPUNIT_ASSERT(found_node);
              }
            else
              {
                CPPUNIT_ASSERT(!found_node);
              }
          }
      }
  }

  void testNodesOnSide()
  {
    LOG_UNIT_TEST;

    for (auto s : make_range(indexbegin, indexend))
      {
        std::unique_ptr<Elem> side = elem.build_side_ptr(s);
        std::vector<unsigned int> side_nodes = elem.nodes_on_side(s);

        CPPUNIT_ASSERT_EQUAL(side_nodes.size(), std::size_t(side->n_nodes()));

        for (auto sn : side->node_index_range())
          {
            const Node * node = side->node_ptr(sn);
            bool found_node = false;
            for (auto si : side_nodes)
              if (node == elem.node_ptr(si))
                {
                  found_node = true;
                  break;
                }
            CPPUNIT_ASSERT(found_node);
          }
      }
  }

  void testSidePtr()
  {
    LOG_UNIT_TEST;

    for (auto s : make_range(indexbegin, indexend))
      {
        std::unique_ptr<Elem> side = elem.side_ptr(s);

        CPPUNIT_ASSERT(side->type() ==
                       Elem::first_order_equivalent_type(side_type));
      }
  }

  void testSidePtrFill()
  {
    LOG_UNIT_TEST;

    std::unique_ptr<Elem> side;

    for (auto s : make_range(indexbegin, indexend))
      {
        elem.side_ptr(side, s);

        CPPUNIT_ASSERT(side->type() ==
                       Elem::first_order_equivalent_type(side_type));
      }
  }

  void testBuildSidePtr()
  {
    LOG_UNIT_TEST;

    for (auto s : make_range(indexbegin, indexend))
      {
        std::unique_ptr<Elem> side = elem.build_side_ptr(s);

        CPPUNIT_ASSERT(side->type() == side_type);
        CPPUNIT_ASSERT(side->subdomain_id() == elem.subdomain_id());

#ifdef LIBMESH_ENABLE_AMR
        // p-refinement is only available if LIBMESH_ENABLE_AMR is defined.
        CPPUNIT_ASSERT(side->p_level() == elem.p_level());
#endif
      }
  }

  void testBuildSidePtrFill()
  {
    LOG_UNIT_TEST;

    std::unique_ptr<Elem> side;

    for (auto s : make_range(indexbegin, indexend))
      {
        elem.build_side_ptr(side, s);
        std::unique_ptr<Elem> side_new = elem.build_side_ptr(s);

        CPPUNIT_ASSERT(side->type() == side_type);
        CPPUNIT_ASSERT(*side == *side_new);
      }
  }

};


#define INSTANTIATE_SIDETEST(elemclass, sidetype, indexbegin, indexend)                \
  class SideTest_##elemclass##_##sidetype##_##indexbegin##_##indexend :                \
    public SideTest<elemclass, sidetype, indexbegin, indexend> {                       \
  public:                                                                              \
  SideTest_##elemclass##_##sidetype##_##indexbegin##_##indexend() :                    \
    SideTest<elemclass,sidetype,indexbegin,indexend>() {                               \
    if (unitlog->summarized_logs_enabled())                                            \
      this->libmesh_suite_name = "SideTest";                                           \
    else                                                                               \
      this->libmesh_suite_name = "SideTest_" #elemclass"_" #sidetype "_" #indexbegin "_" #indexend; \
  }                                                                                    \
  CPPUNIT_TEST_SUITE( SideTest_##elemclass##_##sidetype##_##indexbegin##_##indexend ); \
  SIDETEST                                                                             \
  CPPUNIT_TEST_SUITE_END();                                                            \
  };                                                                                   \
                                                                                       \
  CPPUNIT_TEST_SUITE_REGISTRATION( SideTest_##elemclass##_##sidetype##_##indexbegin##_##indexend );

INSTANTIATE_SIDETEST(Hex20,     QUAD8, 0, 6);
INSTANTIATE_SIDETEST(Hex27,     QUAD9, 0, 6);
INSTANTIATE_SIDETEST(Hex8,      QUAD4, 0, 6);
INSTANTIATE_SIDETEST(Prism15,   TRI6,  0, 1);
INSTANTIATE_SIDETEST(Prism15,   QUAD8, 1, 4);
INSTANTIATE_SIDETEST(Prism15,   TRI6,  4, 5);
INSTANTIATE_SIDETEST(Prism18,   TRI6,  0, 1);
INSTANTIATE_SIDETEST(Prism18,   QUAD9, 1, 4);
INSTANTIATE_SIDETEST(Prism18,   TRI6,  4, 5);
INSTANTIATE_SIDETEST(Prism20,   TRI7,  0, 1);
INSTANTIATE_SIDETEST(Prism20,   QUAD9, 1, 4);
INSTANTIATE_SIDETEST(Prism20,   TRI7,  4, 5);
INSTANTIATE_SIDETEST(Prism21,   TRI7,  0, 1);
INSTANTIATE_SIDETEST(Prism21,   QUAD9, 1, 4);
INSTANTIATE_SIDETEST(Prism21,   TRI7,  4, 5);
INSTANTIATE_SIDETEST(Prism6,    TRI3,  0, 1);
INSTANTIATE_SIDETEST(Prism6,    QUAD4, 1, 4);
INSTANTIATE_SIDETEST(Prism6,    TRI3,  4, 5);
INSTANTIATE_SIDETEST(Pyramid13, TRI6,  0, 4);
INSTANTIATE_SIDETEST(Pyramid13, QUAD8, 4, 5);
INSTANTIATE_SIDETEST(Pyramid14, TRI6,  0, 4);
INSTANTIATE_SIDETEST(Pyramid14, QUAD9, 4, 5);
INSTANTIATE_SIDETEST(Pyramid5,  TRI3,  0, 4);
INSTANTIATE_SIDETEST(Pyramid5,  QUAD4, 4, 5);
INSTANTIATE_SIDETEST(Tet10,     TRI6,  0, 4);
INSTANTIATE_SIDETEST(Tet14,     TRI7,  0, 4);
INSTANTIATE_SIDETEST(Tet4,      TRI3,  0, 4);
INSTANTIATE_SIDETEST(Edge2, NODEELEM,  0, 2);
INSTANTIATE_SIDETEST(Edge3, NODEELEM,  0, 2);
INSTANTIATE_SIDETEST(Edge4, NODEELEM,  0, 2);
INSTANTIATE_SIDETEST(Quad4,     EDGE2, 0, 4);
INSTANTIATE_SIDETEST(Quad8,     EDGE3, 0, 4);
INSTANTIATE_SIDETEST(Quad9,     EDGE3, 0, 4);
INSTANTIATE_SIDETEST(Tri3,      EDGE2, 0, 3);
INSTANTIATE_SIDETEST(Tri6,      EDGE3, 0, 3);
INSTANTIATE_SIDETEST(Tri7,      EDGE3, 0, 3);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
INSTANTIATE_SIDETEST(InfHex16,   QUAD8,    0, 1);
INSTANTIATE_SIDETEST(InfHex16,   INFQUAD6, 1, 5);
INSTANTIATE_SIDETEST(InfHex18,   QUAD9,    0, 1);
INSTANTIATE_SIDETEST(InfHex18,   INFQUAD6, 1, 5);
INSTANTIATE_SIDETEST(InfHex8,    QUAD4,    0, 1);
INSTANTIATE_SIDETEST(InfHex8,    INFQUAD4, 1, 5);
INSTANTIATE_SIDETEST(InfPrism12, TRI6,     0, 1);
INSTANTIATE_SIDETEST(InfPrism12, INFQUAD6, 1, 4);
INSTANTIATE_SIDETEST(InfPrism6,  TRI3,     0, 1);
INSTANTIATE_SIDETEST(InfPrism6,  INFQUAD4, 1, 4);
INSTANTIATE_SIDETEST(InfEdge2,   NODEELEM, 0, 1);
INSTANTIATE_SIDETEST(InfQuad4,   EDGE2,    0, 1);
INSTANTIATE_SIDETEST(InfQuad4,   INFEDGE2, 1, 3);
INSTANTIATE_SIDETEST(InfQuad6,   EDGE3,    0, 1);
INSTANTIATE_SIDETEST(InfQuad6,   INFEDGE2, 1, 3);
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
