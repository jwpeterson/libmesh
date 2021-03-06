// libmesh includes
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh.h>
#include <libmesh/reference_elem.h>
#include <libmesh/node.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/tensor_value.h>

// unit test includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class VolumeTest : public CppUnit::TestCase
{

public:
  CPPUNIT_TEST_SUITE( VolumeTest );
  CPPUNIT_TEST( testEdge3Volume );
  CPPUNIT_TEST( testEdge3Invertible );
  CPPUNIT_TEST( testEdge4Invertible );
  CPPUNIT_TEST( testQuad4Invertible );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void testEdge3Volume()
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_line (mesh, /*nelem=*/1, 0., 1., EDGE3);
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(3), mesh.n_nodes());

    auto edge3 = mesh.query_elem_ptr(0);
    if (!edge3) // We may be on a distributed mesh
      return;

    // Check unperturbed, straight edge case
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Get references to the individual Edge3 nodes
    auto & middle_node = edge3->node_ref(2);
    auto & right_node = edge3->node_ref(1);

    // Check middle node perturbed in +x direction case. This should
    // not change the volume because it's still a straight line
    // element.
    middle_node = Point(0.5 + 1.e-3, 0., 0.);
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Check middle node perturbed in -x direction case. This should
    // not change the volume because it's still a straight line
    // element.
    middle_node = Point(0.5 - 1.e-3, 0., 0.);
    LIBMESH_ASSERT_FP_EQUAL(1.0, edge3->volume(), TOLERANCE*TOLERANCE);

    // Check volume of actual curved element against pre-computed value.
    middle_node = Point(0.5, 0.25, 0.);
    right_node = Point(1., 1., 0.);
    LIBMESH_ASSERT_FP_EQUAL(1.4789428575446, edge3->volume(), TOLERANCE);

    // Compare with volume computed by base class Elem::volume() call
    // which uses quadrature.  We don't expect this to have full
    // floating point accuracy.
    middle_node = Point(0.5, 0.1, 0.);
    right_node = Point(1., 0., 0.);
    LIBMESH_ASSERT_FP_EQUAL(edge3->Elem::volume(), edge3->volume(), std::sqrt(TOLERANCE));
  }

  void testEdge3Invertible()
  {
    // 1.) This is the original test which started the investigation
    // of determining invertibility.  In this test, the actual
    // midpoint of nodes 0 and 1 is 0.5*(1.100328e2 + 1.176528e2) =
    // 113.8428, so we can see that the middle node is closer to the
    // left endpoint. In this case, it is too close and the element is
    // not invertible.
    bool invertible = test_elem({
      Point(-3.566160e1, -6.690970e-1, 1.100328e2),
      Point(-3.566160e1, -6.690970e-1, 1.176528e2),
      Point(-3.566160e1, -6.690970e-1, 1.115568e2)}, EDGE3);
    CPPUNIT_ASSERT(!invertible);

    // 2.) Just like case 1, but now node 2 is at the midpoint, so
    // this case is invertible.
    invertible = test_elem({
      Point(-3.566160e1, -6.690970e-1, 1.100328e2),
      Point(-3.566160e1, -6.690970e-1, 1.176528e2),
      Point(-3.566160e1, -6.690970e-1, 113.8428)}, EDGE3);
    CPPUNIT_ASSERT(invertible);

    // 3.) Non-collinear case where the mid-edge node is "above" and "way
    // past" the right endpoint. This case is not invertible
    invertible = test_elem({Point(0, 0, 0), Point(1, 0, 0), Point(3.5, 1.5, 0)}, EDGE3);
    CPPUNIT_ASSERT(!invertible);
  }

  void testEdge4Invertible()
  {
    // Reference Elem should be invertible
    {
      const Elem & edge4 = ReferenceElem::get(EDGE4);
      CPPUNIT_ASSERT(edge4.has_invertible_map());
    }

    // If node 2 goes to the left past -5/9 = -.555, the element becomes non-invertible
    {
      // x2 > -5/9, the map is still invertible
      bool invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(-0.5, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(invertible);

      // x2 < -5/9, it is too close to x0 now
      invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(-0.57, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(!invertible);
    }

    // If node 2 goes to the right past 5/21 ~ 0.2381, the element becomes non-invertible
    {
      // x2 < 5/21, the map should still be invertible
      bool invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(Real(3)/21, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(invertible);

      // x2 > 5/21, x2 is too close to x3 now
      invertible =
        test_elem({Point(-1, 0, 0), Point(1, 0, 0), Point(Real(6)/21, 0, 0), Point(Real(1)/3, 0, 0)},
                  EDGE4);
      CPPUNIT_ASSERT(!invertible);
    }
  }

  void testQuad4Invertible()
  {
    // Case 1: Test that rigid body rotations have no effect on the
    // invertibility of the reference element
    {
      // 1a) The reference element rotated into various different different planes.
      std::vector<Point> pts = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
      bool invertible = test_elem(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1b) Rotate all points about x-axis by 90 degrees
      Real cost = std::cos(.5*libMesh::pi);
      Real sint = std::sin(.5*libMesh::pi);
      RealTensorValue Rx(1, 0, 0,
                         0, cost, sint,
                         0, sint, cost);

      for (auto & pt : pts)
        pt = Rx * pt;

      invertible = test_elem(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1c) Rotate all points about z-axis by 90 degrees
      RealTensorValue Rz(cost, -sint, 0,
                         sint,  cost, 0,
                         0,        0, 1);

      for (auto & pt : pts)
        pt = Rz * pt;

      invertible = test_elem(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);

      // 1d) Rotate all points about y-axis by 270 degrees
      RealTensorValue Ry(cost,  0, sint,
                         0,     1, 0,
                         -sint, 0, cost);

      for (int cnt=0; cnt<3; ++cnt)
        for (auto & pt : pts)
          pt = Ry * pt;

      invertible = test_elem(pts, QUAD4);
      CPPUNIT_ASSERT(invertible);
    }

    // Case 2: Planar quad with top right vertex displaced to the position
    // (alpha, alpha). Some different cases are described below.
    // .) alpha==1: affine case, always invertible
    // .) 1/2 < alpha < 1: planar case, invertible
    // .) alpha<=1/2: planar case but node is now at center of the
    //    element, should give a zero/negative Jacobian on the displaced
    //    Node -> not invertible.
    {
      const Real alpha = .5;

      bool invertible =
        test_elem({Point(0, 0, 0), Point(1, 0, 0), Point(alpha, alpha, 0), Point(0, 1, 0)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }

    // Case 3) Top right corner is moved to (alpha, 1, 0). Element
    // becomes non-invertible when alpha < 0.
    {
      const Real alpha = -0.25;

      bool invertible =
        test_elem({Point(0, 0, 0), Point(1, 0, 0), Point(alpha, 1, 0), Point(0, 1, 0)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }

    // Case 4) Degenerate case - all 4 points at same location. This
    // zero-volume element does not have an invertible map.
    {
      const Real alpha = std::log(2);

      bool invertible =
        test_elem({Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha),
                   Point(alpha, alpha, alpha)}, QUAD4);

      CPPUNIT_ASSERT(!invertible);
    }
  }

protected:

  // Helper function that builds the specified type of Elem from a
  // vector of Points and returns the value of has_invertible_map()
  // for that Elem.
  bool test_elem(const std::vector<Point> & pts,
                 ElemType elem_type)
  {
    const unsigned int n_points = pts.size();

    // Create Nodes
    std::vector<std::unique_ptr<Node>> nodes(n_points);
    for (unsigned int i=0; i<n_points; i++)
      nodes[i] = Node::build(pts[i], /*id*/ i);

    // Create Elem, assign nodes
    std::unique_ptr<Elem> elem = Elem::build(elem_type, /*parent*/ nullptr);

    // Make sure we were passed consistent input to build this type of Elem
    libmesh_error_msg_if(elem->n_nodes() != n_points,
                         "Wrong number of points "
                         << n_points
                         << " provided to build a "
                         << Utility::enum_to_string(elem_type));

    for (unsigned int i=0; i<n_points; i++)
      elem->set_node(i) = nodes[i].get();

    // Return whether or not this Elem has an invertible map
    return elem->has_invertible_map();
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( VolumeTest );
