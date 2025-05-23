#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/face_quad4.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/dof_map.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh_refinement.h>

#include <libmesh/discontinuity_measure.h>
#include <libmesh/error_vector.h>
#include <libmesh/overlap_coupling.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class SlitFunc : public FEMFunctionBase<Number>
{
public:

  SlitFunc() {}

  ~SlitFunc () {}

  virtual void init_context (const FEMContext &) override {}

  virtual std::unique_ptr<FEMFunctionBase<Number>>
  clone () const override
  {
    return std::make_unique<SlitFunc>();
  }

  virtual Number operator() (const FEMContext & c,
                             const Point & p,
                             const Real /*time*/ = 0.) override
  {
    using std::abs;

    const Real & x = p(0);
    const Real & y = p(1);
    const Point centroid = c.get_elem().vertex_average();
    const Real sign = centroid(1)/std::abs(centroid(1));

    // For testing we want something discontinuous on the slit,
    // continuous everywhere else, and bilinear on all coarse quads
    return (abs(x) + abs(2-x) - 2*abs(1-x)) * (1-abs(y)) * sign;
  }

  virtual void operator() (const FEMContext & c,
                           const Point & p,
                           const Real time,
                           DenseVector<Number> & output) override
  {
    for (unsigned int i=0; i != output.size(); ++i)
      output(i) = (*this)(c, p, time);
  }
};







class SlitMeshTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that a 2D mesh with nodes overlapping
   * on opposite sides of an internal, "slit" edge is usable.  The
   * mesh has to be connected at more than one node on each side of
   * the slit, however, to ensure that we can find point neighbors of
   * each node.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( SlitMeshTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

  std::unique_ptr<Mesh> _mesh;

  void build_mesh()
  {
    _mesh = std::make_unique<Mesh>(*TestCommWorld);

    // (-1,1)     (0,1)      (1,1)      (2,1)      (3,1)
    // o----------o----------o----------o----------o
    // |          |          |          |          |
    // |          |          |          |          |
    // |          |          |          |          |
    // |          |          |          |          |
    // o----------o==========8==========o----------o
    // (-1,0)     (0,0)      (1,0)      (2,0)      (3,0)
    // |          |          |          |          |
    // |          |          |          |          |
    // |          |          |          |          |
    // o----------o----------o----------o----------o
    // (-1,-1)    (0,-1)     (1,-1)     (2,-1)     (3,-1)

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0, 0.0), 0 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(1.0, 0.0), 6 ); // Doubled!
    _mesh->add_point( Point(2.0, 0.0), 7 );
    _mesh->add_point( Point(2.0, 1.0), 8 );
    _mesh->add_point( Point(2.0,-1.0), 9 );
    _mesh->add_point( Point(-1.0,-1.0), 10);
    _mesh->add_point( Point(-1.0, 0.0), 11);
    _mesh->add_point( Point(-1.0, 1.0), 12);
    _mesh->add_point( Point(3.0,-1.0), 13);
    _mesh->add_point( Point(3.0, 0.0), 14);
    _mesh->add_point( Point(3.0, 1.0), 15);

    {
      Elem * elem_top_left = _mesh->add_elem(Elem::build_with_id(QUAD4, 0));
      elem_top_left->set_node(0, _mesh->node_ptr(0));
      elem_top_left->set_node(1, _mesh->node_ptr(1));
      elem_top_left->set_node(2, _mesh->node_ptr(2));
      elem_top_left->set_node(3, _mesh->node_ptr(3));

      Elem * elem_bottom_left = _mesh->add_elem(Elem::build_with_id(QUAD4, 1));
      elem_bottom_left->set_node(0, _mesh->node_ptr(4));
      elem_bottom_left->set_node(1, _mesh->node_ptr(5));
      elem_bottom_left->set_node(2, _mesh->node_ptr(6));
      elem_bottom_left->set_node(3, _mesh->node_ptr(0));

      Elem * elem_top_right = _mesh->add_elem(Elem::build_with_id(QUAD4, 2));
      elem_top_right->set_node(0, _mesh->node_ptr(1));
      elem_top_right->set_node(1, _mesh->node_ptr(7));
      elem_top_right->set_node(2, _mesh->node_ptr(8));
      elem_top_right->set_node(3, _mesh->node_ptr(2));

      Elem * elem_bottom_right = _mesh->add_elem(Elem::build_with_id(QUAD4, 3));
      elem_bottom_right->set_node(0, _mesh->node_ptr(5));
      elem_bottom_right->set_node(1, _mesh->node_ptr(9));
      elem_bottom_right->set_node(2, _mesh->node_ptr(7));
      elem_bottom_right->set_node(3, _mesh->node_ptr(6));

      Elem * elem_top_leftleft = _mesh->add_elem(Elem::build_with_id(QUAD4, 4));
      elem_top_leftleft->set_node(0, _mesh->node_ptr(11));
      elem_top_leftleft->set_node(1, _mesh->node_ptr(0));
      elem_top_leftleft->set_node(2, _mesh->node_ptr(3));
      elem_top_leftleft->set_node(3, _mesh->node_ptr(12));

      Elem * elem_bottom_leftleft = _mesh->add_elem(Elem::build_with_id(QUAD4, 5));
      elem_bottom_leftleft->set_node(0, _mesh->node_ptr(10));
      elem_bottom_leftleft->set_node(1, _mesh->node_ptr(4));
      elem_bottom_leftleft->set_node(2, _mesh->node_ptr(0));
      elem_bottom_leftleft->set_node(3, _mesh->node_ptr(11));

      Elem * elem_top_rightright = _mesh->add_elem(Elem::build_with_id(QUAD4, 6));
      elem_top_rightright->set_node(0, _mesh->node_ptr(7));
      elem_top_rightright->set_node(1, _mesh->node_ptr(14));
      elem_top_rightright->set_node(2, _mesh->node_ptr(15));
      elem_top_rightright->set_node(3, _mesh->node_ptr(8));

      Elem * elem_bottom_rightright = _mesh->add_elem(Elem::build_with_id(QUAD4, 7));
      elem_bottom_rightright->set_node(0, _mesh->node_ptr(9));
      elem_bottom_rightright->set_node(1, _mesh->node_ptr(13));
      elem_bottom_rightright->set_node(2, _mesh->node_ptr(14));
      elem_bottom_rightright->set_node(3, _mesh->node_ptr(7));
    }

    // libMesh shouldn't renumber, or our based-on-initial-id
    // assertions later may fail.
    _mesh->allow_renumbering(false);

    _mesh->prepare_for_use();
  }

public:
  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();
#endif
  }

  void tearDown() {}

  void testMesh()
  {
    LOG_UNIT_TEST;

    // There'd better be 8 elements
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(8), _mesh->n_elem());

    // There'd better still be a full 16 nodes
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(16), _mesh->n_nodes());

    /* The middle nodes should still be distinct between the top and
     * bottom elements */
    if (_mesh->query_elem_ptr(0) && _mesh->query_elem_ptr(1))
      CPPUNIT_ASSERT( _mesh->elem_ref(0).node_id(1) != _mesh->elem_ref(1).node_id(2) );
    if (_mesh->query_elem_ptr(2) && _mesh->query_elem_ptr(3))
      CPPUNIT_ASSERT( _mesh->elem_ref(2).node_id(0) != _mesh->elem_ref(3).node_id(3) );

    /* The middle nodes should still be shared between left and right
     * elements on top and bottom */
    if (_mesh->query_elem_ptr(0) && _mesh->query_elem_ptr(2))
      CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(0).node_id(1),
                            _mesh->elem_ref(2).node_id(0) );
    if (_mesh->query_elem_ptr(1) && _mesh->query_elem_ptr(3))
      CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(1).node_id(2),
                            _mesh->elem_ref(3).node_id(3) );
  }

};

class SlitMeshRefinedMeshTest : public SlitMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we do a
   * uniform refinement and make sure the result mesh is consistent. i.e.
   * the new node shared between the 1D elements is the same as the node
   * shared on the underlying quads, and so on.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( SlitMeshRefinedMeshTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );
#endif

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement(*_mesh).uniformly_refine(1);
#endif
#endif
  }

  void testMesh()
  {
    LOG_UNIT_TEST;

#ifdef LIBMESH_ENABLE_AMR
    // We should have 40 total and 32 active elements.
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(40), _mesh->n_elem());
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(32), _mesh->n_active_elem());

    // We should have 48 nodes, not 45 or 46
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(48), _mesh->n_nodes());
#endif
  }
};

class SlitMeshRefinedSystemTest : public SlitMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we
   * create a system and set dof values to make sure they are properly
   * interpolated after refinement.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( SlitMeshRefinedSystemTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST( testSystem );

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  CPPUNIT_TEST( testRestart );
#endif
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

  System* _sys;
  std::unique_ptr<EquationSystems> _es;

public:

  void setUp()
  {
#if LIBMESH_DIM > 1
    this->build_mesh();

    // libMesh *should* renumber now, or a DistributedMesh might not
    // have contiguous ids, which is a requirement to write xda files.
    _mesh->allow_renumbering(true);

    _es = std::make_unique<EquationSystems>(*_mesh);
    _sys = &_es->add_system<System> ("SimpleSystem");
    _sys->add_variable("u", FIRST);

    // We're going to be integrating across the slit in the mesh, so
    // let's make sure we can *see* elements and data across the slit.
    _mesh->allgather();
    _sys->get_dof_map().add_algebraic_ghosting_functor
      (std::make_shared<OverlapCoupling>());
    _mesh->delete_remote_elements();

    _es->init();
    SlitFunc slitfunc;
    _sys->project_solution(&slitfunc);

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement(*_mesh).uniformly_refine(1);
    _es->reinit();
    MeshRefinement(*_mesh).uniformly_refine(1);
    _es->reinit();
#endif
#endif
  }

  void tearDown() {}

  void testMesh()
  {
    LOG_UNIT_TEST;

#ifdef LIBMESH_ENABLE_AMR
    // We should have 168 total and 128 active elements.
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(8+32+128), _mesh->n_elem());
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(128), _mesh->n_active_elem());

    // We should have 160 nodes
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(160), _mesh->n_nodes());
#endif
  }

  void testSystem()
  {
    LOG_UNIT_TEST;

    SlitFunc slitfunc;

    unsigned int dim = 2;

    CPPUNIT_ASSERT_EQUAL( _sys->n_vars(), 1u );

    FEMContext context(*_sys);
    FEBase * fe = NULL;
    context.get_element_fe( 0, fe, dim );
    const std::vector<Point> & xyz = fe->get_xyz();
    fe->get_phi();

    for (const auto & elem : _mesh->active_local_element_ptr_range())
      {
        context.pre_fe_reinit(*_sys, elem);
        context.elem_fe_reinit();

        const unsigned int n_qp = xyz.size();

        for (unsigned int qp=0; qp != n_qp; ++qp)
          {
            const Number exact_val = slitfunc(context, xyz[qp]);

            const Number discrete_val = context.interior_value(0, qp);

            LIBMESH_ASSERT_NUMBERS_EQUAL
              (exact_val, discrete_val, TOLERANCE*TOLERANCE);
          }
      }

    // We should have no discontinuities (beyond floating-point error)
    // between topologically connected elements
    DiscontinuityMeasure connected_dm;
    ErrorVector connected_err;
    connected_dm.estimate_error(*_sys, connected_err);
    const Real mean_connected_disc = connected_err.mean();
    CPPUNIT_ASSERT_LESS(Real(1e-14), mean_connected_disc);

    // We should be able to see the discontinuity along the slit
    DiscontinuityMeasure slit_dm;
    slit_dm.integrate_slits = true;
    ErrorVector slit_disc;
    slit_dm.estimate_error(*_sys, slit_disc);
    const Real mean_slit_disc = slit_disc.mean();
    CPPUNIT_ASSERT_GREATER(Real(1e-3), mean_slit_disc);
  }

  void testRestart()
  {
    LOG_UNIT_TEST;

    SlitFunc slitfunc;

    _mesh->write("slit_mesh.xda");
    _es->write("slit_solution.xda",
               EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_SERIAL_FILES);

    Mesh mesh2(*TestCommWorld);
    mesh2.read("slit_mesh.xda");
    EquationSystems es2(mesh2);
    es2.read("slit_solution.xda");

    System & sys2 = es2.get_system<System> ("SimpleSystem");

    unsigned int dim = 2;

    CPPUNIT_ASSERT_EQUAL( sys2.n_vars(), 1u );

    FEMContext context(sys2);
    FEBase * fe = NULL;
    context.get_element_fe( 0, fe, dim );
    const std::vector<Point> & xyz = fe->get_xyz();
    fe->get_phi();

    // While we're in the middle of a unique id based test case, let's
    // make sure our unique ids were all read in correctly too.
    std::unique_ptr<PointLocatorBase> locator = _mesh->sub_point_locator();

    if (!_mesh->is_serial())
      locator->enable_out_of_mesh_mode();

    for (const auto & elem : mesh2.active_local_element_ptr_range())
      {
        const Elem * mesh1_elem = (*locator)(elem->vertex_average());
        if (mesh1_elem)
          {
            CPPUNIT_ASSERT_EQUAL( elem->unique_id(),
                                  mesh1_elem->unique_id() );

            for (unsigned int n=0; n != elem->n_nodes(); ++n)
              {
                const Node & node       = elem->node_ref(n);
                const Node & mesh1_node = mesh1_elem->node_ref(n);
                CPPUNIT_ASSERT_EQUAL( node.unique_id(),
                                      mesh1_node.unique_id() );
              }
          }

        context.pre_fe_reinit(sys2, elem);
        context.elem_fe_reinit();

        const unsigned int n_qp = xyz.size();

        for (unsigned int qp=0; qp != n_qp; ++qp)
          {
            const Number exact_val = slitfunc(context, xyz[qp]);

            const Number discrete_val = context.interior_value(0, qp);

            LIBMESH_ASSERT_NUMBERS_EQUAL
              (exact_val, discrete_val, TOLERANCE*TOLERANCE);
          }
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshRefinedMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshRefinedSystemTest );
