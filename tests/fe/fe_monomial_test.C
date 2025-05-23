
#include "fe_test.h"

INSTANTIATE_FETEST(FIRST, MONOMIAL, EDGE2);
INSTANTIATE_FETEST(SECOND, MONOMIAL, EDGE3);
INSTANTIATE_FETEST(THIRD, MONOMIAL, EDGE3);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, EDGE3);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, EDGE3);

#if LIBMESH_DIM > 1
INSTANTIATE_FETEST(FIRST, MONOMIAL, TRI3);
INSTANTIATE_FETEST(SECOND, MONOMIAL, TRI6);
INSTANTIATE_FETEST(THIRD, MONOMIAL, TRI6);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, TRI6);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, TRI6);

INSTANTIATE_FETEST(SECOND, MONOMIAL, TRI7);
INSTANTIATE_FETEST(THIRD, MONOMIAL, TRI7);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, TRI7);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, TRI7);

INSTANTIATE_FETEST(FIRST, MONOMIAL, QUAD4);
INSTANTIATE_FETEST(SECOND, MONOMIAL, QUAD8);
INSTANTIATE_FETEST(SECOND, MONOMIAL, QUAD9);
INSTANTIATE_FETEST(THIRD, MONOMIAL, QUAD8);
INSTANTIATE_FETEST(THIRD, MONOMIAL, QUAD9);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, QUAD8);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, QUAD9);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, QUAD8);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, QUAD9);

// This only is an exact projection because we're using a test
// C0Polygon with an affine map
INSTANTIATE_FETEST(FIRST, MONOMIAL, C0POLYGON);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FETEST(FIRST, MONOMIAL, TET4);
INSTANTIATE_FETEST(SECOND, MONOMIAL, TET10);
INSTANTIATE_FETEST(THIRD, MONOMIAL, TET10);
INSTANTIATE_FETEST(THIRD, MONOMIAL, TET14);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, TET14);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, TET14);

INSTANTIATE_FETEST(FIRST, MONOMIAL, HEX8);
INSTANTIATE_FETEST(SECOND, MONOMIAL, HEX20);
INSTANTIATE_FETEST(SECOND, MONOMIAL, HEX27);
INSTANTIATE_FETEST(THIRD, MONOMIAL, HEX20);
INSTANTIATE_FETEST(THIRD, MONOMIAL, HEX27);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, HEX20);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, HEX27);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, HEX20);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, HEX27);

INSTANTIATE_FETEST(FIRST, MONOMIAL, PRISM6);
INSTANTIATE_FETEST(SECOND, MONOMIAL, PRISM15);
INSTANTIATE_FETEST(SECOND, MONOMIAL, PRISM18);
INSTANTIATE_FETEST(THIRD, MONOMIAL, PRISM15);
INSTANTIATE_FETEST(THIRD, MONOMIAL, PRISM18);
INSTANTIATE_FETEST(THIRD, MONOMIAL, PRISM20);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, PRISM15);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, PRISM18);
INSTANTIATE_FETEST(FOURTH, MONOMIAL, PRISM21);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, PRISM15);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, PRISM18);
INSTANTIATE_FETEST(FIFTH, MONOMIAL, PRISM21);

INSTANTIATE_FETEST(FIRST, MONOMIAL, C0POLYHEDRON);
#endif
