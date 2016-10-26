
#include "fe_test.h"

INSTANTIATE_FETEST(FIRST, L2_LAGRANGE, EDGE2);
INSTANTIATE_FETEST(FIRST, L2_LAGRANGE, TRI3);
INSTANTIATE_FETEST(FIRST, L2_LAGRANGE, QUAD4);
INSTANTIATE_FETEST(FIRST, L2_LAGRANGE, TET4);
INSTANTIATE_FETEST(FIRST, L2_LAGRANGE, HEX8);
INSTANTIATE_FETEST(FIRST, L2_LAGRANGE, PRISM6);

INSTANTIATE_FETEST(SECOND, L2_LAGRANGE, EDGE3);
INSTANTIATE_FETEST(SECOND, L2_LAGRANGE, TRI6);
INSTANTIATE_FETEST(SECOND, L2_LAGRANGE, QUAD9);
INSTANTIATE_FETEST(SECOND, L2_LAGRANGE, TET10);
INSTANTIATE_FETEST(SECOND, L2_LAGRANGE, HEX27);
INSTANTIATE_FETEST(SECOND, L2_LAGRANGE, PRISM18);