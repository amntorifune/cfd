Fractional step solver for lid-driven cavity flow.
# ns2d
Central difference, only works for moderate Re (no more than 3500 on a 256x256 mesh).
# ns2d-up
Upwind, works for all Reynolds numbers, but only of 1st-order accuracy.
# ns2dm
Almost the same as ns2d, but monitors some important values such as div(U) and how many iterations it takes for the pressure equation to converge.
# ns2dcc
Uses 'collocated at cell center' mesh, not bad in accuracy but VERY unstable.
# ns2dcv
Uses 'collocated at vertex' mesh, a bad example because it cannot give the right result.
# ns2dn
Uses non-conservative forms of the equations, although for incompressible fluids the non-conservative form should be equivalent to the conservative forms, the result turned out to be different, maybe due to the numerical scheme.
