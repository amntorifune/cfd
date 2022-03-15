#include "topo.h"
#include "flag.h"
#include <math.h>

void sor(
    double   P[NX + 4][NY + 4][NZ + 4],
    int    FFP[NX + 4][NY + 4][NZ + 4][3],
    double BBP[NX + 4][NY + 4][NZ + 4][3],
    double DIV[NX + 4][NY + 4][NZ + 4],
    double   C[NX + 4][NY + 4][NZ + 4][6],
    double   OMEGA,
    double   DT,
    double  &R
) {
    double ERR = 0;
    int    CNT = 0;
}