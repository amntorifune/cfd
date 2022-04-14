#ifndef FUNC_H
#define FUNC_H

#include "topo.h"

extern void p0gradient(
    double   P[NX + 4][NY + 4],
    double BBP[NX + 4][NY + 4][2]
);

extern void inflow(
    double BBU[NX + 4][NY + 4][2][2],
    double t
);

extern void outflow(
    double   U[NX + 4][NY + 4][2],
    double  UU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double   X[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4],
    double   DT
);

extern void slip(
    double   U[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2]
);

extern void contra(
    double   U[NX + 4][NY + 4][2],
    double  UC[NX + 4][NY + 4][2],
    double  UU[NX + 4][NY + 4][2],
    int    FFU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double   X[NX + 4][NY + 4][2],
    double  KX[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4]
);

extern void diver(
    double  UU[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4],
    double DIV[NX + 4][NY + 4],
    int      F[NX + 4][NY + 4],
    double  &AD
);

extern void sor(
    double   P[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double PSI[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   OMEGA,
    double   DT,
    double  &R
);

extern void jacobi(
    double   P[NX + 4][NY + 4],
    double  PD[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double PSI[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   OMEGA,
    double   DT,
    double  &R
);

extern void hsmac(
    double   P[NX + 4][NY + 4],
    double PSI[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   DT
);

extern void pseudo(
    double   U[NX + 4][NY + 4][2],
    double  UA[NX + 4][NY + 4][2],
    double  UU[NX + 4][NY + 4][2],
    int    FFU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double SGS[NX + 4][NY + 4],
    double  KX[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   ALPHA,
    double   DT,
    double   RI
);

extern void correction1(
    double   U[NX + 4][NY + 4][2],
    double  UD[NX + 4][NY + 4][2],
    double   P[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double  KX[NX + 4][NY + 4][2],
    int      F[NX + 4][NY + 4],
    double   DT
);

extern void correction2(
    double  UU[NX + 4][NY + 4][2],
    double UUD[NX + 4][NY + 4][2],
    int    FFU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double   P[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double   X[NX + 4][NY + 4][2],
    double  KX[NX + 4][NY + 4][2],
    double   G[NX + 4][NY + 4][2],
    double   DT
);

extern void correction3(
    double  P[NX + 4][NY + 4],
    double Pi[NX + 4][NY + 4]
);

#endif
