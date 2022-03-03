#ifndef FUNCLIST_H
#define FUNCLIST_H

#include "topo.h"

extern void geom(
    double   X[NX + 2][NY + 2][NZ + 2][3]
);

extern void init(
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double  UD[NX + 2][NY + 2][NZ + 2][3],
    double  UC[NX + 2][NY + 2][NZ + 2][3],
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double   P[NX + 2][NY + 2][NZ + 2],
    double SGS[NX + 2][NY + 2][NZ + 2],
    double   X[NX + 2][NY + 2][NZ + 2][3],
    double  KX[NX + 2][NY + 2][NZ + 2][3],
    double   J[NX + 2][NY + 2][NZ + 2],
    double   C[NX + 2][NY + 2][NZ + 2][6],
    double   G[NX + 2][NY + 2][NZ + 2][3]
);

extern void bcu(
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double  UD[NX + 2][NY + 2][NZ + 2][3],
    double  KX[NX + 2][NY + 2][NZ + 2][3],
    double   DT
);

extern void bcp(
    double   P[NX + 2][NY + 2][NZ + 2],
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double   C[NX + 2][NY + 2][NZ + 2][6],
    double  KX[NX + 2][NY + 2][NZ + 2][3],
    double   RI
);

extern void solid(
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double  UU[NX + 2][NY + 2][NZ + 2][3]
);

extern void fs1(
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double  UD[NX + 2][NY + 2][NZ + 2][3],
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double SGS[NX + 2][NY + 2][NZ + 2],
    double  KX[NX + 2][NY + 2][NZ + 2][3],
    double   J[NX + 2][NY + 2][NZ + 2],
    double   C[NX + 2][NY + 2][NZ + 2][6],
    double   ALPHA,
    double   DT,
    double   RI
);

extern void contra(
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double  UC[NX + 2][NY + 2][NZ + 2][3],
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double  KX[NX + 2][NY + 2][NZ + 2][3],
    double   J[NX + 2][NY + 2][NZ + 2]
);

extern void diver(
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double   J[NX + 2][NY + 2][NZ + 2],
    double DIV[NX + 2][NY + 2][NZ + 2],
    double  &AD
);

extern void sor(
    double   P[NX + 2][NY + 2][NZ + 2],
    double DIV[NX + 2][NY + 2][NZ + 2],
    double   C[NX + 2][NY + 2][NZ + 2][6],
    double   OMEGA,
    double   E,
    int      MAXIT,
    double   DT,
    int     &IT,
    double  &R
);

extern void avp0(
    double   P[NX + 2][NY + 2][NZ + 2]
);

extern void fs2 (
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double   P[NX + 2][NY + 2][NZ + 2],
    double  KX[NX + 2][NY + 2][NZ + 2][3],
    double   G[NX + 2][NY + 2][NZ + 2][3],
    double   DT
);

extern void fio(
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double   P[NX + 2][NY + 2][NZ + 2],
    double   X[NX + 2][NY + 2][NZ + 2][3],
    char*    fname
);

#endif