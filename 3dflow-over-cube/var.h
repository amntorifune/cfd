#ifndef VAR_H
#define VAR_H

#include "topo.h"

double   U[NX + 4][NY + 4][NZ + 4][3];
double  UA[NX + 4][NY + 4][NZ + 4][3];
double  UC[NX + 4][NY + 4][NZ + 4][3];
double  UU[NX + 4][NY + 4][NZ + 4][3];
double UUA[NX + 4][NY + 4][NZ + 4][3];
int    FFU[NX + 4][NY + 4][NZ + 4][3];
double BBU[NX + 4][NY + 4][NZ + 4][3][3];
double   P[NX + 4][NY + 4][NZ + 4];
int    FFP[NX + 4][NY + 4][NZ + 4][3];
double BBP[NX + 4][NY + 4][NZ + 4][3];
double DIV[NX + 4][NY + 4][NZ + 4];
double   X[NX + 4][NY + 4][NZ + 4][3];
double  KX[NX + 4][NY + 4][NZ + 4][3];
double   J[NX + 4][NY + 4][NZ + 4];
double   G[NX + 4][NY + 4][NZ + 4][3];
double   C[NX + 4][NY + 4][NZ + 4][6];
double SGS[NX + 4][NY + 4][NZ + 4];
int      F[NX + 4][NY + 4][NZ + 4];

#endif
