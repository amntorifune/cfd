#ifndef VARLIST_H
#define VARLIST_H

#include "topo.h"

double   X[NX + 2][NY + 2][NZ + 2][3];
double   U[NX + 2][NY + 2][NZ + 2][3];
double  UD[NX + 2][NY + 2][NZ + 2][3];
double  UC[NX + 2][NY + 2][NZ + 2][3];
double  UU[NX + 2][NY + 2][NZ + 2][3];
double   P[NX + 2][NY + 2][NZ + 2];
double DIV[NX + 2][NY + 2][NZ + 2];
double  KX[NX + 2][NY + 2][NZ + 2][3];
double   J[NX + 2][NY + 2][NZ + 2];
double   G[NX + 2][NY + 2][NZ + 2][3];
double   C[NX + 2][NY + 2][NZ + 2][6];
double SGS[NX + 2][NY + 2][NZ + 2];

const double   DT    = 0.001;
const double   RE    = 10000;
const double   RI    = 1 / RE;
const double   EPOI  = 1E-3;
const double   EDIV  = 1E-2;
const double   OMEGA = 1.2;
const double   ALPHA = 0.5;
const int      MAXIT = 100;
const int      NSTEP = 50000;
double         AD    = 0;
double         R     = 0;
int            IT    = 0;
int            STEP  = 0;

#endif
