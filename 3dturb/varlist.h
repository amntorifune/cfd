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

double   DT    = 0.001;
double   RE    = 10000;
double   RI    = 1 / RE;
double   EPOI  = 1E-3;
double   DIV   = 0;
double   EDIV  = 1E-5;
double   OMEGA = 1.2;
double   ALPHA = 0.5;
double   AD    = 1E9;
double   R     = 0;
int      MAXIT = 500;
int      IT    = 0;
int      NSTEP = 50000;
int      STEP  = 0;

#endif