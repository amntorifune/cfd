#ifndef VAR_H
#define VAR_H

#include "topo.h"

//  physical variables

double   U[NX + 4][NY + 4][2];
double  UA[NX + 4][NY + 4][2];
double  UP[NX + 4][NY + 4][2];
double  UD[NX + 4][NY + 4][2];
double  UC[NX + 4][NY + 4][2];
double  UU[NX + 4][NY + 4][2];
double UUA[NX + 4][NY + 4][2];
double UUP[NX + 4][NY + 4][2];
double UUD[NX + 4][NY + 4][2];
int    FFU[NX + 4][NY + 4][2];
double BBU[NX + 4][NY + 4][2][2];
double   P[NX + 4][NY + 4];
double  PD[NX + 4][NY + 4];
int    FFP[NX + 4][NY + 4][2];
double BBP[NX + 4][NY + 4][2];
double  Pi[NX + 4][NY + 4];
double PiD[NX + 4][NY + 4];
double BPi[NX + 4][NY + 4][2];
double DVR[NX + 4][NY + 4];
double DVA[NX + 4][NY + 4];
double DVP[NX + 4][NY + 4];
double   X[NX + 4][NY + 4][2];
double  KX[NX + 4][NY + 4][2];
double   J[NX + 4][NY + 4];
double   G[NX + 4][NY + 4][2];
double   C[NX + 4][NY + 4][4];
double SGS[NX + 4][NY + 4];
int      F[NX + 4][NY + 4];

//  parameters

const double   DT    = 0.005;
const double   RE    = 1000;
const double   RI    = 1 / RE;
const double   EPOI  = 1E-6;
const double   EDIV  = 1E-4;
const double   OMEGA = 1.2;
const double   ALPHA = 1.0 / 24.0;
const int      MAXIT = 1000;
const int      NSTEP = 20000;
      double   AD    = 0;
      double   R     = 0;
      int      ITER  = 0;
      int      STEP  = 0;

#endif
