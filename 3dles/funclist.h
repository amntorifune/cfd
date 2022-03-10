#ifndef FUNCLIST_H
#define FUNCLIST_H

#include "topo.h"

extern void init(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double  UD[N1 + 4][N2 + 4][N3 + 4][3],
    double  UU[N1 + 4][N2 + 4][N3 + 4][3],
    double  JU[N1 + 4][N2 + 4][N3 + 4][3],
    double BBU[N1 + 4][N2 + 4][N3 + 4][3][3],
    double   P[N1 + 4][N2 + 4][N3 + 4],
    double BBP[N1 + 4][N2 + 4][N3 + 4][3],
    double DIV[N1 + 4][N2 + 4][N3 + 4],
    double   X[N1 + 4][N2 + 4][N3 + 4][3],
    double  KX[N1 + 4][N2 + 4][N3 + 4][3],
    double   J[N1 + 4][N2 + 4][N3 + 4],
    double   G[N1 + 4][N2 + 4][N3 + 4][3],
    double   C[N1 + 4][N2 + 4][N3 + 4][6],
    double SGS[N1 + 4][N2 + 4][N3 + 4],
    int      F[N1 + 4][N2 + 4][N3 + 4],
    int    FBU[N1 + 4][N2 + 4][N3 + 4][3],
    int    FBP[N1 + 4][N2 + 4][N3 + 4][3]
);

#endif