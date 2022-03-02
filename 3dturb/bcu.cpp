#include "topo.h"

void bcu(
    double  U[NX + 2][NY + 2][NZ + 2][3],
    double UU[NX + 2][NY + 2][NZ + 2][3],
    double UD[NX + 2][NY + 2][NZ + 2][3],
    double KX[NX + 2][NY + 2][NZ + 2][3],
    double  DT
) {
    int    i, j, k;
    double K1X1;
    double DUX, DVX, DWX;
    double UC0, UW1, UW2;
    double VC0, VW1, VW2;
    double WC0, WW1, WW2;

//  right boundary : inlet, {u,v,w} = {1,0,0}
//  left boundary : outlet, special treatment

    for (j = 1; j <= NY; j ++) {
        for (k = 1; k <= NZ; k ++) {

            K1X1 = KX[NX    ][j][k][0];
            UC0  = UD[NX    ][j][k][0];
            VC0  = UD[NX    ][j][k][1];
            WC0  = UD[NX    ][j][k][2];
            UW1  = UD[NX - 1][j][k][0];
            VW1  = UD[NX - 1][j][k][1];
            WW1  = UD[NX - 1][j][k][2];
            UW2  = UD[NX - 2][j][k][0];
            VW2  = UD[NX - 2][j][k][1];
            WW2  = UD[NX - 2][j][k][2];
            DUX  = 0.5 * K1X1 * (3 * UC0 - 4 * UW1 + UW2);
            DVX  = 0.5 * K1X1 * (3 * VC0 - 4 * VW1 + VW2);
            DWX  = 0.5 * K1X1 * (3 * WC0 - 4 * WW1 + WW2);

            U[NX    ][j][k][0] =    UD[NX][j][k][0] - DT * DUX;
            U[NX    ][j][k][1] =    UD[NX][j][k][1] - DT * DVX;
            U[NX    ][j][k][2] =    UD[NX][j][k][2] - DT * DWX;
            U[NX + 1][j][k][0] = 2 * U[NX][j][k][0] - U[NX - 1][j][k][0];
            U[NX + 1][j][k][1] = 2 * U[NX][j][k][1] - U[NX - 1][j][k][1];
            U[NX + 1][j][k][2] = 2 * U[NX][j][k][2] - U[NX - 1][j][k][2];

            U[0     ][j][k][0] = 1.0;
            U[0     ][j][k][1] = 0.0;
            U[0     ][j][k][2] = 0.0;
            U[1     ][j][k][0] = 1.0;
            U[1     ][j][k][1] = 0.0;
            U[1     ][j][k][2] = 0.0;
        }
    }

//  front and back boundaries : free, reflective, slip wall

    for (i = 1; i <= NY; i ++) {
        for (k = 1; k <= NZ; k ++) {
            U[i][1     ][k][0] =   U[i][2     ][k][0];
            U[i][1     ][k][1] =   0.0;
            U[i][1     ][k][2] =   U[i][2     ][k][2];
            U[i][0     ][k][0] =   U[i][2     ][k][0];
            U[i][0     ][k][1] = - U[i][2     ][k][1];
            U[i][0     ][k][2] =   U[i][2     ][k][2];

            U[i][NY    ][k][0] =   U[i][NY - 1][k][0];
            U[i][NY    ][k][1] =   0.0;
            U[i][NY    ][k][2] =   U[i][NY - 1][k][2];
            U[i][NY + 1][k][0] =   U[i][NY - 1][k][0];
            U[i][NY + 1][k][1] = - U[i][NY - 1][k][1];
            U[i][NY + 1][k][2] =   U[i][NY - 1][k][2];
        }
    }

//  upper boundary : free, reflective, slip wall
//  lower boundary : non-slip wall, {u, v, w} = {0, 0, 0}

    for (i = 1; i <= NX; i ++) {
        for (j = 1; j <= NY; j ++) {
            U[i][j][1     ][0] =   0.0;
            U[i][j][1     ][1] =   0.0;
            U[i][j][1     ][2] =   0.0;
            U[i][j][0     ][0] = - U[i][j][2     ][0];
            U[i][j][0     ][1] = - U[i][j][2     ][1];
            U[i][j][0     ][2] = - U[i][j][2     ][2];

            U[i][j][NZ    ][0] =   U[i][j][NZ - 1][0];
            U[i][j][NZ    ][1] =   U[i][j][NZ - 1][1];
            U[i][j][NZ    ][2] =   0;
            U[i][j][NZ + 1][0] =   U[i][j][NZ - 1][0];
            U[i][j][NZ + 1][1] =   U[i][j][NZ - 1][1];
            U[i][j][NZ + 1][2] = - U[i][j][NZ - 1][2];
        }
    }

}