#include "topo.h"

void bcp(
    double  P[NX + 2][NY + 2][NZ + 2],
    double  U[NX + 2][NY + 2][NZ + 2][3],
    double  C[NX + 2][NY + 2][NZ + 2][6],
    double KX[NX + 2][NY + 2][NZ + 2][3],
    double  RI
) {
    double LAP, DK1, DK2, DK3, D2K1, D2K2, D2K3;
    double WC0, WE1, WW1, WN1, WS1, WT1, WB1;
    double C1, C2, C3, C7, C8, C9;
    double X3K3;

//  right and left boundaries : inlet and outlet, zero gradient

    for (int j = 1; j <= NY ; j ++) {
        for (int k = 1; k <= NZ; k ++) {
            P[1 ][j][k] = P[2     ][j][k];
            P[NX][j][k] = P[NX - 1][j][k];
        }
    }

//  front and back boundaries : free, zero gradient

    for (int i = 1; i <= NX; i ++) {
        for (int k = 1; k <= NZ; k ++) {
            P[i][1 ][k] = P[i][2     ][k];
            P[i][NY][k] = P[i][NY - 1][k];
        }
    }

//  upper boundary : free, zero gradient
//  lower boundary : wall, use a wall NS equation

    for (int i = 1; i <= NX; i ++) {
        for (int j = 1; j <= NY; j ++) {
            C1   = C[i    ][j    ][1    ][0];
            C2   = C[i    ][j    ][1    ][1];
            C3   = C[i    ][j    ][1    ][2];
            C7   = C[i    ][j    ][1    ][3];
            C8   = C[i    ][j    ][1    ][4];
            C9   = C[i    ][j    ][1    ][5];
            WC0  = U[i    ][j    ][1    ][2];
            WE1  = U[i + 1][j    ][1    ][2];
            WW1  = U[i - 1][j    ][1    ][2];
            WN1  = U[i    ][j + 1][1    ][2];
            WS1  = U[i    ][j - 1][1    ][2];
            WT1  = U[i    ][j    ][1 + 1][2];
            WB1  = U[i    ][j    ][1 - 1][2];
            X3K3 = 1 / KX[i][j][1][2];
            DK1  = 0.5 * (WE1 - WW1);
            DK2  = 0.5 * (WN1 - WS1);
            DK3  = 0.5 * (WT1 - WB1);
            D2K1 = WE1 - 2 * WC0 + WW1;
            D2K2 = WN1 - 2 * WC0 + WS1;
            D2K3 = WT1 - 2 * WC0 + WB1;
            LAP  = C1 * D2K1 + C2 * D2K2 + C3 * D2K3 - C7 * DK1 - C8 * DK2 - C9 * DK3;

            P[i][j][1 ] = P[i][j][1  + 1] - X3K3 * RI * LAP;
            P[i][j][NZ] = P[i][j][NZ - 1];
        }
    }
}