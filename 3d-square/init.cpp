#include "topo.h"

void init(
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
) {
    int    i, j, k;
    double XC0, XE1, XW1;
    double YC0, YN1, YS1;
    double ZC0, ZT1, ZB1;
    double X1K1, X2K2, X3K3, X1KK1, X2KK2, X3KK3, DET;
    double K1X1, K2X2, K3X3;
    double G11, G22, G33;
    double C1, C2, C3, C7, C8, C9;

    for (i = 1; i <= NX; i ++) {
        for (j = 1; j <= NY; j ++) {
            for (k = 1; k <= NZ ; k ++) {
                XC0   = X[i    ][j    ][k    ][0];
                YC0   = X[i    ][j    ][k    ][1];
                ZC0   = X[i    ][j    ][k    ][2];
                XE1   = X[i + 1][j    ][k    ][0];
                XW1   = X[i - 1][j    ][k    ][0];
                YN1   = X[i    ][j + 1][k    ][1];
                YS1   = X[i    ][j - 1][k    ][1];
                ZT1   = X[i    ][j    ][k + 1][2];
                ZB1   = X[i    ][j    ][k - 1][2];

                X1K1  = 0.5 * (XE1 - XW1);
                X2K2  = 0.5 * (YN1 - YS1);
                X3K3  = 0.5 * (ZT1 - ZB1);
                X1KK1 = XE1 - 2 * XC0 + XW1;
                X2KK2 = YN1 - 2 * YC0 + YS1;
                X3KK3 = ZT1 - 2 * ZC0 + ZB1;
                K1X1  = 1 / X1K1;
                K2X2  = 1 / X2K2;
                K3X3  = 1 / X3K3;
                DET   = X1K1 * X2K2 * X3K3;

                C1    = 1 / (X1K1 * X1K1);
                C2    = 1 / (X2K2 * X2K2);
                C3    = 1 / (X3K3 * X3K3);
                C7    = X1KK1 / (X1K1 * X1K1 * X1K1);
                C8    = X2KK2 / (X2K2 * X2K2 * X2K2);
                C9    = X3KK3 / (X3K3 * X3K3 * X3K3);
                
                G11   = DET * K1X1 * K1X1;
                G22   = DET * K2X2 * K2X2;
                G33   = DET * K3X3 * K3X3;

                KX[i][j][k][0] = K1X1;
                KX[i][j][k][1] = K2X2;
                KX[i][j][k][2] = K3X3;
                
                J[i][j][k] = DET;

                C[i][j][k][0] = C1;
                C[i][j][k][1] = C2;
                C[i][j][k][2] = C3;
                C[i][j][k][3] = C7;
                C[i][j][k][4] = C8;
                C[i][j][k][5] = C9;

                G[i][j][k][0] = G11;
                G[i][j][k][1] = G22;
                G[i][j][k][2] = G33;

                U[i][j][k][0]  = 0.0;
                U[i][j][k][1]  = 0.0;
                U[i][j][k][2]  = 0.0;
                UD[i][j][k][0] = 0.0;
                UD[i][j][k][1] = 0.0;
                UD[i][j][k][2] = 0.0;
                UU[i][j][k][0] = 0.0;
                UU[i][j][k][1] = 0.0;
                UU[i][j][k][2] = 0.0;
                P[i][j][k]     = 0.0;
                SGS[i][j][k]   = 0.0;
            }
        }
    }
}