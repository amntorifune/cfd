#include "topo.h"

void outflow(
    double   U[NX + 4][NY + 4][NZ + 4][3],
    double BBU[NX + 4][NY + 4][NZ + 4][3][3],
    double   X[NX + 4][NY + 4][NZ + 4][3],
    double   DT
) {
    for (int j = J1; j <= J2; j ++) {
        for (int k = K1; k <= K2; k ++) {
            double X1K1, K1X1;
            double DUX, DVX, DWX;
            double UF0, VF0, WF0;
            double UW1, VW1, WW1;

            X1K1 = X[I2 + 1][j][k][x] - X[I2][j][k][x];
            K1X1 = 1 / X1K1;
            UF0  = BBU[I2][j][k][x][u];
            VF0  = BBU[I2][j][k][x][v];
            WF0  = BBU[I2][j][k][x][w];
            UW1  =   U[I2][j][k][u];
            VW1  =   U[I2][j][k][v];
            WW1  =   U[I2][j][k][w];
            DUX  = K1X1 * 2 * (UF0 - UW1);
            DVX  = K1X1 * 2 * (VF0 - VW1);
            DWX  = K1X1 * 2 * (WF0 - WW1);
            UF0  = UF0 - DT * DUX;
            VF0  = VF0 - DT * DVX;
            WF0  = WF0 - DT * DWX;

            BBU[I2][j][k][x][u] = UF0;
            BBU[I2][j][k][x][v] = VF0;
            BBU[I2][j][k][x][w] = WF0;
        }
    }
}

void slip(
    double   U[NX + 4][NY + 4][NZ + 4][3],
    double BBU[NX + 4][NY + 4][NZ + 4][3][3]
) {
//  top boundary
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            BBU[i][j][K2][z][u] = U[i][j][K2][u];
            BBU[i][j][K2][z][v] = U[i][j][K2][v];
            BBU[i][j][K2][z][w] = 0.0;
        }
    }

//  north boundary
    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            BBU[i][J2][k][y][u] = U[i][J2][k][u];
            BBU[i][J2][k][y][v] = 0.0;
            BBU[i][J2][k][y][w] = U[i][J2][k][w];
        }
    }

//  south boundary
    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            BBU[i][J1 - 1][k][y][u] = U[i][J1][k][u];
            BBU[i][J1 - 1][k][y][v] = 0.0;
            BBU[i][J1 - 1][k][y][w] = U[i][J1][k][w];
        }
    }
}