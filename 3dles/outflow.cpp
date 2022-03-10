#include "topo.h"

void outflow(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double BBU[N1 + 4][N2 + 4][N3 + 4][3][3],
    double   X[N1 + 4][N2 + 4][N3 + 4][3],
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