#include "topo.h"

void slip(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double BBU[N1 + 4][N2 + 4][N3 + 4][3][3]
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