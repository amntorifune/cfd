#include "topo.h"

void p0g(
    double   P[N1 + 4][N2 + 4][N3 + 4],
    double BBP[N1 + 4][N2 + 4][N3 + 4][3]
) {
//  top boundary

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            BBP[i][j][K2][z] = P[i][j][K2];
        }
    }

//  east boundary

    for (int j = J1; j <= J2; j ++) {
        for (int k = K1; k <= K2; k ++) {
            BBP[I2][j][k][x] = P[I2][j][k];
        }
    }

//  west boundary

    for (int j = J1; j <= J2; j ++) {
        for (int k = K1; k <= K2; k ++) {
            BBP[I1 - 1][j][k][x] = P[I1][j][k];
        }
    }

//  north boundary

    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            BBP[i][J2][k][y] = P[i][J2][k];
        }
    }

//  south face

    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            BBP[i][J1 - 1][k][y] = P[i][J1][k];
        }
    }
}