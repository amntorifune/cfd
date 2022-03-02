#include "topo.h"

void solid(
    double  U[NX + 2][NY + 2][NZ + 2][3],
    double UU[NX + 2][NY + 2][NZ + 2][3]
) {
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            for (int k = K1; k <= K2; k ++) {
                U[i][j][k][0] = 0.0;
                U[i][j][k][1] = 0.0;
                U[i][j][k][2] = 0.0;

                if (i < I2) {
                    UU[i][j][k][0] = 0.0;
                }
                if (j < J2) {
                    UU[i][j][k][1] = 0.0;
                }
                if (k < K2) {
                    UU[i][j][k][2] = 0.0;
                }
            }
        }
    }
}