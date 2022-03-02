#include "topo.h"

void geom(
    double X[NX + 2][NY + 2][NZ + 2][3]
) {
    double DX, DY, DZ;
    int    i, j, k;

    DX = LX / (NX - 1);
    DY = LY / (NY - 1);
    DZ = LZ / (NZ - 1);

    for (i = 0; i <= NX + 1; i ++) {
        for (j = 0; j <= NY + 1; j ++) {
            for (k = 0; k <= NZ + 1; k ++) {
                X[i][j][k][0] = DX * (i - 1);
                X[i][j][k][1] = DY * (j - 1);
                X[i][j][k][2] = DZ * (k - 1);
            }
        }
    }
}