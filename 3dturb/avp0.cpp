#include "topo.h"

void avp0(
    double P[NX + 2][NY + 2][NZ + 2]
) {
    int    i, j, k;
    double AVG = 0;

    for (i = 1; i <= NX; i ++) {
        for (j = 1; j <= NY; j ++) {
            for (k = 1; k <= NZ; k ++) {
                AVG += P[i][j][k];
            }
        }
    }
    AVG = AVG / NXYZ;
    for (i = 1; i <= NX; i ++) {
        for (j = 1; j <= NY; j ++) {
            for (k = 1; k <= NZ; k ++) {
                P[i][j][k] -= AVG;
            }
        }
    }
}