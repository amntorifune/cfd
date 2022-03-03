#include "topo.h"

void avp0(
    double P[NX + 2][NY + 2][NZ + 2]
) {
    double AVG = 0;

    for (int i = 1; i <= NX; i ++) {
        for (int j = 1; j <= NY; j ++) {
            for (int k = 1; k <= NZ; k ++) {
                AVG += P[i][j][k];
            }
        }
    }

    AVG = AVG / NXYZ;
    
    for (int i = 1; i <= NX; i ++) {
        for (int j = 1; j <= NY; j ++) {
            for (int k = 1; k <= NZ; k ++) {
                P[i][j][k] -= AVG;
            }
        }
    }
}