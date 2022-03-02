#include "topo.h"

void contra(
    double  U[NX + 2][NY + 2][NZ + 2][3],
    double UC[NX + 2][NY + 2][NZ + 2][3],
    double UU[NX + 2][NY + 2][NZ + 2][3],
    double KX[NX + 2][NY + 2][NZ + 2][3],
    double  J[NX + 2][NY + 2][NZ + 2]
) {
    int    i, j, k;
    double JC0, K1X1, K2X2, K3X3;
    double UC0U, UC0V, UC0W;
    double UUC0, UUE1, UUN1, UUT1;
    double UUE, UUN, UUT;

//  contravariant velocity at cell center

    for (i = 1; i <= NX; i ++) {
        for (j = 1; j <= NY; j ++) {
            for (k = 1; k <= NZ; k ++) {
                UC0U =  U[i][j][k][0];
                UC0V =  U[i][j][k][1];
                UC0W =  U[i][j][k][2];
                JC0  =  J[i][j][k];
                K1X1 = KX[i][j][k][1];
                K2X2 = KX[i][j][k][2];
                K3X3 = KX[i][j][k][3];

                UC[i][j][k][0] = UC0U * JC0 * K1X1;
                UC[i][j][k][1] = UC0V * JC0 * K2X2;
                UC[i][j][k][2] = UC0W * JC0 * K3X3;
            }
        }
    }

//  contravariant velocity at cell faces

    for (i = 1; i <= NX; i ++) {
        for (j = 1; j <= NY; j ++) {
            for (k = 1; k <= NZ; k ++) {
                UC0U = UC[i    ][j    ][k    ][0];
                UC0V = UC[i    ][j    ][k    ][1];
                UC0W = UC[i    ][j    ][k    ][2];
                UUE1 = UC[i + 1][j    ][k    ][0];
                UUN1 = UC[i    ][j + 1][k    ][1];
                UUT1 = UC[i    ][j    ][k + 1][2];
                UUE  = 0.5 * (UC0U + UUE1);
                UUN  = 0.5 * (UC0V + UUN1);
                UUT  = 0.5 * (UC0W + UUT1);
                
                UU[i][j][k][0] = UUE;
                UU[i][j][k][1] = UUN;
                UU[i][j][k][2] = UUT;
            }
        }
    }
}