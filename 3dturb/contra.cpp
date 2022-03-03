#include "topo.h"

void contra(
    double  U[NX + 2][NY + 2][NZ + 2][3],
    double UC[NX + 2][NY + 2][NZ + 2][3],
    double UU[NX + 2][NY + 2][NZ + 2][3],
    double KX[NX + 2][NY + 2][NZ + 2][3],
    double  J[NX + 2][NY + 2][NZ + 2]
) {
//  contravariant velocity at cell center

    #pragma acc kernels loop independent collapse(2) present(U, KX, J, UC)
    for (int i = 1; i <= NX; i ++) {
        for (int j = 1; j <= NY; j ++) {
            for (int k = 1; k <= NZ; k ++) {
                double JC0;
                double K1X1, K2X2, K3X3;
                double UC0U, UC0V, UC0W;

                UC0U =  U[i][j][k][0];
                UC0V =  U[i][j][k][1];
                UC0W =  U[i][j][k][2];
                K1X1 = KX[i][j][k][0];
                K2X2 = KX[i][j][k][1];
                K3X3 = KX[i][j][k][2];
                JC0  =  J[i][j][k];

                UC[i][j][k][0] = UC0U * JC0 * K1X1;
                UC[i][j][k][1] = UC0V * JC0 * K2X2;
                UC[i][j][k][2] = UC0W * JC0 * K3X3;
            }
        }
    }

//  contravariant velocity at cell faces

    #pragma acc kernels loop independent collapse(2) present(UC, UU)
    for (int i = 1; i <= NX; i ++) {
        for (int j = 1; j <= NY; j ++) {
            for (int k = 1; k <= NZ; k ++) {
                double UC0U, UC0V, UC0W;
                double UE1U, UN1V, UT1W;
                double UUE, UUN, UUT;

                UC0U = UC[i    ][j    ][k    ][0];
                UC0V = UC[i    ][j    ][k    ][1];
                UC0W = UC[i    ][j    ][k    ][2];
                UE1U = UC[i + 1][j    ][k    ][0];
                UN1V = UC[i    ][j + 1][k    ][1];
                UT1W = UC[i    ][j    ][k + 1][2];
                UUE  = 0.5 * (UC0U + UE1U);
                UUN  = 0.5 * (UC0V + UN1V);
                UUT  = 0.5 * (UC0W + UT1W);
                
                UU[i][j][k][0] = UUE;
                UU[i][j][k][1] = UUN;
                UU[i][j][k][2] = UUT;
            }
        }
    }
}