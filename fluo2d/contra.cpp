#include "topo.h"
#include "flag.h"

void contra(
    double   U[NX + 4][NY + 4][2],
    double  UC[NX + 4][NY + 4][2],
    double  UU[NX + 4][NY + 4][2],
    int    FFU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double   X[NX + 4][NY + 4][2],
    double  KX[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4]
) {
    #pragma acc kernels loop independent collapse(2) present(U, UC, KX, J)
    for (int i = I1 - 1; i <= I2 + 1; i ++) {
        for (int j = J1 - 1; j <= J2 + 1; j ++) {
            double J0J;
            double K1X1, K2X2;
            double U0U, U0V;

            U0U  =  U[i][j][0];
            U0V  =  U[i][j][1];
            K1X1 = KX[i][j][0];
            K2X2 = KX[i][j][1];
            J0J  =  J[i][j];

            UC[i][j][0] = J0J * K1X1 * U0U;
            UC[i][j][1] = J0J * K2X2 * U0V;
            
        }
    }

    #pragma acc kernels loop independent collapse(2) present(UC, UU, FFU, BBU, X, KX)
    for (int i = I1 - 1; i <= I2; i ++) {
        for (int j = J1 - 1; j <= J2; j ++) {
            double UU0, UU1, UUF, UF;
            double JF;
            double X1K1, X2K2, K1X1, K2X2;
            int    FF;
    
    //  east face

            FF = FFU[i][j][0];
            if (FF == BOUNDARY) {
                X1K1 =         X[i + 1][j][0] -  X[i][j][0];
                K2X2 = 0.5 * (KX[i + 1][j][1] + KX[i][j][1]);
                K1X1 = 1 / X1K1;
                X2K2 = 1 / K2X2;
                JF   = X1K1 * X2K2;
                UF   = BBU[i][j][0][0];
                UUF  = JF * K1X1 * UF;
            }
            else {
                UU0 = UC[i    ][j][0];
                UU1 = UC[i + 1][j][0];
                UUF = 0.5 * (UU0 + UU1);
            }
            UU[i][j][0] = UUF;
    
    //  north face

            FF = FFU[i][j][1];
            if (FF == BOUNDARY) {
                K1X1 = 0.5 * (KX[i][j + 1][0] + KX[i][j][0]);
                X2K2 =         X[i][j + 1][1] -  X[i][j][1];
                X1K1 = 1 / K1X1;
                K2X2 = 1 / X2K2;
                JF   = X1K1 * X2K2;
                UF   = BBU[i][j][1][1];
                UUF  = JF * K2X2 * UF;
            }
            else {
                UU0 = UC[i][j    ][1];
                UU1 = UC[i][j + 1][1];
                UUF = 0.5 * (UU0 + UU1);
            }
            UU[i][j][1] = UUF;
        }
    }
}
