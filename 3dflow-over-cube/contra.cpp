#include "topo.h"
#include "flag.h"

void contra(
    double   U[NX + 4][NY + 4][NZ + 4][3],
    double  UC[NX + 4][NY + 4][NZ + 4][3],
    double  UU[NX + 4][NY + 4][NZ + 4][3],
    int    FFU[NX + 4][NY + 4][NZ + 4][3],
    double BBU[NX + 4][NY + 4][NZ + 4][3][3],
    double   X[NX + 4][NY + 4][NZ + 4][3],
    double  KX[NX + 4][NY + 4][NZ + 4][3],
    double   J[NX + 4][NY + 4][NZ + 4]
) {
    for (int i = I1 - 1; i <= I2 + 1; i ++) {
        for (int j = J1 - 1; j <= J2 + 1; j ++) {
            for (int k = K1 - 1; k <= K2 + 1; k ++) {
                double J0J;
                double K1X1, K2X2, K3X3;
                double U0U, U0V, U0W;

                U0U  =  U[i][j][k][0];
                U0V  =  U[i][j][k][1];
                U0W  =  U[i][j][k][2];
                K1X1 = KX[i][j][k][0];
                K2X2 = KX[i][j][k][1];
                K3X3 = KX[i][j][k][2];
                J0J  =  J[i][j][k];

                UC[i][j][k][0] = J0J * K1X1 * U0U;
                UC[i][j][k][1] = J0J * K2X2 * U0V;
                UC[i][j][k][2] = J0J * K3X3 * U0W;
            }
        }
    }

    for (int i = I1 - 1; i <= I2; i ++) {
        for (int j = J1 - 1; j <= J2; j ++) {
            for (int k = K1 - 1; k <= K2; k ++) {
                double UU0, UU1, UUF, UF;
                double UF, JF;
                double X1K1, X2K2, X3K3, K1X1, K2X2, K3X3;
                int    FF;
    
    //  east face

                FF = FFU[i][j][k][0];
                if (FF) {
                    X1K1 =         X[i + 1][j][k][0] -  X[i][j][k][0];
                    K2X2 = 0.5 * (KX[i + 1][j][k][1] + KX[i][j][k][1]);
                    K3X3 = 0.5 * (KX[i + 1][j][k][2] + KX[i][j][k][2]); 
                    K1X1 = 1 / X1K1;
                    X2K2 = 1 / K2X2;
                    X3K3 = 1 / K3X3;
                    JF   = X1K1 * X2K2 * X3K3;
                    UF   = BBU[i][j][k][0][0];
                    UUF  = JF * K1X1 * UF;
                }
                else {
                    UU0 = UC[i    ][j][k][0];
                    UU1 = UC[i + 1][j][k][0];
                    UUF = 0.5 * (UU0 + UU1);
                }
                UU[i][j][k][0] = UUF;
    
    //  north face

                FF = FFU[i][j][k][1];
                if (FF) {
                    K1X1 = 0.5 * (KX[i][j + 1][k][0] + KX[i][j][k][0]);
                    X2K2 =         X[i][j + 1][k][1] -  X[i][j][k][1];
                    K3X3 = 0.5 * (KX[i][j + 1][k][2] + KX[i][j][k][2]);
                    X1K1 = 1 / K1X1;
                    K2X2 = 1 / X2K2;
                    X3K3 = 1 / K3X3;
                    JF   = X1K1 * X2K2 * X3K3;
                    UF   = BBU[i][j][k][1][1];
                    UUF  = JF * K2X2 * UF;
                }
                else {
                    UU0 = UC[i][j    ][k][1];
                    UU1 = UC[i][j + 1][k][1];
                    UUF = 0.5 * (UU0 + UU1);
                }
                UU[i][j][k][1] = UUF;
    
    //  top face

                FF = FFU[i][j][k][2];
                if (FF) {
                    K1X1 = 0.5 * (KX[i][j][k + 1][0] + KX[i][j][k][0]);
                    K2X2 = 0.5 * (KX[i][j][k + 1][1] + KX[i][j][k][1]);
                    X3K3 =         X[i][j][k + 1][2] -  X[i][j][k][2];
                    X1K1 = 1 / K1X1;
                    X2K2 = 1 / K2X2;
                    K3X3 = 1 / X3K3;
                    JF   = X1K1 * X2K2 * X3K3;
                    UF   = BBU[i][j][k][2][2];
                    UUF  = JF * K3X3 * UF;
                }
                else {
                    UU0 = UC[i][j][k    ][2];
                    UU1 = UC[i][j][k + 1][2];
                    UUF = 0.5 * (UU0 + UU1);
                }
                UU[i][j][k][2] = UUF;
            }
        }
    }
}