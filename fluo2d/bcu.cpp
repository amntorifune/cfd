#include "topo.h"
#include "flag.h"

void inflow(
    double BBU[NX + 4][NY + 4][2][2],
    double t
) {
    double inflowu = (t < 5.0)? 1.0 * (t / 5.0) : 1.0;
    #pragma acc kernels loop independent present(BBU) copyin(inflowu)
    for (int j = J1; j <= J2; j ++) {
        BBU[I1 - 1][j][0][0] = inflowu;
        BBU[I1 - 1][j][0][1] = 0.0;
    }
}

void outflow(
    double   U[NX + 4][NY + 4][2],
    double  UU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double   X[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4],
    double   DT
) {
    double UM = 1.0;
    #pragma acc kernels loop independent present(U, UU, BBU, X, J) copyin(DT, UM)
    for (int j = J1; j <= J2; j ++) {
        double XKF, KXF, UF, UUF, UOF, JF;
        JF  = 0.5 * ( J[I2][j]    +  J[I2 - 1][j]   );
        XKF =         X[I2][j][0] -  X[I2 - 1][j][0];
        KXF = 1 / XKF;
        // KXF = 0.5 * (KX[I2][j][0] + KX[I2 - 1][j][0]);
        UUF = UU[I2 - 1][j][0];
        UF  = UUF / (KXF * JF);
        UOF = BBU[I2][j][0][0];
        
        XKF = X[I2 + 1][j][0] - X[I2][j][0];
        KXF = 1 / XKF;
        UOF = UOF - DT * UM * KXF * (UOF - UF);
        BBU[I2][j][0][0] = UOF;

        UOF = BBU[I2][j][0][1];
        UF  =   U[I2][j][1];
        UOF = UOF - DT * UM * KXF * 2 * (UOF - UF);
        BBU[I2][j][0][1] = UOF;
    }
}

void slip(
    double   U[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2]
) {
//  north boundary

    #pragma acc kernels loop independent present(BBU, U)
    for (int i = I1; i <= I2; i ++) {
            BBU[i][J2][1][0] = U[i][J2][0];
            BBU[i][J2][1][1] = 0.0;
    }

//  south boundary

    #pragma acc kernels loop independent present(BBU, U)
    for (int i = I1; i <= I2; i ++) {
            BBU[i][J1 - 1][1][0] = U[i][J1][0];
            BBU[i][J1 - 1][1][1] = 0.0;
    }
}
