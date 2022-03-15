#include "topo.h"
#include "flag.h"

void setupf(
    int F[NX + 4][NY + 4][NZ + 4]
) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                F[i][j][k] = NONFLUID;
            }
        }
    }

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            for (int k = K1; k <= K2; k ++) {
                F[i][j][k] = FLUID;
            }
        }
    }

    for (int i = I3; i <= I4; i ++) {
        for (int j = J3; j <= J4; j ++) {
            for (int k = K3; k <= K4; k ++) {
                F[i][j][k] = NONFLUID;
            }
        }
    }
}

void setupbb(
    int    FFU[NX + 4][NY + 4][NZ + 4][3],
    int    FFP[NX + 4][NY + 4][NZ + 4][3],
    double BBU[NX + 4][NY + 4][NZ + 4][3][3],
    double BBP[NX + 4][NY + 4][NZ + 4][3]
) {

    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                FFU[i][j][k][0]    = NOTHING;
                FFU[i][j][k][1]    = NOTHING;
                FFU[i][j][k][2]    = NOTHING;
                BBU[i][j][k][0][0] = 0.0;
                BBU[i][j][k][0][1] = 0.0;
                BBU[i][j][k][0][2] = 0.0;
                BBU[i][j][k][1][0] = 0.0;
                BBU[i][j][k][1][1] = 0.0;
                BBU[i][j][k][1][2] = 0.0;
                BBU[i][j][k][2][0] = 0.0;
                BBU[i][j][k][2][1] = 0.0;
                BBU[i][j][k][2][2] = 0.0;

                FFP[i][j][k][0]    = NOTHING;
                FFP[i][j][k][1]    = NOTHING;
                FFP[i][j][k][2]    = NOTHING;
                BBP[i][j][k][0]    = 0.0;
                BBP[i][j][k][1]    = 0.0;
                BBP[i][j][k][2]    = 0.0;
            }
        }
    }

//  inflow boundary : (u,v,w) = (1,0,0), dp/dx = 0

    for (int j = J1; j <= J2; j ++) {
        for (int k = K1; k <= K2; k ++) {
            FFU[I1 - 1][j][k][0]    = BOUNDARY;
            BBU[I1 - 1][j][k][0][0] = 1.0;
            BBU[I1 - 1][j][k][0][1] = 0.0;
            BBU[I1 - 1][j][k][0][2] = 0.0;

            FFP[I1 - 1][j][k][0]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  outflow boundary : (u,v,w) = (outflow), dp/dx = 0

    for (int j = J1; j <= J2; j ++) {
        for (int k = K1; k <= K2; k ++) {
            FFU[I2][j][k][0] = BOUNDARY;

    //  apply in outflow correction to velocity

            FFP[I2][j][k][0] = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  top boundary : (u,v,w) = (slip wall), dp/dz = 0

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            FFU[i][j][K2][2]    = BOUNDARY;

    //  apply slip wall condition to velocity

            FFP[i][j][K2][2]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  bottom boundary : (u,v,w) = (0,0,0), dp/dz = [dp/dz]wall

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (i >= I3 && i <= I4 && j >= J3 && j <= J4) {
                continue;
            }
            FFU[i][j][K1 - 1][2]    = BOUNDARY;
            BBU[i][j][K1 - 1][2][0] = 0.0;
            BBU[i][j][K1 - 1][2][1] = 0.0;
            BBU[i][j][K1 - 1][2][2] = 0.0;

            FFP[i][j][K1 - 1][2]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  north boundary : (du/dy,v,dw/dy) = (0,0,0), dp/dy = 0

    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            FFU[i][J2][k][1]    = BOUNDARY;

    //  apply slip wall condition to velocity

            FFP[i][J2][k][1]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  south boundary : (du/dy,v,dw/dy) = (0,0,0), dp/dy = 0

    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            FFU[i][J1 - 1][k][1]    = BOUNDARY;

    //  apply slip wall condition to velocity

            FFP[i][J1 - 1][k][1]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  cube top face : (u,v,w) = (0,0,0), dp/dz = [dp/dz]wall

    for (int i = I3; i <= I4; i ++) {
        for (int j = J3; j <= J4; j ++) {
            FFU[i][j][K4][2]    = BOUNDARY;
            BBU[i][j][K4][2][0] = 0.0;
            BBU[i][j][K4][2][1] = 0.0;
            BBU[i][j][K4][2][2] = 0.0;

            FFP[i][j][K4][2]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube east face : (u,v,w) = (0,0,0), dp/dx = [dp/dx]wall

    for (int j = J3; j <= J4; j ++) {
        for (int k = K3; k <= K4; k ++) {
            FFU[I4][j][k][0]    = BOUNDARY;
            BBU[I4][j][k][0][0] = 0.0;
            BBU[I4][j][k][0][1] = 0.0;
            BBU[I4][j][k][0][2] = 0.0;

            FFP[I4][j][k][0]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube west face : (u,v,w) = (0,0,0), dp/dx = [dp/dx]wall

    for (int j = J3; j <= J4; j ++) {
        for (int k = K3; k <= K4; k ++) {
            FFU[I3 - 1][j][k][0]    = BOUNDARY;
            BBU[I3 - 1][j][k][0][0] = 0.0;
            BBU[I3 - 1][j][k][0][1] = 0.0;
            BBU[I3 - 1][j][k][0][2] = 0.0;

            FFP[I3 - 1][j][k][0]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube north face : (u,v,w) = (0,0,0), dp/dy = [dp/dy]wall

    for (int i = I3; i <= I4; i ++) {
        for (int k = K3; k <= K4; k ++) {
            FFU[i][J4][k][1]    = BOUNDARY;
            BBU[i][J4][k][1][0] = 0.0;
            BBU[i][J4][k][1][1] = 0.0;
            BBU[i][J4][k][1][2] = 0.0;

            FFP[i][J4][k][1]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube south face : (u,v,w) = (0,0,0), dp/dy = [dp/dy]wall

    for (int i = I3; i <= I4; i ++) {
        for (int k = K3; k <= K4; k ++) {
            FFU[i][J3 - 1][k][1]    = BOUNDARY;
            BBU[i][J3 - 1][k][1][0] = 0.0;
            BBU[i][J3 - 1][k][1][1] = 0.0;
            BBU[i][J3 - 1][k][1][2] = 0.0;

            FFP[i][J3 - 1][k][1]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

}

void setupx(
    double   X[NX + 4][NY + 4][NZ + 4][3],
    double  KX[NX + 4][NY + 4][NZ + 4][3],
    double   J[NX + 4][NY + 4][NZ + 4],
    double   G[NX + 4][NY + 4][NZ + 4][3],
    double   C[NX + 4][NY + 4][NZ + 4][6]
) {
    double DX, DY, DZ;

    DX = (X2 - X1) / NX;
    DY = (Y2 - Y1) / NY;
    DZ = (Z2 - Z1) / NZ;

    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                double XX, YY, ZZ;

                XX = (i - I1) * DX + 0.5 * DX;
                YY = (j - J1) * DY + 0.5 * DY;
                ZZ = (k - K1) * DZ + 0.5 * DZ;

                X[i][j][k][0] = XX;
                X[i][j][k][1] = YY;
                X[i][j][k][2] = ZZ;
            }
        }
    }

    for (int i = I1 - 1; i <= I2 + 1; i ++) {
        for (int j = J1 - 1; j <= J2 + 1; j ++) {
            for (int k = K1 - 1; k <= K2 + 1; k ++) {
                double XC0, XE1, XW1;
                double YC0, YN1, YS1;
                double ZC0, ZT1, ZB1;
                double X1K1, X2K2, X3K3;
                double X1K11, X2K22, X3K33;
                double K1X1, K2X2, K3X3;
                double C1, C2, C3, C7, C8, C9;
                double G11, G22, G33;
                double DET;

                XC0   = X[i    ][j    ][k    ][0];
                YC0   = X[i    ][j    ][k    ][1];
                ZC0   = X[i    ][j    ][k    ][2];
                XE1   = X[i + 1][j    ][k    ][0];
                XW1   = X[i - 1][j    ][k    ][0];
                YN1   = X[i    ][j + 1][k    ][1];
                YS1   = X[i    ][j - 1][k    ][1];
                ZT1   = X[i    ][j    ][k + 1][2];
                ZB1   = X[i    ][j    ][k - 1][2];

                X1K1  = 0.5 * (XE1 - XW1);
                X2K2  = 0.5 * (YN1 - YS1);
                X3K3  = 0.5 * (ZT1 - ZB1);
                X1K11 = XE1 - 2 * XC0 + XW1;
                X2K22 = YN1 - 2 * YC0 + YS1;
                X3K33 = ZT1 - 2 * ZC0 + ZB1;
                K1X1  = 1 / X1K1;
                K2X2  = 1 / X2K2;
                K3X3  = 1 / X3K3;
                DET   = X1K1 * X2K2 * X3K3;

                C1    = 1 / (X1K1 * X1K1);
                C2    = 1 / (X2K2 * X2K2);
                C3    = 1 / (X3K3 * X3K3);
                C7    = X1K11 / (X1K1 * X1K1 * X1K1);
                C8    = X2K22 / (X2K2 * X2K2 * X2K2);
                C9    = X3K33 / (X3K3 * X3K3 * X3K3);
                
                G11   = DET * K1X1 * K1X1;
                G22   = DET * K2X2 * K2X2;
                G33   = DET * K3X3 * K3X3;

                KX[i][j][k][0] = K1X1;
                KX[i][j][k][1] = K2X2;
                KX[i][j][k][2] = K3X3;
                C[i][j][k][0]  = C1;
                C[i][j][k][1]  = C2;
                C[i][j][k][2]  = C3;
                C[i][j][k][3]  = C7;
                C[i][j][k][4]  = C8;
                C[i][j][k][5]  = C9;
                G[i][j][k][0]  = G11;
                G[i][j][k][1]  = G22;
                G[i][j][k][2]  = G33;
                J[i][j][k]     = DET;
            }
        }
    }
}


void init(
    double   U[NX + 4][NY + 4][NZ + 4][3],
    double  UA[NX + 4][NY + 4][NZ + 4][3],
    double  UC[NX + 4][NY + 4][NZ + 4][3],
    double  UU[NX + 4][NY + 4][NZ + 4][3],
    double UUA[NX + 4][NY + 4][NZ + 4][3],
    int    FFU[NX + 4][NY + 4][NZ + 4][3],
    double BBU[NX + 4][NY + 4][NZ + 4][3][3],
    double   P[NX + 4][NY + 4][NZ + 4],
    int    FFP[NX + 4][NY + 4][NZ + 4][3],
    double BBP[NX + 4][NY + 4][NZ + 4][3],
    double DIV[NX + 4][NY + 4][NZ + 4],
    double   X[NX + 4][NY + 4][NZ + 4][3],
    double  KX[NX + 4][NY + 4][NZ + 4][3],
    double   J[NX + 4][NY + 4][NZ + 4],
    double   G[NX + 4][NY + 4][NZ + 4][3],
    double   C[NX + 4][NY + 4][NZ + 4][6],
    double SGS[NX + 4][NY + 4][NZ + 4],
    int      F[NX + 4][NY + 4][NZ + 4]
) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                U[i][j][k][0]      = 0;
                U[i][j][k][1]      = 0;
                U[i][j][k][2]      = 0;
                UA[i][j][k][0]     = 0;
                UA[i][j][k][1]     = 0;
                UA[i][j][k][2]     = 0;
                UC[i][j][k][0]     = 0;
                UC[i][j][k][1]     = 0;
                UC[i][j][k][2]     = 0;
                UU[i][j][k][0]     = 0;
                UU[i][j][k][1]     = 0;
                UU[i][j][k][2]     = 0;
                UUA[i][j][k][0]    = 0;
                UUA[i][j][k][1]    = 0;
                UUA[i][j][k][2]    = 0;
                FFU[i][j][k][0]    = 0;
                FFU[i][j][k][1]    = 0;
                FFU[i][j][k][2]    = 0;
                BBU[i][j][k][0][0] = 0;
                BBU[i][j][k][0][1] = 0;
                BBU[i][j][k][0][2] = 0;
                BBU[i][j][k][1][0] = 0;
                BBU[i][j][k][1][1] = 0;
                BBU[i][j][k][1][2] = 0;
                BBU[i][j][k][2][0] = 0;
                BBU[i][j][k][2][1] = 0;
                BBU[i][j][k][2][2] = 0;

                P[i][j][k]         = 0;
                FFP[i][j][k][0]    = 0;
                FFP[i][j][k][1]    = 0;
                FFP[i][j][k][2]    = 0;
                BBP[i][j][k][0]    = 0;
                BBP[i][j][k][1]    = 0;
                BBP[i][j][k][2]    = 0;
                DIV[i][j][k]       = 0;
                X[i][j][k][0]      = 0;
                X[i][j][k][1]      = 0;
                X[i][j][k][2]      = 0;
                KX[i][j][k][0]     = 0;
                KX[i][j][k][1]     = 0;
                KX[i][j][k][2]     = 0;
                J[i][j][k]         = 0;
                G[i][j][k][0]      = 0;
                G[i][j][k][1]      = 0;
                G[i][j][k][2]      = 0;
                C[i][j][k][0]      = 0;
                C[i][j][k][1]      = 0;
                C[i][j][k][2]      = 0;
                C[i][j][k][3]      = 0;
                C[i][j][k][4]      = 0;
                C[i][j][k][5]      = 0;
                SGS[i][j][k]       = 0;
                F[i][j][k]         = 0;
            }
        }
    }

    setupx(X, KX, J, G, C);
    setupf(F);
    setupbb(FFU, FFP, BBU, BBP);
}
