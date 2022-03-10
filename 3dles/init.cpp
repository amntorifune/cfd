#include "topo.h"
#include "flist.h"

void setupx(
    double X[N1 + 4][N2 + 4][N3 + 4][3]
) {
    double DX, DY, DZ;

    DX = (X2 - X1) / N1;
    DY = (Y2 - Y1) / N2;
    DZ = (Z2 - Z1) / N3;

    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                double XX, YY, ZZ;

                XX = (i - I1) * DX + 0.5 * DX;
                YY = (j - J1) * DY + 0.5 * DY;
                ZZ = (k - K1) * DZ + 0.5 * DZ;

                X[i][j][k][x] = XX;
                X[i][j][k][y] = YY;
                X[i][j][k][z] = ZZ;
            }
        }
    }
}

void setupkx(
    double   X[N1 + 4][N2 + 4][N3 + 4][3],
    double  KX[N1 + 4][N2 + 4][N3 + 4][3],
    double   J[N1 + 4][N2 + 4][N3 + 4],
    double   G[N1 + 4][N2 + 4][N3 + 4][3],
    double   C[N1 + 4][N2 + 4][N3 + 4][6]
) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                KX[i][j][k][x] = 0.0;
                KX[i][j][k][y] = 0.0;
                KX[i][j][k][z] = 0.0;
                J[i][j][k]     = 0.0;
                G[i][j][k][x]  = 0.0;
                G[i][j][k][y]  = 0.0;
                G[i][j][k][z]  = 0.0;
                C[i][j][k][0]  = 0.0;
                C[i][j][k][1]  = 0.0;
                C[i][j][k][2]  = 0.0;
                C[i][j][k][3]  = 0.0;
                C[i][j][k][4]  = 0.0;
                C[i][j][k][5]  = 0.0;
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

                XC0   = X[i    ][j    ][k    ][x];
                YC0   = X[i    ][j    ][k    ][y];
                ZC0   = X[i    ][j    ][k    ][z];
                XE1   = X[i + 1][j    ][k    ][x];
                XW1   = X[i - 1][j    ][k    ][x];
                YN1   = X[i    ][j + 1][k    ][y];
                YS1   = X[i    ][j - 1][k    ][y];
                ZT1   = X[i    ][j    ][k + 1][z];
                ZB1   = X[i    ][j    ][k - 1][z];

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

                KX[i][j][k][x] = K1X1;
                KX[i][j][k][y] = K2X2;
                KX[i][j][k][z] = K3X3;
                C[i][j][k][0]  = C1;
                C[i][j][k][1]  = C2;
                C[i][j][k][2]  = C3;
                C[i][j][k][3]  = C7;
                C[i][j][k][4]  = C8;
                C[i][j][k][5]  = C9;
                G[i][j][k][x]  = G11;
                G[i][j][k][y]  = G22;
                G[i][j][k][z]  = G33;
                J[i][j][k]     = DET;
            }
        }
    }
}

void setupf(
    int F[N1 + 4][N2 + 4][N3 + 4]
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
    int    FBU[N1 + 4][N2 + 4][N3 + 4][3],
    int    FBP[N1 + 4][N2 + 4][N3 + 4][3],
    double BBU[N1 + 4][N2 + 4][N3 + 4][3][3],
    double BBP[N1 + 4][N2 + 4][N3 + 4][3]
) {

    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                FBU[i][j][k][x]    = NOTHING;
                FBU[i][j][k][y]    = NOTHING;
                FBU[i][j][k][z]    = NOTHING;
                BBU[i][j][k][x][u] = 0.0;
                BBU[i][j][k][x][v] = 0.0;
                BBU[i][j][k][x][w] = 0.0;
                BBU[i][j][k][y][u] = 0.0;
                BBU[i][j][k][y][v] = 0.0;
                BBU[i][j][k][y][w] = 0.0;
                BBU[i][j][k][z][u] = 0.0;
                BBU[i][j][k][z][v] = 0.0;
                BBU[i][j][k][z][w] = 0.0;

                FBP[i][j][k][x]    = NOTHING;
                FBP[i][j][k][y]    = NOTHING;
                FBP[i][j][k][z]    = NOTHING;
                BBP[i][j][k][x]    = 0.0;
                BBP[i][j][k][y]    = 0.0;
                BBP[i][j][k][z]    = 0.0;
            }
        }
    }

//  inflow boundary : (u,v,w) = (1,0,0), dp/dx = 0

    for (int j = J1; j <= J2; j ++) {
        for (int k = K1; k <= K2; k ++) {
            FBU[I1 - 1][j][k][x]    = BOUNDARY;
            BBU[I1 - 1][j][k][x][u] = 1.0;
            BBU[I1 - 1][j][k][x][v] = 0.0;
            BBU[I1 - 1][j][k][x][w] = 0.0;

            FBP[I1 - 1][j][k][x]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  outflow boundary : (u,v,w) = (outflow), dp/dx = 0

    for (int j = J1; j <= J2; j ++) {
        for (int k = K1; k <= K2; k ++) {
            FBU[I2][j][k][x] = BOUNDARY;

    //  apply in outflow correction to velocity

            FBP[I2][j][k][x] = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  top boundary : (u,v,w) = (slip wall), dp/dz = 0

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            FBU[i][j][K2][z]    = BOUNDARY;

    //  apply slip wall condition to velocity

            FBP[i][j][K2][z]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  bottom boundary : (u,v,w) = (0,0,0), dp/dz = [dp/dz]wall

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (i >= I3 && i <= I4 && j >= J3 && j <= J4) {
                continue;
            }
            FBU[i][j][K1 - 1][z]    = BOUNDARY;
            BBU[i][j][K1 - 1][z][u] = 0.0;
            BBU[i][j][K1 - 1][z][v] = 0.0;
            BBU[i][j][K1 - 1][z][w] = 0.0;

            FBP[i][j][K1 - 1][z]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  north boundary : (du/dy,v,dw/dy) = (0,0,0), dp/dy = 0

    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            FBU[i][J2][k][y]    = BOUNDARY;

    //  apply slip wall condition to velocity

            FBP[i][J2][k][y]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  south boundary : (du/dy,v,dw/dy) = (0,0,0), dp/dy = 0

    for (int i = I1; i <= I2; i ++) {
        for (int k = K1; k <= K2; k ++) {
            FBU[i][J1 - 1][k][y]    = BOUNDARY;

    //  apply slip wall condition to velocity

            FBP[i][J1 - 1][k][y]    = BOUNDARY;

    //  apply zero gradient condition to pressure

        }
    }

//  cube top face : (u,v,w) = (0,0,0), dp/dz = [dp/dz]wall

    for (int i = I3; i <= I4; i ++) {
        for (int j = J3; j <= J4; j ++) {
            FBU[i][j][K4][z]    = BOUNDARY;
            BBU[i][j][K4][z][u] = 0.0;
            BBU[i][j][K4][z][v] = 0.0;
            BBU[i][j][K4][z][w] = 0.0;

            FBP[i][j][K4][z]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube east face : (u,v,w) = (0,0,0), dp/dx = [dp/dx]wall

    for (int j = J3; j <= J4; j ++) {
        for (int k = K3; k <= K4; k ++) {
            FBU[I4][j][k][x]    = BOUNDARY;
            BBU[I4][j][k][x][u] = 0.0;
            BBU[I4][j][k][x][v] = 0.0;
            BBU[I4][j][k][x][w] = 0.0;

            FBP[I4][j][k][x]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube west face : (u,v,w) = (0,0,0), dp/dx = [dp/dx]wall

    for (int j = J3; j <= J4; j ++) {
        for (int k = K3; k <= K4; k ++) {
            FBU[I3 - 1][j][k][x]    = BOUNDARY;
            BBU[I3 - 1][j][k][x][u] = 0.0;
            BBU[I3 - 1][j][k][x][v] = 0.0;
            BBU[I3 - 1][j][k][x][w] = 0.0;

            FBP[I3 - 1][j][k][x]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube north face : (u,v,w) = (0,0,0), dp/dy = [dp/dy]wall

    for (int i = I3; i <= I4; i ++) {
        for (int k = K3; k <= K4; k ++) {
            FBU[i][J4][k][y]    = BOUNDARY;
            BBU[i][J4][k][y][u] = 0.0;
            BBU[i][J4][k][y][v] = 0.0;
            BBU[i][J4][k][y][w] = 0.0;

            FBP[i][J4][k][y]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

//  cube south face : (u,v,w) = (0,0,0), dp/dy = [dp/dy]wall

    for (int i = I3; i <= I4; i ++) {
        for (int k = K3; k <= K4; k ++) {
            FBU[i][J3 - 1][k][y]    = BOUNDARY;
            BBU[i][J3 - 1][k][y][u] = 0.0;
            BBU[i][J3 - 1][k][y][v] = 0.0;
            BBU[i][J3 - 1][k][y][w] = 0.0;

            FBP[i][J3 - 1][k][y]    = BOUNDARY;

    //  apply wall NS equation to pressure

        }
    }

}

void clearvar(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double  UD[N1 + 4][N2 + 4][N3 + 4][3],
    double  UU[N1 + 4][N2 + 4][N3 + 4][3],
    double  JU[N1 + 4][N2 + 4][N3 + 4][3],
    double   P[N1 + 4][N2 + 4][N3 + 4],
    double DIV[N1 + 4][N2 + 4][N3 + 4],
    double SGS[N1 + 4][N2 + 4][N3 + 4]
) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                U[i][j][k][u]  = 0.0;
                U[i][j][k][v]  = 0.0;
                U[i][j][k][w]  = 0.0;
                UD[i][j][k][u] = 0.0;
                UD[i][j][k][v] = 0.0;
                UD[i][j][k][w] = 0.0;
                UU[i][j][k][u] = 0.0;
                UU[i][j][k][v] = 0.0;
                UU[i][j][k][w] = 0.0;
                JU[i][j][k][u] = 0.0;
                JU[i][j][k][v] = 0.0;
                JU[i][j][k][w] = 0.0;
                P[i][j][k]     = 0.0;
                DIV[i][j][k]   = 0.0;
                SGS[i][j][k]   = 0.0;
            }
        }
    }
}

void init(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double  UD[N1 + 4][N2 + 4][N3 + 4][3],
    double  UU[N1 + 4][N2 + 4][N3 + 4][3],
    double  JU[N1 + 4][N2 + 4][N3 + 4][3],
    double BBU[N1 + 4][N2 + 4][N3 + 4][3][3],
    double   P[N1 + 4][N2 + 4][N3 + 4],
    double BBP[N1 + 4][N2 + 4][N3 + 4][3],
    double DIV[N1 + 4][N2 + 4][N3 + 4],
    double   X[N1 + 4][N2 + 4][N3 + 4][3],
    double  KX[N1 + 4][N2 + 4][N3 + 4][3],
    double   J[N1 + 4][N2 + 4][N3 + 4],
    double   G[N1 + 4][N2 + 4][N3 + 4][3],
    double   C[N1 + 4][N2 + 4][N3 + 4][6],
    double SGS[N1 + 4][N2 + 4][N3 + 4],
    int      F[N1 + 4][N2 + 4][N3 + 4],
    int    FBU[N1 + 4][N2 + 4][N3 + 4][3],
    int    FBP[N1 + 4][N2 + 4][N3 + 4][3]
) {
    clearvar(U, UD, UU, JU, P, DIV, SGS);
    setupx(X);
    setupkx(X, KX, J, G, C);
    setupf(F);
    setupbb(FBU, FBP, BBU, BBP);
}
