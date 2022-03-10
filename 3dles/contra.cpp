#include "topo.h"
#include "flist.h"

void uu(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double  UU[N1 + 4][N2 + 4][N3 + 4][3],
    double  KX[N1 + 4][N2 + 4][N3 + 4][3],
    double   J[N1 + 4][N2 + 4][N3 + 4]
) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                double U0U, U0V, U0W;
                double KX0X, KX0Y, KX0Z;
                double J0J;

                U0U  =  U[i][j][k][u];
                U0V  =  U[i][j][k][v];
                U0W  =  U[i][j][k][w];
                KX0X = KX[i][j][k][x];
                KX0Y = KX[i][j][k][y];
                KX0Z = KX[i][j][k][z];
                J0J  =  J[i][j][k];

                UU[i][j][k][u] = KX0X * J0J * U0U;
                UU[i][j][k][v] = KX0Y * J0J * U0V;
                UU[i][j][k][w] = KX0Z * J0J * U0W;
            }
        }
    }
}

void ju(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double  UD[N1 + 4][N2 + 4][N3 + 4][3],
    double  UU[N1 + 4][N2 + 4][N3 + 4][3],
    double  JU[N1 + 4][N2 + 4][N3 + 4][3],
    double BBU[N1 + 4][N2 + 4][N3 + 4][3][3],
    int    FBU[N1 + 4][N2 + 4][N3 + 4][3],
    double  KX[N1 + 4][N2 + 4][N3 + 4][3],
    double   J[N1 + 4][N2 + 4][N3 + 4],
    bool     Ph
) {
    for (int i = I1 - 1; i <= I2 + 1; i ++) {
        for (int j = J1 - 1; j <= J2 + 1; j ++) {
            for (int k = K1 - 1; k <= K2 + 1; k ++) {
                double U0, U1, UU0, UU1;
                double UF;
                int    FF;
                double KX1, J1J;

//  Ξ[ξ] direction

                UU0 =  UU[i    ][j][k][u];
                UU1 =  UU[i + 1][j][k][u];
                FF  = FBU[i    ][j][k][x];
                if (FF & BOUNDARY) {
                    U0  = (Ph)? (U[i][j][k][u]) : (UD[i][j][k][u]);
                    UF  = BBU[i][j][k][x][u];
                    U1  = 2 * UF - U0;
                    KX1 = KX[i + 1][j][k][x];
                    J1J =  J[i + 1][j][k];
                    UU1 = KX1 * J1J * U1;
                }
                JU[i][j][k][u] = 0.5 * (UU0 + UU1);

//  Ξ[η] direction

                UU0 =  UU[i][j    ][k][v];
                UU1 =  UU[i][j + 1][k][v];
                FF  = FBU[i][j    ][k][y];
                if (FF & BOUNDARY) {
                    U0  = (Ph)? (U[i][j][k][v]) : (UD[i][j][k][v]);
                    UF  = BBU[i][j][k][y][v];
                    U1  = 2 * UF - U0;
                    KX1 = KX[i][j + 1][k][y];
                    J1J =  J[i][j + 1][k];
                    UU1 = KX1 * J1J * U1;
                }
                JU[i][j][k][v] = 0.5 * (UU0 + UU1);

//  Ξ[ζ] direction

                UU0 =  UU[i][j][k    ][w];
                UU1 =  UU[i][j][k + 1][w];
                FF  = FBU[i][j][k    ][z];
                if (FF & BOUNDARY) {
                    U0  = (Ph)? (U[i][j][k][w]) : (UD[i][j][k][w]);
                    UF  = BBU[i][j][k][z][w];
                    U1  = 2 * UF - U0;
                    KX1 = KX[i][j][k + 1][z];
                    J1J =  J[i][j][k + 1];
                    UU1 = KX1 * J1J * U1;
                }
                JU[i][j][k][w] = 0.5 * (UU0 + UU1);
            }
        }
    }
}