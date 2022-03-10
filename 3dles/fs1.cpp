#include "topo.h"
#include "flist.h"
#include <math.h>

//  advection term

double adv(
    double UC0,
    double UE1,
    double UE2,
    double UW1,
    double UW2,
    double UN1,
    double UN2,
    double US1,
    double US2,
    double UT1,
    double UT2,
    double UB1,
    double UB2,
    double JUE,
    double JUW,
    double JUN,
    double JUS,
    double JUT,
    double JUB,
    double ABX,
    double ABY,
    double ABZ,
    double JC0,
    double ALPHA
) {
    double ADVX, ADVY, ADVZ;

    ADVX = JUE * (- UE2 + 27 * UE1 - 27 * UC0 + UW1);
    ADVX = JUW * (- UE1 + 27 * UC0 - 27 * UW1 + UW2) + ADVX;
    ADVX = ADVX / (48.0 * JC0);
    ADVX = ADVX + ALPHA * ABX * (UE2 - 4 * UE1 + 6 * UC0 - 4 * UW1 + UW2);

    ADVY = JUN * (- UN2 + 27 * UN1 - 27 * UC0 + US1);
    ADVY = JUS * (- UN1 + 27 * UC0 - 27 * US1 + US2) + ADVY;
    ADVY = ADVY / (48.0 * JC0);
    ADVY = ADVY + ALPHA * ABY * (UN2 - 4 * UN1 + 6 * UC0 - 4 * US1 + US2);

    ADVZ = JUT * (- UT2 + 27 * UT1 - 27 * UC0 + UB1);
    ADVZ = JUB * (- UT1 + 27 * UC0 - 27 * UB1 + UB2) + ADVZ;
    ADVZ = ADVZ / (48.0 * JC0);
    ADVZ = ADVZ + ALPHA * ABZ * (UT2 - 4 * UT1 + 6 * UC0 - 4 * UB1 + UB2);

    return (ADVX + ADVY + ADVZ);
}

// viscoucity term

double vis(
    double UC0,
    double UE1,
    double UW1,
    double UN1,
    double US1,
    double UT1,
    double UB1,
    double NUE,
    double RI,
    double C1,
    double C2,
    double C3,
    double C7,
    double C8,
    double C9
) {
    double LAP, DK1, DK2, DK3, D2K1, D2K2, D2K3;

    DK1  = 0.5 * (UE1 - UW1);
    DK2  = 0.5 * (UN1 - US1);
    DK3  = 0.5 * (UT1 - UB1);
    D2K1 = UE1 - 2 * UC0 + UW1;
    D2K2 = UN1 - 2 * UC0 + US1;
    D2K3 = UT1 - 2 * UC0 + UB1;
    LAP  = C1 * D2K1 + C2 * D2K2 + C3 * D2K3 - C7 * DK1 - C8 * DK2 - C9 * DK3;

    return (NUE + RI) * LAP;
}

//  pseudo velocity

void fs1(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double  UD[N1 + 4][N2 + 4][N3 + 4][3],
    double  UU[N1 + 4][N2 + 4][N3 + 4][3],
    double  JU[N1 + 4][N2 + 4][N3 + 4][3],
    int    FBU[N1 + 4][N2 + 4][N3 + 4][3],
    double BBU[N1 + 4][N2 + 4][N3 + 4][3][3],
    double  KX[N1 + 4][N2 + 4][N3 + 4][3],
    double   J[N1 + 4][N2 + 4][N3 + 4],
    double   C[N1 + 4][N2 + 4][N3 + 4][6],
    double SGS[N1 + 4][N2 + 4][N3 + 4],
    int      F[N1 + 4][N2 + 4][N3 + 4],
    double   ALPHA,
    double   DT,
    double   RI
) {
//  backup values of previous step

    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int k = K1 - 2; k <= K2 + 2; k ++) {
                UD[i][j][k][u] = U[i][j][k][u];
                UD[i][j][k][v] = U[i][j][k][v];
                UD[i][j][k][w] = U[i][j][k][w];
            }
        }
    }

//  pseudo velocity for fluid cells

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            for (int k = K1; k <= K2; k ++) {
                if (!(F[i][j][k] & FLUID)) {
                    continue;
                }
                double ABX, ABY, ABZ;
                double ADVX, ADVY, ADVZ;
                double VISX, VISY, VISZ;
                double UE1, UE2, UW1, UW2;
                double UN1, UN2, US1, US2;
                double UT1, UT2, UB1, UB2;
                double UC0U, UC0V, UC0W, UC0;
                double JUE, JUW, JUN, JUS, JUT, JUB;
                double JC0;
                double C1, C2, C3, C7, C8, C9;
                double NUE;
                double DUDX, DUDY, DUDZ;
                double DVDX, DVDY, DVDZ;
                double DWDX, DWDY, DWDZ;
                double K1X1, K2X2, K3X3;
                double DNDX, DNDY, DNDZ;
                double RSTX, RSTY, RSTZ;
                int    FE1, FE2, FW1, FW2;
                int    FN1, FN2, FS1, FS2;
                int    FT1, FT2, FB1, FB2;

                FE1  = FBU[i    ][j    ][k    ][x];
                FW1  = FBU[i - 1][j    ][k    ][x];
                FN1  = FBU[i    ][j    ][k    ][y];
                FS1  = FBU[i    ][j - 1][k    ][y];
                FT1  = FBU[i    ][j    ][k    ][z];
                FB1  = FBU[i    ][j    ][k - 1][z];
                FE2  = FBU[i + 1][j    ][k    ][x];
                FW2  = FBU[i - 2][j    ][k    ][x];
                FN2  = FBU[i    ][j + 1][k    ][y];
                FS2  = FBU[i    ][j - 2][k    ][y];
                FT2  = FBU[i    ][j    ][k + 1][z];
                FB2  = FBU[i    ][j    ][k - 2][z];
                K1X1 =  KX[i    ][j    ][k    ][x];
                K2X2 =  KX[i    ][j    ][k    ][y];
                K3X3 =  KX[i    ][j    ][k    ][z];
                UC0U =  UD[i    ][j    ][k    ][u];
                UC0V =  UD[i    ][j    ][k    ][v];
                UC0W =  UD[i    ][j    ][k    ][w];
                JUE  =  UU[i    ][j    ][k    ][u];
                JUN  =  UU[i    ][j    ][k    ][v];
                JUT  =  UU[i    ][j    ][k    ][w];
                JUW  =  UU[i - 1][j    ][k    ][u];
                JUS  =  UU[i    ][j - 1][k    ][v];
                JUB  =  UU[i    ][j    ][k - 1][w];
                JC0  =   J[i    ][j    ][k    ];
                C1   =   C[i    ][j    ][k    ][0];
                C2   =   C[i    ][j    ][k    ][1];
                C3   =   C[i    ][j    ][k    ][2];
                C7   =   C[i    ][j    ][k    ][3];
                C8   =   C[i    ][j    ][k    ][4];
                C9   =   C[i    ][j    ][k    ][5];
                NUE  = SGS[i    ][j    ][k    ];
                ABX  = abs(UC0U * K1X1);
                ABY  = abs(UC0V * K2X2);
                ABZ  = abs(UC0W * K3X3);

//  ADV and VIS terms for u[u]

                UC0  = UC0U;
                UE1  = UD[i + 1][j    ][k    ][u];
                UE2  = UD[i + 2][j    ][k    ][u];
                UW1  = UD[i - 1][j    ][k    ][u];
                UW2  = UD[i - 2][j    ][k    ][u];
                UN1  = UD[i    ][j + 1][k    ][u];
                UN2  = UD[i    ][j + 2][k    ][u];
                US1  = UD[i    ][j - 1][k    ][u];
                US2  = UD[i    ][j - 2][k    ][u];
                UT1  = UD[i    ][j    ][k + 1][u];
                UT2  = UD[i    ][j    ][k + 2][u];
                UB1  = UD[i    ][j    ][k - 1][u];
                UB2  = UD[i    ][j    ][k - 2][u];
                UE1  = (FE1)? (2 * BBU[i    ][j    ][k    ][x][u] - UC0) : (UE1);
                UW1  = (FW1)? (2 * BBU[i - 1][j    ][k    ][x][u] - UC0) : (UW1);
                UN1  = (FN1)? (2 * BBU[i    ][j    ][k    ][y][u] - UC0) : (UN1);
                US1  = (FS1)? (2 * BBU[i    ][j - 1][k    ][y][u] - UC0) : (US1);
                UT1  = (FT1)? (2 * BBU[i    ][j    ][k    ][z][u] - UC0) : (UT1);
                UB1  = (FB1)? (2 * BBU[i    ][j    ][k - 1][z][u] - UC0) : (UB1);
                UE2  = (FE2)? (2 * BBU[i + 1][j    ][k    ][x][u] - UE1) : (UE2);
                UW2  = (FW2)? (2 * BBU[i - 2][j    ][k    ][x][u] - UW1) : (UW2);
                UN2  = (FN2)? (2 * BBU[i    ][j + 1][k    ][y][u] - UN1) : (UN2);
                US2  = (FS2)? (2 * BBU[i    ][j - 2][k    ][y][u] - US1) : (US2);
                UT2  = (FT2)? (2 * BBU[i    ][j    ][k + 1][z][u] - UT1) : (UT2);
                UB2  = (FB2)? (2 * BBU[i    ][j    ][k - 2][z][u] - UB1) : (UB2);
                UE2  = (FE1)? (2 * BBU[i    ][j    ][k    ][x][u] - UW1) : (UE2);
                UW2  = (FW1)? (2 * BBU[i - 1][j    ][k    ][x][u] - UE1) : (UW2);
                UN2  = (FN1)? (2 * BBU[i    ][j    ][k    ][y][u] - US1) : (UN2);
                US2  = (FS1)? (2 * BBU[i    ][j - 1][k    ][y][u] - UN1) : (US2);
                UT2  = (FT1)? (2 * BBU[i    ][j    ][k    ][z][u] - UB1) : (UT2);
                UB2  = (FB1)? (2 * BBU[i    ][j    ][k - 1][z][u] - UT1) : (UB2);
                ADVX = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UT1, UT2, UB1, UB2, JUE, JUW, JUN, JUS, JUT, JUB, ABX, ABY, ABZ, JC0, ALPHA);
                VISX = vis(UC0, UE1, UW1, UN1, US1, UT1, UB1, NUE, RI, C1, C2, C3, C7, C8, C9);
                DUDX = 0.5 * (UE1 - UW1) * K1X1;
                DUDY = 0.5 * (UN1 - US1) * K2X2;
                DUDZ = 0.5 * (UT1 - UB1) * K3X3;

//  ADV and VIS terms for u[v]

                UC0  = UC0U;
                UE1  = UD[i + 1][j    ][k    ][v];
                UE2  = UD[i + 2][j    ][k    ][v];
                UW1  = UD[i - 1][j    ][k    ][v];
                UW2  = UD[i - 2][j    ][k    ][v];
                UN1  = UD[i    ][j + 1][k    ][v];
                UN2  = UD[i    ][j + 2][k    ][v];
                US1  = UD[i    ][j - 1][k    ][v];
                US2  = UD[i    ][j - 2][k    ][v];
                UT1  = UD[i    ][j    ][k + 1][v];
                UT2  = UD[i    ][j    ][k + 2][v];
                UB1  = UD[i    ][j    ][k - 1][v];
                UB2  = UD[i    ][j    ][k - 2][v];
                UE1  = (FE1)? (2 * BBU[i    ][j    ][k    ][x][v] - UC0) : (UE1);
                UW1  = (FW1)? (2 * BBU[i - 1][j    ][k    ][x][v] - UC0) : (UW1);
                UN1  = (FN1)? (2 * BBU[i    ][j    ][k    ][y][v] - UC0) : (UN1);
                US1  = (FS1)? (2 * BBU[i    ][j - 1][k    ][y][v] - UC0) : (US1);
                UT1  = (FT1)? (2 * BBU[i    ][j    ][k    ][z][v] - UC0) : (UT1);
                UB1  = (FB1)? (2 * BBU[i    ][j    ][k - 1][z][v] - UC0) : (UB1);
                UE2  = (FE2)? (2 * BBU[i + 1][j    ][k    ][x][v] - UE1) : (UE2);
                UW2  = (FW2)? (2 * BBU[i - 2][j    ][k    ][x][v] - UW1) : (UW2);
                UN2  = (FN2)? (2 * BBU[i    ][j + 1][k    ][y][v] - UN1) : (UN2);
                US2  = (FS2)? (2 * BBU[i    ][j - 2][k    ][y][v] - US1) : (US2);
                UT2  = (FT2)? (2 * BBU[i    ][j    ][k + 1][z][v] - UT1) : (UT2);
                UB2  = (FB2)? (2 * BBU[i    ][j    ][k - 2][z][v] - UB1) : (UB2);
                UE2  = (FE1)? (2 * BBU[i    ][j    ][k    ][x][v] - UW1) : (UE2);
                UW2  = (FW1)? (2 * BBU[i - 1][j    ][k    ][x][v] - UE1) : (UW2);
                UN2  = (FN1)? (2 * BBU[i    ][j    ][k    ][y][v] - US1) : (UN2);
                US2  = (FS1)? (2 * BBU[i    ][j - 1][k    ][y][v] - UN1) : (US2);
                UT2  = (FT1)? (2 * BBU[i    ][j    ][k    ][z][v] - UB1) : (UT2);
                UB2  = (FB1)? (2 * BBU[i    ][j    ][k - 1][z][v] - UT1) : (UB2);
                ADVY = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UT1, UT2, UB1, UB2, JUE, JUW, JUN, JUS, JUT, JUB, ABX, ABY, ABZ, JC0, ALPHA);
                VISY = vis(UC0, UE1, UW1, UN1, US1, UT1, UB1, NUE, RI, C1, C2, C3, C7, C8, C9);
                DVDX = 0.5 * (UE1 - UW1) * K1X1;
                DVDY = 0.5 * (UN1 - US1) * K2X2;
                DVDZ = 0.5 * (UT1 - UB1) * K3X3;

//  ADV and VIS terms for u[w]

                UC0  = UC0U;
                UE1  = UD[i + 1][j    ][k    ][w];
                UE2  = UD[i + 2][j    ][k    ][w];
                UW1  = UD[i - 1][j    ][k    ][w];
                UW2  = UD[i - 2][j    ][k    ][w];
                UN1  = UD[i    ][j + 1][k    ][w];
                UN2  = UD[i    ][j + 2][k    ][w];
                US1  = UD[i    ][j - 1][k    ][w];
                US2  = UD[i    ][j - 2][k    ][w];
                UT1  = UD[i    ][j    ][k + 1][w];
                UT2  = UD[i    ][j    ][k + 2][w];
                UB1  = UD[i    ][j    ][k - 1][w];
                UB2  = UD[i    ][j    ][k - 2][w];
                UE1  = (FE1)? (2 * BBU[i    ][j    ][k    ][x][w] - UC0) : (UE1);
                UW1  = (FW1)? (2 * BBU[i - 1][j    ][k    ][x][w] - UC0) : (UW1);
                UN1  = (FN1)? (2 * BBU[i    ][j    ][k    ][y][w] - UC0) : (UN1);
                US1  = (FS1)? (2 * BBU[i    ][j - 1][k    ][y][w] - UC0) : (US1);
                UT1  = (FT1)? (2 * BBU[i    ][j    ][k    ][z][w] - UC0) : (UT1);
                UB1  = (FB1)? (2 * BBU[i    ][j    ][k - 1][z][w] - UC0) : (UB1);
                UE2  = (FE2)? (2 * BBU[i + 1][j    ][k    ][x][w] - UE1) : (UE2);
                UW2  = (FW2)? (2 * BBU[i - 2][j    ][k    ][x][w] - UW1) : (UW2);
                UN2  = (FN2)? (2 * BBU[i    ][j + 1][k    ][y][w] - UN1) : (UN2);
                US2  = (FS2)? (2 * BBU[i    ][j - 2][k    ][y][w] - US1) : (US2);
                UT2  = (FT2)? (2 * BBU[i    ][j    ][k + 1][z][w] - UT1) : (UT2);
                UB2  = (FB2)? (2 * BBU[i    ][j    ][k - 2][z][w] - UB1) : (UB2);
                UE2  = (FE1)? (2 * BBU[i    ][j    ][k    ][x][w] - UW1) : (UE2);
                UW2  = (FW1)? (2 * BBU[i - 1][j    ][k    ][x][w] - UE1) : (UW2);
                UN2  = (FN1)? (2 * BBU[i    ][j    ][k    ][y][w] - US1) : (UN2);
                US2  = (FS1)? (2 * BBU[i    ][j - 1][k    ][y][w] - UN1) : (US2);
                UT2  = (FT1)? (2 * BBU[i    ][j    ][k    ][z][w] - UB1) : (UT2);
                UB2  = (FB1)? (2 * BBU[i    ][j    ][k - 1][z][w] - UT1) : (UB2);
                ADVZ = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UT1, UT2, UB1, UB2, JUE, JUW, JUN, JUS, JUT, JUB, ABX, ABY, ABZ, JC0, ALPHA);
                VISZ = vis(UC0, UE1, UW1, UN1, US1, UT1, UB1, NUE, RI, C1, C2, C3, C7, C8, C9);
                DWDX = 0.5 * (UE1 - UW1) * K1X1;
                DWDY = 0.5 * (UN1 - US1) * K2X2;
                DWDZ = 0.5 * (UT1 - UB1) * K3X3;

//  SGS vortex viscousity terms 

                DNDX = K1X1 * 0.5 * (SGS[i + 1][j    ][k    ] - SGS[i - 1][j    ][k    ]);
                DNDY = K2X2 * 0.5 * (SGS[i    ][j + 1][k    ] - SGS[i    ][j - 1][k    ]);
                DNDZ = K3X3 * 0.5 * (SGS[i    ][j    ][k + 1] - SGS[i    ][j    ][k - 1]);
                RSTX = (DUDX + DUDX) * DNDX + (DUDY + DVDX) * DNDY + (DUDZ + DWDX) * DNDZ;
                RSTY = (DVDX + DUDY) * DNDX + (DVDY + DVDY) * DNDY + (DVDZ + DWDY) * DNDZ;
                RSTZ = (DWDX + DUDZ) * DNDX + (DWDY + DVDZ) * DNDY + (DWDZ + DWDZ) * DNDZ;

// pseudo velocity vector

                U[i][j][k][u] = UC0U + DT * (- ADVX + VISX + RSTX);
                U[i][j][k][v] = UC0V + DT * (- ADVY + VISY + RSTY);
                U[i][j][k][w] = UC0W + DT * (- ADVZ + VISZ + RSTZ);
            }
        }
    }
}