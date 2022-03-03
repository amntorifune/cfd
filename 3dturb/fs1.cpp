#include "topo.h"
#include <math.h>

//  calculate advection term

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
    double UUE,
    double UUW,
    double UUN,
    double UUS,
    double UUT,
    double UUB,
    double ABX,
    double ABY,
    double ABZ,
    double JC0,
    double ALPHA
) {
    double ADVX, ADVY, ADVZ;

    ADVX = UUE * (- UE2 + 27 * UE1 - 27 * UC0 + UW1);
    ADVX = UUW * (- UE1 + 27 * UC0 - 27 * UW1 + UW2) + ADVX;
    ADVX = ADVX / (48.0 * JC0);
    ADVX = ADVX + ALPHA * ABX * (UE2 - 4 * UE1 + 6 * UC0 - 4 * UW1 + UW2);

    ADVY = UUN * (- UN2 + 27 * UN1 - 27 * UC0 + US1);
    ADVY = UUS * (- UN1 + 27 * UC0 - 27 * US1 + US2) + ADVY;
    ADVY = ADVY / (48.0 * JC0);
    ADVY = ADVY + ALPHA * ABY * (UN2 - 4 * UN1 + 6 * UC0 - 4 * US1 + US2);

    ADVZ = UUT * (- UT2 + 27 * UT1 - 27 * UC0 + UB1);
    ADVZ = UUB * (- UT1 + 27 * UC0 - 27 * UB1 + UB2) + ADVZ;
    ADVZ = ADVZ / (48.0 * JC0);
    ADVZ = ADVZ + ALPHA * ABZ * (UT2 - 4 * UT1 + 6 * UC0 - 4 * UB1 + UB2);

    return (ADVX + ADVY + ADVZ);
}

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

//  calculate intermediate velocity

void fs1(
    double   U[NX + 2][NY + 2][NZ + 2][3],
    double  UD[NX + 2][NY + 2][NZ + 2][3],
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double SGS[NX + 2][NY + 2][NZ + 2],
    double  KX[NX + 2][NY + 2][NZ + 2][3],
    double   J[NX + 2][NY + 2][NZ + 2],
    double   C[NX + 2][NY + 2][NZ + 2][6],
    double   ALPHA,
    double   DT,
    double   RI
) {
    double ABX, ABY, ABZ; // abs of contravariant velocity at cell center
    double ADVX, ADVY, ADVZ;    // advection terms of each component
    double VISX, VISY, VISZ;    // viscoucity stress in each direction
    double UE1, UE2, UW1, UW2;
    double UN1, UN2, US1, US2;
    double UT1, UT2, UB1, UB2;
    double UC0U, UC0V, UC0W, UC0;
    double UUE, UUW, UUN, UUS, UUT, UUB;
    double JC0;
    double C1, C2, C3, C7, C8, C9;
    double NUE;
    double DUDX, DUDY, DUDZ;
    double DVDX, DVDY, DVDZ;
    double DWDX, DWDY, DWDZ;
    double K1X1, K2X2, K3X3;
    double DNDX, DNDY, DNDZ;
    double RSTX, RSTY, RSTZ;

//  store the values of previous step

    for (int i = 0; i <= NX + 1; i ++) {
        for (int j = 0; j <= NY + 1; j ++) {
            for (int k = 0; k <= NZ + 1; k ++) {
                UD[i][j][k][0] = U[i][j][k][0];
                UD[i][j][k][1] = U[i][j][k][1];
                UD[i][j][k][2] = U[i][j][k][2];
            }
        }
    }

// calculate values for inner points, i.e. [2:NX)[2:NY)[2:NZ)

    for (int i = 2; i <= NX - 1; i ++) {
        for (int j = 2; j <= NY - 1; j ++) {
            for (int k = 2; k <= NZ - 1; k ++) {
                K1X1 =  KX[i][j][k][0];
                K2X2 =  KX[i][j][k][1];
                K3X3 =  KX[i][j][k][2];
                UC0U =  UD[i][j][k][0];
                UC0V =  UD[i][j][k][1];
                UC0W =  UD[i][j][k][2];
                UUE  =  UU[i    ][j    ][k    ][0];
                UUN  =  UU[i    ][j    ][k    ][1];
                UUT  =  UU[i    ][j    ][k    ][2];
                UUW  =  UU[i - 1][j    ][k    ][0];
                UUS  =  UU[i    ][j - 1][k    ][1];
                UUB  =  UU[i    ][j    ][k - 1][2];
                JC0  =   J[i][j][k];
                C1   =   C[i][j][k][0];
                C2   =   C[i][j][k][1];
                C3   =   C[i][j][k][2];
                C7   =   C[i][j][k][3];
                C8   =   C[i][j][k][4];
                C9   =   C[i][j][k][5];
                NUE  = SGS[i][j][k];
                ABX  = abs(UC0U * K1X1);
                ABY  = abs(UC0V * K2X2);
                ABZ  = abs(UC0W * K3X3);

//  ADV and VIS terms for u[u]

                UC0  = UC0U;
                UE1  = UD[i + 1][j    ][k    ][0];
                UE2  = UD[i + 2][j    ][k    ][0];
                UW1  = UD[i - 1][j    ][k    ][0];
                UW2  = UD[i - 2][j    ][k    ][0];
                UN1  = UD[i    ][j + 1][k    ][0];
                UN2  = UD[i    ][j + 2][k    ][0];
                US1  = UD[i    ][j - 1][k    ][0];
                US2  = UD[i    ][j - 2][k    ][0];
                UT1  = UD[i    ][j    ][k + 1][0];
                UT2  = UD[i    ][j    ][k + 2][0];
                UB1  = UD[i    ][j    ][k - 1][0];
                UB2  = UD[i    ][j    ][k - 2][0];
                ADVX = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UT1, UT2, UB1, UB2, UUE, UUW, UUN, UUS, UUT, UUB, ABX, ABY, ABZ, JC0, ALPHA);
                VISX = vis(UC0, UE1, UW1, UN1, US1, UT1, UB1, NUE, RI, C1, C2, C3, C7, C8, C9);
                DUDX = 0.5 * (UE1 - UW1) * K1X1;
                DUDY = 0.5 * (UN1 - US1) * K2X2;
                DUDZ = 0.5 * (UT1 - UB1) * K3X3;

//  ADV and VIS terms for u[v]

                UC0  = UC0V;
                UE1  = UD[i + 1][j    ][k    ][1];
                UE2  = UD[i + 2][j    ][k    ][1];
                UW1  = UD[i - 1][j    ][k    ][1];
                UW2  = UD[i - 2][j    ][k    ][1];
                UN1  = UD[i    ][j + 1][k    ][1];
                UN2  = UD[i    ][j + 2][k    ][1];
                US1  = UD[i    ][j - 1][k    ][1];
                US2  = UD[i    ][j - 2][k    ][1];
                UT1  = UD[i    ][j    ][k + 1][1];
                UT2  = UD[i    ][j    ][k + 2][1];
                UB1  = UD[i    ][j    ][k - 1][1];
                UB2  = UD[i    ][j    ][k - 2][1];
                ADVY = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UT1, UT2, UB1, UB2, UUE, UUW, UUN, UUS, UUT, UUB, ABX, ABY, ABZ, JC0, ALPHA);
                VISY = vis(UC0, UE1, UW1, UN1, US1, UT1, UB1, NUE, RI, C1, C2, C3, C7, C8, C9);
                DVDX = 0.5 * (UE1 - UW1) * K1X1;
                DVDY = 0.5 * (UN1 - US1) * K2X2;
                DVDZ = 0.5 * (UT1 - UB1) * K3X3;

//  ADV and VIS terms for u[w]

                UC0  = UC0W;
                UE1  = UD[i + 1][j    ][k    ][2];
                UE2  = UD[i + 2][j    ][k    ][2];
                UW1  = UD[i - 1][j    ][k    ][2];
                UW2  = UD[i - 2][j    ][k    ][2];
                UN1  = UD[i    ][j + 1][k    ][2];
                UN2  = UD[i    ][j + 2][k    ][2];
                US1  = UD[i    ][j - 1][k    ][2];
                US2  = UD[i    ][j - 2][k    ][2];
                UT1  = UD[i    ][j    ][k + 1][2];
                UT2  = UD[i    ][j    ][k + 2][2];
                UB1  = UD[i    ][j    ][k - 1][2];
                UB2  = UD[i    ][j    ][k - 2][2];
                ADVZ = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UT1, UT2, UB1, UB2, UUE, UUW, UUN, UUS, UUT, UUB, ABX, ABY, ABZ, JC0, ALPHA);
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

// intermediate velocity vector

                U[i][j][k][0] = UC0U + DT * (- ADVX + VISX + RSTX);
                U[i][j][k][1] = UC0V + DT * (- ADVY + VISY + RSTY);
                U[i][j][k][2] = UC0W + DT * (- ADVZ + VISZ + RSTZ);

            }
        }
    }
}