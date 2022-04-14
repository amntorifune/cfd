#include "topo.h"
#include "flag.h"
#include <math.h>

double max(double a, double b) {
    return ((a > b)? a : b);
}

double min(double a, double b) {
    return ((a < b)? a : b);
}

double minmod(double r) {
    return max(0.0, min(1.0, r));
}

double leer(double r) {
    return (r + fabs(r)) / (1.0 + fabs(r));
}

double superbee(double r) {
    return max(0.0, max(min(2.0 * r, 1.0), min(r, 2.0)));
}

double mc(double r) {
    return max(0.0, min(2.0 * r, min(0.5 * (1.0 + r), 2.0)));
}

double rr(double d1, double d2) {
    return (d1 + copysign(1.0e-24, d1)) / (d2 + copysign(1.0e-24, d2));
}

double muscl(
    double UC0,
    double UE1,
    double UE2,
    double UW1,
    double UW2,
    double UN1,
    double UN2,
    double US1,
    double US2,
    double UUE,
    double UUW,
    double UUN,
    double UUS,
    double J0J
) {
    double D1, D2, D3, D4;
    double URR, URL, ULR, ULL;
    double FR, FL;
    double A1, A2;

    double kappa = 1.0 / 3.0;
    double beta  = (3.0 - kappa) / (1.0 - kappa);
    double m1    = 1.0 - kappa;
    double m2    = 1.0 + kappa;

    D4  = UE2 - UE1;
    D3  = UE1 - UC0;
    D2  = UC0 - UW1;
    D1  = UW1 - UW2;
    URR = UE1 - 0.25 * mc(rr(D3, D4)) * (m1 * D4 + m2 * D3);
    URL = UC0 + 0.25 * mc(rr(D2, D3)) * (m1 * D2 + m2 * D3);
    ULR = UC0 - 0.25 * mc(rr(D2, D3)) * (m1 * D3 + m2 * D2);
    ULL = UW1 + 0.25 * mc(rr(D1, D2)) * (m1 * D1 + m2 * D2);
    FR  = 0.5 * (UUE * (URR + URL) - fabs(UUE) * (URR - URL));
    FL  = 0.5 * (UUW * (ULR + ULL) - fabs(UUW) * (ULR - ULL));
    A1  = (FR - FL) / J0J;

    D4  = UN2 - UN1;
    D3  = UN1 - UC0;
    D2  = UC0 - US1;
    D1  = US1 - US2;
    URR = UN1 - 0.25 * mc(rr(D3, D4)) * (m1 * D4 + m2 * D3);
    URL = UC0 + 0.25 * mc(rr(D2, D3)) * (m1 * D2 + m2 * D3);
    ULR = UC0 - 0.25 * mc(rr(D2, D3)) * (m1 * D3 + m2 * D2);
    ULL = US1 + 0.25 * mc(rr(D1, D2)) * (m1 * D1 + m2 * D2);
    FR  = 0.5 * (UUN * (URR + URL) - fabs(UUN) * (URR - URL));
    FL  = 0.5 * (UUS * (ULR + ULL) - fabs(UUS) * (ULR - ULL));
    A2  = (FR - FL) / J0J;

    return A1 + A2;
}

double ffvcmuscl(
    double UC0,
    double UE1,
    double UE2,
    double UW1,
    double UW2,
    double UN1,
    double UN2,
    double US1,
    double US2,
    double UUE,
    double UUW,
    double UUN,
    double UUS,
    double J0J,
    int    EPS
) {
    double D1, D2, D3, D4;
    double S1, S2, S3, S4;
    double G1, G2, G3, G4, G5, G6;
    double URR, URL, ULR, ULL;
    double FR, FL;
    double A1, A2;

    double kappa = 1.0 / 3.0;
    double beta  = (3.0 - kappa) / (1.0 - kappa);
    double m1    = 1.0 - kappa;
    double m2    = 1.0 + kappa;

    D4  = UE2 - UE1;
    D3  = UE1 - UC0;
    D2  = UC0 - UW1;
    D1  = UW1 - UW2;
    S4  = copysign(1.0, D4);
    S3  = copysign(1.0, D3);
    S2  = copysign(1.0, D2);
    S1  = copysign(1.0, D1);
    G6  = S4 * max(0.0, min(fabs(D4), S4 * beta * D3));
    G5  = S3 * max(0.0, min(fabs(D3), S3 * beta * D4));
    G4  = S3 * max(0.0, min(fabs(D3), S3 * beta * D2));
    G3  = S2 * max(0.0, min(fabs(D2), S2 * beta * D3));
    G2  = S2 * max(0.0, min(fabs(D2), S2 * beta * D1));
    G1  = S1 * max(0.0, min(fabs(D1), S1 * beta * D2));
    URR = UE1 - 0.25 * (m1 * G6 + m2 * G5) * EPS;
    URL = UC0 + 0.25 * (m1 * G3 + m2 * G4) * EPS;
    ULR = UC0 - 0.25 * (m1 * G4 + m2 * G3) * EPS;
    ULL = UW1 + 0.25 * (m1 * G1 + m2 * G2) * EPS;
    FR  = 0.5 * (UUE * (URR + URL) - fabs(UUE) * (URR - URL));
    FL  = 0.5 * (UUW * (ULR + ULL) - fabs(UUW) * (ULR - ULL));
    A1  = (FR - FL) / J0J;

    D4  = UN2 - UN1;
    D3  = UN1 - UC0;
    D2  = UC0 - US1;
    D1  = US1 - US2;
    S4  = copysign(1.0, D4);
    S3  = copysign(1.0, D3);
    S2  = copysign(1.0, D2);
    S1  = copysign(1.0, D1);
    G6  = S4 * max(0.0, min(fabs(D4), S4 * beta * D3));
    G5  = S3 * max(0.0, min(fabs(D3), S3 * beta * D4));
    G4  = S3 * max(0.0, min(fabs(D3), S3 * beta * D2));
    G3  = S2 * max(0.0, min(fabs(D2), S2 * beta * D3));
    G2  = S2 * max(0.0, min(fabs(D2), S2 * beta * D1));
    G1  = S1 * max(0.0, min(fabs(D1), S1 * beta * D2));
    URR = UN1 - 0.25 * (m1 * G6 + m2 * G5) * EPS;
    URL = UC0 + 0.25 * (m1 * G3 + m2 * G4) * EPS;
    ULR = UC0 - 0.25 * (m1 * G4 + m2 * G3) * EPS;
    ULL = US1 + 0.25 * (m1 * G1 + m2 * G2) * EPS;
    FR  = 0.5 * (UUN * (URR + URL) - fabs(UUN) * (URR - URL));
    FL  = 0.5 * (UUS * (ULR + ULL) - fabs(UUS) * (ULR - ULL));
    A2  = (FR - FL) / J0J;

    return A1 + A2;
}

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
    double UUE,
    double UUW,
    double UUN,
    double UUS,
    double ABX,
    double ABY,
    double J0J,
    double ALP
) {
    double A1, A2;

    A1  = UUE * (- UE2 + 27 * UE1 - 27 * UC0 + UW1);
    A1 += UUW * (- UE1 + 27 * UC0 - 27 * UW1 + UW2);
    A1  = A1 / (48 * J0J);
    A1 += ALP * ABX * (UE2 - 4 * UE1 + 6 * UC0 - 4 * UW1 + UW2);

    A2  = UUN * (- UN2 + 27 * UN1 - 27 * UC0 + US1);
    A2 += UUS * (- UN1 + 27 * UC0 - 27 * US1 + US2);
    A2  = A2 / (48 * J0J);
    A2 += ALP * ABY * (UN2 - 4 * UN1 + 6 * UC0 - 4 * US1 + US2);

    return A1 + A2;
}

double vis(
    double UC0,
    double UE1,
    double UW1,
    double UN1,
    double US1,
    double NUE,
    double RI,
    double C1,
    double C2,
    double C7,
    double C8
) {
    double LAP, DK1, DK2, D2K1, D2K2;

    DK1  = 0.5 * (UE1 - UW1);
    DK2  = 0.5 * (UN1 - US1);
    D2K1 = UE1 - 2 * UC0 + UW1;
    D2K2 = UN1 - 2 * UC0 + US1;
    LAP  = C1 * D2K1 + C2 * D2K2 - C7 * DK1 - C8 * DK2;

    return (NUE + RI) * LAP;
}

void pseudo(
    double   U[NX + 4][NY + 4][2],
    double  UA[NX + 4][NY + 4][2],
    double  UU[NX + 4][NY + 4][2],
    int    FFU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double SGS[NX + 4][NY + 4],
    double  KX[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   ALPHA,
    double   DT,
    double   RI
) {
    #pragma acc kernels loop independent collapse(2) present(U, UA, UU, FFU, BBU, SGS, KX, J, C, F) copyin(ALPHA, DT, RI)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (F[i][j] == FLUID) {
                double ABX , ABY ;
                double ADV1, ADV2;
                double VISX, VISY;
                double UE1, UE2, UW1, UW2;
                double UN1, UN2, US1, US2;
                double U0U, U0V, UC0;
                double UUE, UUW, UUN, UUS;
                int    FE1, FE2, FW1, FW2;
                int    FN1, FN2, FS1, FS2;
                double J0J;
                double C1, C2, C7, C8;
                double NUE;
                double DUDX, DUDY;
                double DVDX, DVDY;
                double K1X1, K2X2;
                double DNDX, DNDY;
                double RSTX, RSTY;

                FE1  = FFU[i    ][j    ][0];
                FW1  = FFU[i - 1][j    ][0];
                FN1  = FFU[i    ][j    ][1];
                FS1  = FFU[i    ][j - 1][1];
                FE2  = FFU[i + 1][j    ][0];
                FW2  = FFU[i - 2][j    ][0];
                FN2  = FFU[i    ][j + 1][1];
                FS2  = FFU[i    ][j - 2][1];
                K1X1 =  KX[i    ][j    ][0];
                K2X2 =  KX[i    ][j    ][1];
                U0U  =   U[i    ][j    ][0];
                U0V  =   U[i    ][j    ][1];
                UUE  =  UU[i    ][j    ][0];
                UUN  =  UU[i    ][j    ][1];
                UUW  =  UU[i - 1][j    ][0];
                UUS  =  UU[i    ][j - 1][1];
                J0J  =   J[i    ][j    ];
                C1   =   C[i    ][j    ][0];
                C2   =   C[i    ][j    ][1];
                C7   =   C[i    ][j    ][2];
                C8   =   C[i    ][j    ][3];
                NUE  = SGS[i][j];
                ABX  = fabs(U0U * K1X1);
                ABY  = fabs(U0V * K2X2);

//  ADV and VIS terms for u[u]

                UC0  = U0U;
                UE1  = U[i + 1][j    ][0];
                UE2  = U[i + 2][j    ][0];
                UW1  = U[i - 1][j    ][0];
                UW2  = U[i - 2][j    ][0];
                UN1  = U[i    ][j + 1][0];
                UN2  = U[i    ][j + 2][0];
                US1  = U[i    ][j - 1][0];
                US2  = U[i    ][j - 2][0];
                UE1  = (FE1 == BOUNDARY)? (2 * BBU[i    ][j    ][0][0] - UC0) : (UE1);
                UW1  = (FW1 == BOUNDARY)? (2 * BBU[i - 1][j    ][0][0] - UC0) : (UW1);
                UN1  = (FN1 == BOUNDARY)? (2 * BBU[i    ][j    ][1][0] - UC0) : (UN1);
                US1  = (FS1 == BOUNDARY)? (2 * BBU[i    ][j - 1][1][0] - UC0) : (US1);
                UE2  = (FE2 == BOUNDARY)? (2 * BBU[i + 1][j    ][0][0] - UE1) : (UE2);
                UW2  = (FW2 == BOUNDARY)? (2 * BBU[i - 2][j    ][0][0] - UW1) : (UW2);
                UN2  = (FN2 == BOUNDARY)? (2 * BBU[i    ][j + 1][1][0] - UN1) : (UN2);
                US2  = (FS2 == BOUNDARY)? (2 * BBU[i    ][j - 2][1][0] - US1) : (US2);
                UE2  = (FE1 == BOUNDARY)? (2 * BBU[i    ][j    ][0][0] - UW1) : (UE2);
                UW2  = (FW1 == BOUNDARY)? (2 * BBU[i - 1][j    ][0][0] - UE1) : (UW2);
                UN2  = (FN1 == BOUNDARY)? (2 * BBU[i    ][j    ][1][0] - US1) : (UN2);
                US2  = (FS1 == BOUNDARY)? (2 * BBU[i    ][j - 1][1][0] - UN1) : (US2);
                // ADV1 = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UUE, UUW, UUN, UUS, ABX, ABY, J0J, ALPHA);
                ADV1 = ffvcmuscl(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UUE, UUW, UUN, UUS, J0J, 1);
                VISX = vis(UC0, UE1, UW1, UN1, US1, NUE, RI, C1, C2, C7, C8);
                DUDX = 0.5 * (UE1 - UW1) * K1X1;
                DUDY = 0.5 * (UN1 - US1) * K2X2;

//  ADV and VIS terms for u[v]

                UC0  = U0V;
                UE1  = U[i + 1][j    ][1];
                UE2  = U[i + 2][j    ][1];
                UW1  = U[i - 1][j    ][1];
                UW2  = U[i - 2][j    ][1];
                UN1  = U[i    ][j + 1][1];
                UN2  = U[i    ][j + 2][1];
                US1  = U[i    ][j - 1][1];
                US2  = U[i    ][j - 2][1];
                UE1  = (FE1 == BOUNDARY)? (2 * BBU[i    ][j    ][0][1] - UC0) : (UE1);
                UW1  = (FW1 == BOUNDARY)? (2 * BBU[i - 1][j    ][0][1] - UC0) : (UW1);
                UN1  = (FN1 == BOUNDARY)? (2 * BBU[i    ][j    ][1][1] - UC0) : (UN1);
                US1  = (FS1 == BOUNDARY)? (2 * BBU[i    ][j - 1][1][1] - UC0) : (US1);
                UE2  = (FE2 == BOUNDARY)? (2 * BBU[i + 1][j    ][0][1] - UE1) : (UE2);
                UW2  = (FW2 == BOUNDARY)? (2 * BBU[i - 2][j    ][0][1] - UW1) : (UW2);
                UN2  = (FN2 == BOUNDARY)? (2 * BBU[i    ][j + 1][1][1] - UN1) : (UN2);
                US2  = (FS2 == BOUNDARY)? (2 * BBU[i    ][j - 2][1][1] - US1) : (US2);
                UE2  = (FE1 == BOUNDARY)? (2 * BBU[i    ][j    ][0][1] - UW1) : (UE2);
                UW2  = (FW1 == BOUNDARY)? (2 * BBU[i - 1][j    ][0][1] - UE1) : (UW2);
                UN2  = (FN1 == BOUNDARY)? (2 * BBU[i    ][j    ][1][1] - US1) : (UN2);
                US2  = (FS1 == BOUNDARY)? (2 * BBU[i    ][j - 1][1][1] - UN1) : (US2);
                // ADV2 = adv(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UUE, UUW, UUN, UUS, ABX, ABY, J0J, ALPHA);
                ADV2 = ffvcmuscl(UC0, UE1, UE2, UW1, UW2, UN1, UN2, US1, US2, UUE, UUW, UUN, UUS, J0J, 1);
                VISY = vis(UC0, UE1, UW1, UN1, US1, NUE, RI, C1, C2, C7, C8);
                DVDX = 0.5 * (UE1 - UW1) * K1X1;
                DVDY = 0.5 * (UN1 - US1) * K2X2;

//  SGS vortex viscousity terms 

                DNDX = K1X1 * 0.5 * (SGS[i + 1][j    ] - SGS[i - 1][j    ]);
                DNDY = K2X2 * 0.5 * (SGS[i    ][j + 1] - SGS[i    ][j - 1]);
                RSTX = (DUDX + DUDX) * DNDX + (DUDY + DVDX) * DNDY;
                RSTY = (DVDX + DUDY) * DNDX + (DVDY + DVDY) * DNDY;

// intermediate velocity vector

                UA[i][j][0] = U0U + DT * (- ADV1 + VISX);
                UA[i][j][1] = U0V + DT * (- ADV2 + VISY);
            }
        }
    }
}

void correction1(
    double   U[NX + 4][NY + 4][2],
    double  UD[NX + 4][NY + 4][2],
    double   P[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double  KX[NX + 4][NY + 4][2],
    int      F[NX + 4][NY + 4],
    double   DT
) {
    #pragma acc kernels loop independent collapse(2) present(U, UD, P, FFP, BBP, KX, F) copyin(DT)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (F[i][j] == FLUID) {
                double D1, D2;
                double PC0, PE1, PW1, PN1, PS1;
                int    FE1, FW1, FN1, FS1;
                double K1X1, K2X2;

                PC0  =   P[i    ][j    ];
                PE1  =   P[i + 1][j    ];
                PW1  =   P[i - 1][j    ];
                PN1  =   P[i    ][j + 1];
                PS1  =   P[i    ][j - 1];
                K1X1 =  KX[i    ][j    ][0];
                K2X2 =  KX[i    ][j    ][1];
                FE1  = FFP[i    ][j    ][0];
                FW1  = FFP[i - 1][j    ][0];
                FN1  = FFP[i    ][j    ][1];
                FS1  = FFP[i    ][j - 1][1];
                PE1  = (FE1 == BOUNDARY)? (2 * BBP[i    ][j    ][0] - PC0) : (PE1);
                PW1  = (FW1 == BOUNDARY)? (2 * BBP[i - 1][j    ][0] - PC0) : (PW1);
                PN1  = (FN1 == BOUNDARY)? (2 * BBP[i    ][j    ][1] - PC0) : (PN1);
                PS1  = (FS1 == BOUNDARY)? (2 * BBP[i    ][j - 1][1] - PC0) : (PS1);

                D1   = K1X1 * 0.5 * (PE1 - PW1);
                D2   = K2X2 * 0.5 * (PN1 - PS1);

                U[i][j][0] = UD[i][j][0] - DT * D1;
                U[i][j][1] = UD[i][j][1] - DT * D2;
            }
            
        }
    }
}

void correction2(
    double  UU[NX + 4][NY + 4][2],
    double UUD[NX + 4][NY + 4][2],
    int    FFU[NX + 4][NY + 4][2],
    double BBU[NX + 4][NY + 4][2][2],
    double   P[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double   X[NX + 4][NY + 4][2],
    double  KX[NX + 4][NY + 4][2],
    double   G[NX + 4][NY + 4][2],
    double   DT
) {
    #pragma acc kernels loop independent collapse(2) present(UU, UUD, FFU, BBU, P, FFP, BBP, X, KX, G) copyin(DT)
    for (int i = I1 - 1; i <= I2; i ++) {
        for (int j = J1 - 1; j <= J2; j ++) {
            int    FF;
            double G0, G1;
            double P0, P1;
            double JF, X1K1, X2K2, K1X1, K2X2;
            double DF, UF, UUF;

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
                P0   =              P[i    ][j];
                P1   =              P[i + 1][j];
                G0   =              G[i    ][j][0];
                G1   =              G[i + 1][j][0];
                FF   =            FFP[i    ][j][0];
                P1   = (FF == BOUNDARY)? (2 * BBP[i    ][j][0] - P0) : (P1);
                DF   = 0.5 * (G0 + G1) * (P1 - P0);
                UUF  = UUD[i][j][0] - DT * DF;
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
                P0   =              P[i][j    ];
                P1   =              P[i][j + 1];
                G0   =              G[i][j    ][1];
                G1   =              G[i][j + 1][1];
                FF   =            FFP[i][j    ][1];
                P1   = (FF == BOUNDARY)? (2 * BBP[i][j    ][1] - P0) : (P1);
                DF   = 0.5 * (G0 + G1) * (P1 - P0);
                UUF  = UUD[i][j][1] - DT * DF;
            }
            UU[i][j][1] = UUF;
        }
    }
}

void correction3(
    double  P[NX + 4][NY + 4],
    double Pi[NX + 4][NY + 4]
) {
    #pragma acc kernels loop independent collapse(2) present(P, Pi)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            P[i][j] += Pi[i][j];
        }
    }
}
