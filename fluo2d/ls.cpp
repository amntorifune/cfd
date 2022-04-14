#include "topo.h"
#include "flag.h"
#include <math.h>

void hsmac(
    double   P[NX + 4][NY + 4],
    double PSI[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   DT
) {
    #pragma acc kernels loop independent collapse(2) present(P, PSI, C, F) copyin(DT)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (F[i][j] == FLUID) {
                double AC0, RHS;
                RHS = PSI[i][j] / DT;
                AC0 = - 2.0 * (C[i][j][0] + C[i][j][1]);
                P[i][j] = RHS / AC0;
            }
        }
    }
}

void sor(
    double   P[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double PSI[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   OMEGA,
    double   DT,
    double  &R
) {
    double ERR = 0;
    int    CNT = 0;

    #pragma acc kernels loop independent collapse(2) reduction(+:ERR, CNT) present(P, FFP, BBP, PSI, C, F) copyin(OMEGA, DT) copy(ERR, CNT)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if ((F[i][j] == FLUID) && ((i + j) % 2 == 0)) {
                double AC0, AE1, AW1, AN1, AS1;
                double PC0, PE1, PW1, PN1, PS1;
                int    FE1, FW1, FN1, FS1;
                double C1, C2, C7, C8;
                double RHS, RP;

                RHS = PSI[i    ][j    ] / DT;
                C1  =   C[i    ][j    ][0];
                C2  =   C[i    ][j    ][1];
                C7  =   C[i    ][j    ][2];
                C8  =   C[i    ][j    ][3];
                PC0 =   P[i    ][j    ];
                PE1 =   P[i + 1][j    ];
                PW1 =   P[i - 1][j    ];
                PN1 =   P[i    ][j + 1];
                PS1 =   P[i    ][j - 1];
                FE1 = FFP[i    ][j    ][0];
                FW1 = FFP[i - 1][j    ][0];
                FN1 = FFP[i    ][j    ][1];
                FS1 = FFP[i    ][j - 1][1];
                PE1 = (FE1 == BOUNDARY)? (2 * BBP[i    ][j    ][0] - PC0) : (PE1);
                PW1 = (FW1 == BOUNDARY)? (2 * BBP[i - 1][j    ][0] - PC0) : (PW1);
                PN1 = (FN1 == BOUNDARY)? (2 * BBP[i    ][j    ][1] - PC0) : (PN1);
                PS1 = (FS1 == BOUNDARY)? (2 * BBP[i    ][j - 1][1] - PC0) : (PS1);

                AC0 = - 2.0 * (C1 + C2);
                AE1 = C1 - 0.5 * C7;
                AW1 = C1 + 0.5 * C7;
                AN1 = C2 - 0.5 * C8;
                AS1 = C2 + 0.5 * C8;
                RP  = (RHS - AE1 * PE1 - AW1 * PW1 - AN1 * PN1 - AS1 * PS1) / AC0 - PC0;
                ERR = ERR + RP * RP;

                P[i][j] = PC0 + OMEGA * RP;
                CNT += 1;
            }
            
        }
    }

    #pragma acc kernels loop independent collapse(2) reduction(+:ERR, CNT) present(P, FFP, BBP, PSI, C, F) copyin(OMEGA, DT) copy(ERR, CNT)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if ((F[i][j] == FLUID) && ((i + j) % 2 == 1)) {
                double AC0, AE1, AW1, AN1, AS1;
                double PC0, PE1, PW1, PN1, PS1;
                int    FE1, FW1, FN1, FS1;
                double C1, C2, C7, C8;
                double RHS, RP;

                RHS = PSI[i    ][j    ] / DT;
                C1  =   C[i    ][j    ][0];
                C2  =   C[i    ][j    ][1];
                C7  =   C[i    ][j    ][2];
                C8  =   C[i    ][j    ][3];
                PC0 =   P[i    ][j    ];
                PE1 =   P[i + 1][j    ];
                PW1 =   P[i - 1][j    ];
                PN1 =   P[i    ][j + 1];
                PS1 =   P[i    ][j - 1];
                FE1 = FFP[i    ][j    ][0];
                FW1 = FFP[i - 1][j    ][0];
                FN1 = FFP[i    ][j    ][1];
                FS1 = FFP[i    ][j - 1][1];
                PE1 = (FE1 == BOUNDARY)? (2 * BBP[i    ][j    ][0] - PC0) : (PE1);
                PW1 = (FW1 == BOUNDARY)? (2 * BBP[i - 1][j    ][0] - PC0) : (PW1);
                PN1 = (FN1 == BOUNDARY)? (2 * BBP[i    ][j    ][1] - PC0) : (PN1);
                PS1 = (FS1 == BOUNDARY)? (2 * BBP[i    ][j - 1][1] - PC0) : (PS1);

                AC0 = - 2.0 * (C1 + C2);
                AE1 = C1 - 0.5 * C7;
                AW1 = C1 + 0.5 * C7;
                AN1 = C2 - 0.5 * C8;
                AS1 = C2 + 0.5 * C8;
                RP  = (RHS - AE1 * PE1 - AW1 * PW1 - AN1 * PN1 - AS1 * PS1) / AC0 - PC0;
                ERR = ERR + RP * RP;

                P[i][j] = PC0 + OMEGA * RP;
                CNT += 1;
            }
            
        }
    }

    R = sqrt(ERR / CNT);
}

void jacobi(
    double   P[NX + 4][NY + 4],
    double  PD[NX + 4][NY + 4],
    int    FFP[NX + 4][NY + 4][2],
    double BBP[NX + 4][NY + 4][2],
    double PSI[NX + 4][NY + 4],
    double   C[NX + 4][NY + 4][4],
    int      F[NX + 4][NY + 4],
    double   OMEGA,
    double   DT,
    double  &R
) {
    double ERR = 0;
    int    CNT = 0;

    #pragma acc kernels loop independent collapse(2) present(P, PD)
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            PD[i][j] = P[i][j];
        }
    }

    #pragma acc kernels loop independent collapse(2) reduction(+:ERR, CNT) present(P, PD, FFP, BBP, PSI, C, F) copyin(OMEGA, DT) copy(ERR, CNT)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (F[i][j] == FLUID) {
                double AC0, AE1, AW1, AN1, AS1;
                double PC0, PE1, PW1, PN1, PS1;
                int    FE1, FW1, FN1, FS1;
                double C1, C2, C7, C8;
                double RHS, RP;

                RHS = PSI[i    ][j    ] / DT;
                C1  =   C[i    ][j    ][0];
                C2  =   C[i    ][j    ][1];
                C7  =   C[i    ][j    ][2];
                C8  =   C[i    ][j    ][3];
                PC0 =  PD[i    ][j    ];
                PE1 =  PD[i + 1][j    ];
                PW1 =  PD[i - 1][j    ];
                PN1 =  PD[i    ][j + 1];
                PS1 =  PD[i    ][j - 1];
                FE1 = FFP[i    ][j    ][0];
                FW1 = FFP[i - 1][j    ][0];
                FN1 = FFP[i    ][j    ][1];
                FS1 = FFP[i    ][j - 1][1];
                PE1 = (FE1 == BOUNDARY)? (2 * BBP[i    ][j    ][0] - PC0) : (PE1);
                PW1 = (FW1 == BOUNDARY)? (2 * BBP[i - 1][j    ][0] - PC0) : (PW1);
                PN1 = (FN1 == BOUNDARY)? (2 * BBP[i    ][j    ][1] - PC0) : (PN1);
                PS1 = (FS1 == BOUNDARY)? (2 * BBP[i    ][j - 1][1] - PC0) : (PS1);
     
                AC0 = - 2.0 * (C1 + C2);
                AE1 = C1 - 0.5 * C7;
                AW1 = C1 + 0.5 * C7;
                AN1 = C2 - 0.5 * C8;
                AS1 = C2 + 0.5 * C8;
                RP  = (RHS - AE1 * PE1 - AW1 * PW1 - AN1 * PN1 - AS1 * PS1) / AC0 - PC0;
                ERR = ERR + RP * RP;

                P[i][j] = PC0 + OMEGA * RP;
                CNT += 1;
            }
            
        }
    }

    R = sqrt(ERR / CNT);
}
