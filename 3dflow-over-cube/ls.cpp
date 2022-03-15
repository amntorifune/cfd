#include "topo.h"
#include "flag.h"
#include <math.h>

void sor(
    double   P[NX + 4][NY + 4][NZ + 4],
    int    FFP[NX + 4][NY + 4][NZ + 4][3],
    double BBP[NX + 4][NY + 4][NZ + 4][3],
    double DIV[NX + 4][NY + 4][NZ + 4],
    double   C[NX + 4][NY + 4][NZ + 4][6],
    int      F[NX + 4][NY + 4][NZ + 4],
    double   OMEGA,
    double   DT,
    double  &R
) {
    double ERR;
    int    CNT;

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            for (int k = K1; k <= K2; k ++) {
                if (!(((i + j + k) % 2 == 0) && (F[i][j][k] & FLUID))) {
                    continue;
                }
                double AC0, AE1, AW1, AN1, AS1, AT1, AB1;
                double PC0, PE1, PW1, PN1, PS1, PT1, PB1;
                double C1, C2, C3, C7, C8, C9;
                double RHS, RP;
                int    FE1, FW1, FN1, FS1, FT1, FB1;

                RHS = DIV[i    ][j    ][k    ] / DT;
                C1  =   C[i    ][j    ][k    ][0];
                C2  =   C[i    ][j    ][k    ][1];
                C3  =   C[i    ][j    ][k    ][2];
                C7  =   C[i    ][j    ][k    ][3];
                C8  =   C[i    ][j    ][k    ][4];
                C9  =   C[i    ][j    ][k    ][5];
                PC0 =   P[i    ][j    ][k    ];
                PE1 =   P[i + 1][j    ][k    ];
                PW1 =   P[i - 1][j    ][k    ];
                PN1 =   P[i    ][j + 1][k    ];
                PS1 =   P[i    ][j - 1][k    ];
                PT1 =   P[i    ][j    ][k + 1];
                PB1 =   P[i    ][j    ][k - 1];
                FE1 = FFP[i    ][j    ][k    ][x];
                FW1 = FFP[i - 1][j    ][k    ][x];
                FN1 = FFP[i    ][j    ][k    ][y];
                FS1 = FFP[i    ][j - 1][k    ][y];
                FT1 = FFP[i    ][j    ][k    ][z];
                FB1 = FFP[i    ][j    ][k - 1][z];
                PE1 = (FE1)? (2 * BBP[i    ][j    ][k    ][x] - PC0) : (PE1);
                PW1 = (FW1)? (2 * BBP[i - 1][j    ][k    ][x] - PC0) : (PW1);
                PN1 = (FN1)? (2 * BBP[i    ][j    ][k    ][y] - PC0) : (PN1);
                PS1 = (FS1)? (2 * BBP[i    ][j - 1][k    ][y] - PC0) : (PS1);
                PT1 = (FT1)? (2 * BBP[i    ][j    ][k    ][z] - PC0) : (PT1);
                PB1 = (FB1)? (2 * BBP[i    ][j    ][k - 1][z] - PC0) : (PB1);

                AC0 = - 2.0 * (C1 + C2 + C3);
                AE1 = C1 + 0.5 * C7;
                AW1 = C1 - 0.5 * C7;
                AN1 = C2 + 0.5 * C8;
                AS1 = C2 - 0.5 * C8;
                AT1 = C3 + 0.5 * C9;
                AB1 = C3 - 0.5 * C9;
                RP  = RHS - AE1 * PE1 - AW1 * PW1 - AN1 * PN1 - AS1 * PS1 - AT1 * PT1 - AB1 * PB1;
                RP  = RP / AC0;
                RP  = RP - PC0;
                ERR = ERR + RP * RP;

                P[i][j][k] = PC0 + OMEGA * RP;
                CNT += 1;
            }
        }
    }

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            for (int k = K1; k <= K2; k ++) {
                if (!(((i + j + k) % 2 == 1) && (F[i][j][k] & FLUID))) {
                    continue;
                }
                double AC0, AE1, AW1, AN1, AS1, AT1, AB1;
                double PC0, PE1, PW1, PN1, PS1, PT1, PB1;
                double C1, C2, C3, C7, C8, C9;
                double RHS, RP;
                int    FE1, FW1, FN1, FS1, FT1, FB1;

                RHS = DIV[i    ][j    ][k    ] / DT;
                C1  =   C[i    ][j    ][k    ][0];
                C2  =   C[i    ][j    ][k    ][1];
                C3  =   C[i    ][j    ][k    ][2];
                C7  =   C[i    ][j    ][k    ][3];
                C8  =   C[i    ][j    ][k    ][4];
                C9  =   C[i    ][j    ][k    ][5];
                PC0 =   P[i    ][j    ][k    ];
                PE1 =   P[i + 1][j    ][k    ];
                PW1 =   P[i - 1][j    ][k    ];
                PN1 =   P[i    ][j + 1][k    ];
                PS1 =   P[i    ][j - 1][k    ];
                PT1 =   P[i    ][j    ][k + 1];
                PB1 =   P[i    ][j    ][k - 1];
                FE1 = FFP[i    ][j    ][k    ][x];
                FW1 = FFP[i - 1][j    ][k    ][x];
                FN1 = FFP[i    ][j    ][k    ][y];
                FS1 = FFP[i    ][j - 1][k    ][y];
                FT1 = FFP[i    ][j    ][k    ][z];
                FB1 = FFP[i    ][j    ][k - 1][z];
                PE1 = (FE1)? (2 * BBP[i    ][j    ][k    ][x] - PC0) : (PE1);
                PW1 = (FW1)? (2 * BBP[i - 1][j    ][k    ][x] - PC0) : (PW1);
                PN1 = (FN1)? (2 * BBP[i    ][j    ][k    ][y] - PC0) : (PN1);
                PS1 = (FS1)? (2 * BBP[i    ][j - 1][k    ][y] - PC0) : (PS1);
                PT1 = (FT1)? (2 * BBP[i    ][j    ][k    ][z] - PC0) : (PT1);
                PB1 = (FB1)? (2 * BBP[i    ][j    ][k - 1][z] - PC0) : (PB1);

                AC0 = - 2.0 * (C1 + C2 + C3);
                AE1 = C1 + 0.5 * C7;
                AW1 = C1 - 0.5 * C7;
                AN1 = C2 + 0.5 * C8;
                AS1 = C2 - 0.5 * C8;
                AT1 = C3 + 0.5 * C9;
                AB1 = C3 - 0.5 * C9;
                RP  = RHS - AE1 * PE1 - AW1 * PW1 - AN1 * PN1 - AS1 * PS1 - AT1 * PT1 - AB1 * PB1;
                RP  = RP / AC0;
                RP  = RP - PC0;
                ERR = ERR + RP * RP;

                P[i][j][k] = PC0 + OMEGA * RP;
                CNT += 1;
            }
        }
    }

    R = sqrt(ERR / CNT);
}