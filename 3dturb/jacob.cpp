#include "topo.h"
#include <math.h>

void jacob(
    double   P[NX + 2][NY + 2][NZ + 2],
    double DIV[NX + 2][NY + 2][NZ + 2],
    double   C[NX + 2][NY + 2][NZ + 2][6],
    double   OMEGA,
    double   E,
    int      MAXIT,
    double   DT,
    int     &IT,
    double  &R
) {
    double ERR;

    IT = 0;
    R  = 1E9;
    while (IT < MAXIT && R > E) {
        ERR = 0.0;
        #pragma acc kernels loop independent reduction(+:ERR) collapse(2) present(C, P, DIV) copyin(OMEGA, DT) copy(ERR)
        for (int i = 2; i <= NX - 1; i ++) {
            for (int j = 2; j <= NY - 1; j ++) {
                for (int k = 2; k <= NZ - 1; k ++) {
                    double AC0, AE1, AW1, AN1, AS1, AT1, AB1;
                    double PC0, PE1, PW1, PN1, PS1, PT1, PB1;
                    double C1, C2, C3, C7, C8, C9;
                    double RHS, RP;

                    RHS = DIV[i][j][k] / DT;
                    C1  = C[i][j][k][0];
                    C2  = C[i][j][k][1];
                    C3  = C[i][j][k][2];
                    C7  = C[i][j][k][3];
                    C8  = C[i][j][k][4];
                    C9  = C[i][j][k][5];
                    PC0 = P[i    ][j    ][k    ];
                    PE1 = P[i + 1][j    ][k    ];
                    PW1 = P[i - 1][j    ][k    ];
                    PN1 = P[i    ][j + 1][k    ];
                    PS1 = P[i    ][j - 1][k    ];
                    PT1 = P[i    ][j    ][k + 1];
                    PB1 = P[i    ][j    ][k - 1];
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
                }
            }
        }
        R = sqrt(ERR / NXYZ);
        IT += 1;
    }
}
