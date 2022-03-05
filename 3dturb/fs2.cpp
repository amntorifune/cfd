#include "topo.h"

void fs2 (
    double  U[NX + 2][NY + 2][NZ + 2][3],
    double UU[NX + 2][NY + 2][NZ + 2][3],
    double  P[NX + 2][NY + 2][NZ + 2],
    double KX[NX + 2][NY + 2][NZ + 2][3],
    double  G[NX + 2][NY + 2][NZ + 2][3],
    double  DT
) {
    #pragma acc kernels loop independent collapse(2) present(U, UU, P, KX, G) copyin(DT)
    for (int i = 2; i <= NX - 1; i ++) {
        for (int j = 2; j <= NY - 1; j ++) {
            for (int k = 2; k <= NZ - 1; k ++) {
                double D1, D2, D3;
                double PC0, PE1, PW1, PN1, PS1, PT1, PB1;
                double K1X1, K2X2, K3X3;
                double G11C0, G11E1, G22C0, G22N1, G33C0, G33T1;
                
                PC0   =  P[i    ][j    ][k    ];
                PE1   =  P[i + 1][j    ][k    ];
                PW1   =  P[i - 1][j    ][k    ];
                PN1   =  P[i    ][j + 1][k    ];
                PS1   =  P[i    ][j - 1][k    ];
                PT1   =  P[i    ][j    ][k + 1];
                PB1   =  P[i    ][j    ][k - 1];
                K1X1  = KX[i    ][j    ][k    ][0];
                K2X2  = KX[i    ][j    ][k    ][1];
                K3X3  = KX[i    ][j    ][k    ][2];
                G11C0 =  G[i    ][j    ][k    ][0];
                G22C0 =  G[i    ][j    ][k    ][1];
                G33C0 =  G[i    ][j    ][k    ][2];
                G11E1 =  G[i + 1][j    ][k    ][0];
                G22N1 =  G[i    ][j + 1][k    ][1];
                G33T1 =  G[i    ][j    ][k + 1][2];

                D1 = K1X1 * 0.5 * (PE1 - PW1);
                D2 = K2X2 * 0.5 * (PN1 - PS1);
                D3 = K3X3 * 0.5 * (PT1 - PB1);

                U[i][j][k][0] = U[i][j][k][0] - DT * D1;
                U[i][j][k][1] = U[i][j][k][1] - DT * D2;
                U[i][j][k][2] = U[i][j][k][2] - DT * D3;

                D1 = 0.5 * (G11C0 + G11E1) * (PE1 - PC0);
                D2 = 0.5 * (G22C0 + G22N1) * (PN1 - PC0);
                D3 = 0.5 * (G33C0 + G33T1) * (PT1 - PC0);

                UU[i][j][k][0] = UU[i][j][k][0] - DT * D1;
                UU[i][j][k][1] = UU[i][j][k][1] - DT * D2;
                UU[i][j][k][2] = UU[i][j][k][2] - DT * D3;
            }
        }
    }
}
