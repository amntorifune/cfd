#include "topo.h"
#include <math.h>

void diver(
    double  UU[NX + 2][NY + 2][NZ + 2][3],
    double   J[NX + 2][NY + 2][NZ + 2],
    double DIV[NX + 2][NY + 2][NZ + 2],
    double  &AD
) {
    double UUE, UUW, UUN, UUS, UUT, UUB;
    double JC0;
    double D;

    AD = 0;
    for (int i = 2; i <= NX - 1; i ++) {
        for (int j = 2; j < NY - 1; j ++) {
            for (int k = 2; k < NZ - 1; k ++) {
                JC0 =  J[i    ][j    ][k    ];
                UUE = UU[i    ][j    ][k    ][0];
                UUN = UU[i    ][j    ][k    ][1];
                UUT = UU[i    ][j    ][k    ][2];
                UUW = UU[i - 1][j    ][k    ][0];
                UUS = UU[i    ][j - 1][k    ][1];
                UUB = UU[i    ][j    ][k - 1][2];

                D  = (UUE - UUW + UUN - UUS + UUT - UUB) / (JC0);
                AD = AD + D * D;
                DIV[i][j][k] = D;
            }
        }
    }
    AD = sqrt(AD / NXYZ);
}