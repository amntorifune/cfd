#include "topo.h"
#include "flag.h"
#include <math.h>

void diver(
    double  UU[NX + 4][NY + 4][2],
    double   J[NX + 4][NY + 4],
    double DIV[NX + 4][NY + 4],
    int      F[NX + 4][NY + 4],
    double  &AD
) {
    double ADL = 0;
    int    CNT = 0;

    #pragma acc kernels loop independent collapse(2) reduction(+:ADL, CNT) present(UU, J, DIV, F) copy(ADL, CNT)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (F[i][j] == FLUID) {
                double UUE, UUW, UUN, UUS;
                double J0J;
                double D;

                J0J =  J[i    ][j    ];
                UUE = UU[i    ][j    ][0];
                UUN = UU[i    ][j    ][1];
                UUW = UU[i - 1][j    ][0];
                UUS = UU[i    ][j - 1][1];

                D   = (UUE - UUW + UUN - UUS) / (J0J);
                ADL = ADL + D * D;
                DIV[i][j] = D;
                CNT += 1;
            }
        }
    }
    ADL = sqrt(ADL / CNT);
    AD  = ADL;
}
