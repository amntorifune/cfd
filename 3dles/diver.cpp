#include "topo.h"
#include <math.h>

void diver(
    double  JU[N1 + 2][N2 + 2][N3 + 2][3],
    double   J[N1 + 2][N2 + 2][N3 + 2],
    double DIV[N1 + 2][N2 + 2][N3 + 2],
    double  &AD
) {
    double ADL;
    int    CNT;
    ADL = 0;
    CNT = 0;
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            for (int k = K1; k <= K2; k ++) {
                double JUE, JUW, JUN, JUS, JUT, JUB;
                double J0J;
                double D0D;

                J0J =  J[i    ][j    ][k    ];
                JUE = JU[i    ][j    ][k    ][u];
                JUN = JU[i    ][j    ][k    ][v];
                JUT = JU[i    ][j    ][k    ][w];
                JUW = JU[i - 1][j    ][k    ][u];
                JUS = JU[i    ][j - 1][k    ][v];
                JUB = JU[i    ][j    ][k - 1][w];

                D0D = (JUE - JUW + JUN - JUS + JUT - JUB) / (J0J);
                ADL = ADL + D0D * D0D;
                DIV[i][j][k] = D0D;
                CNT += 1;
            }
        }
    }
    ADL = sqrt(ADL / CNT);
}