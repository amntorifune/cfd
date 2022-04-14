#include "topo.h"
#include "flag.h"

void p0gradient(
    double   P[NX + 4][NY + 4],
    double BBP[NX + 4][NY + 4][2]
) {
//  east boundary value designated p=0
    #pragma acc kernels loop independent present(P, BBP)
    for (int j = J1; j <= J2; j ++) {
        BBP[I2][j][0] = P[I2][j];
    }

//  west boundary

    #pragma acc kernels loop independent present(P, BBP)
    for (int j = J1; j <= J2; j ++) {
        BBP[I1 - 1][j][0] = P[I1][j];
    }

//  north boundary

    #pragma acc kernels loop independent present(P, BBP)
    for (int i = I1; i <= I2; i ++) {
        BBP[i][J2][1] = P[i][J2];
    }

//  south boundary

    #pragma acc kernels loop independent present(P, BBP)
    for (int i = I1; i <= I2; i ++) {
        BBP[i][J1 - 1][1] = P[i][J1];
    }

//  square east face

    #pragma acc kernels loop independent present(P, BBP)
    for (int j = J3; j <= J4; j ++) {
        BBP[I4][j][0] = P[I4 + 1][j];
    }

//  square west face

    #pragma acc kernels loop independent present(P, BBP)
    for (int j = J3; j <= J4; j ++) {
        BBP[I3 - 1][j][0] = P[I3 - 1][j];
    }

//  square north face

    #pragma acc kernels loop independent present(P, BBP)
    for (int i = I3; i <= I4; i ++) {
        BBP[i][J4][1] = P[i][J4 + 1];
    }

//  square south face

    #pragma acc kernels loop independent present(P, BBP)
    for (int i = I3; i <= I4; i ++) {
        BBP[i][J3 - 1][1] = P[i][J3 - 1];
    }
}
