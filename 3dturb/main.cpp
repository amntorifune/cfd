#include "funclist.h"
#include "varlist.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    geom(X);
    init(U, UD, UC, UU, P, SGS, X, KX, J, C, G);

    for (STEP = 1; STEP <= NSTEP; STEP ++) {
        fs1(U, UD, UU, SGS, KX, J, C, ALPHA, DT, RI);
        contra(U, UC, UU, KX, J);
        // solid(U, UU);
        diver(UU, J, DIV, AD);
        while (AD > EDIV) {
            sor(P, DIV, C, OMEGA, EPOI, MAXIT, DT, IT, R);
            avp0(P);
            bcp(P, U, C, KX, RI);
            fs2(U, UU, P, KX, G, DT);
            // solid(U, UU);
            bcu(U, UU, UD, KX, DT);
            diver(UU, J, DIV, AD);
            printf("\rstep %5d: p(%3d, %5.8lf), d(%5.8lf)", STEP, IT, R, AD);
        }
    }

    fio(U, P, X, (char*)"o.csv");

    return 0;
}