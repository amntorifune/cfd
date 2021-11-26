#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <algorithm>

using namespace std;

#define u_ 0
#define v_ 1
#define x_ 0
#define y_ 1

const int NCX      = 256;
const int NCY      = 128;
const int GHST     = 2;
const int NNX      = NCX + 2 * GHST;
const int NNY      = NCY + 2 * GHST;
const int WBOUND   = GHST;
const int SBOUND   = GHST;
const int EBOUND   = NCX + GHST;
const int NBOUND   = NCY + GHST;
const int NGX      = NCX + 1;
const int NGY      = NCY + 1;

// other constants
const double dx    = 2.0 / NCX;
const double dy    = 1.0 / NCY;
const double B0    = 2.0 / (dx * dx) + 2.0 / (dy * dy);
const double dtau  = 1.0 / B0;
const double CFL   = 1.0;

// solver parameters
const int EPSILON  = 1;
const double KAPPA = 1.0 / 3.0;
const double BETA  = (3 - KAPPA) / (1 - KAPPA);

double Re;  // Reynolds number
double dt;  // timestep length
double Tr;  // simulation time span

// data vectors
// note that Pi is the potential and P is the pressure
double U[NNX][NNY][2], UN[NNX][NNY][2], UF[NNX][NNY][2], UFN[NNX][NNY][2], UG[NGX][NGY][2]; // velocity
double Pi[NNX][NNY], PiN[NNX][NNY], P[NNX][NNY], PG[NGX][NGY];                              // pressure
int    FF[NNX][NNY];                                                                        // fluid flag
int    BF[NNX][NNY][2];                                                                     // boundary flag

// square cylinder
int cylinder[2][2] = {
    {50, 50},
    {82, 82}
};

// cell types
const int FLUID        = 0;
const int NONFLUID     = 1;

/**********************************************************************************/
/******************************* boundary parameters ******************************/
// boundary indices
const int NONE         = 0;
const int NSLIP        = 1;
const int IFLOW        = 2;
const int OFLOW        = 3;

// boundary condition types
const int DIRICHLET    = 1;
const int NEUMANN      = 2;

// dirichlet boundary values for none-slip and in-flow boundaries for velocity
const double UBCNSLIP  = 0;
const double VBCNSLIP  = 0;
const double UBCIFLOW  = 1;
const double VBCIFLOW  = 0;

// neumann boundary values for out-flow boundaries for velocity
const double DUBCOFLOW = 0;
const double DVBCOFLOW = 0;

// neumann boundary values pressure
const double DPBC      = 0;

// boundary condition types for velocity
const int    BTU[]     = {
    NONE     ,              // no boundary
    DIRICHLET,              // all the solid boundaries
    DIRICHLET,              // the in-flow boundary
    NEUMANN                 // the out-flow boundary
};
// boundary condition types for pressure
const int    BTP[]     = {
    NONE   ,                // no boundary
    NEUMANN,                // all the solid boundaries
    NEUMANN,                // the in-flow boundary
    NEUMANN                 // the out-flow boundary
};
// boundary condition values for velocity
const double BCU[][2]  = {
    {NONE     , NONE     }, // no boundary
    {UBCNSLIP , VBCNSLIP }, // all the solid boundaries
    {UBCIFLOW , VBCIFLOW }, // the in-flow boundary
    {DUBCOFLOW, DVBCOFLOW}  // the out-flow boundary
};
// boundary condition values for pressure
const double BCP[]     = {
    NONE,                   // no boundary
    DPBC,                   // all the solid boundaries
    DPBC,                   // the in-flow boundary
    DPBC                    // the out-flow boundary
};
/**********************************************************************************/

void setupGeom(int ff[NNX][NNY], int bf[NNX][NNY][2]) {
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            ff[i][k] = NONFLUID;
        }
    }
    // the tunnel
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            ff[i][k] = FLUID;
        }
    }
    // the square cylinder
    for (int i = cylinder[0][x_]; i < cylinder[1][x_]; i ++) {
        for (int k = cylinder[0][y_]; k < cylinder[1][y_]; k ++) {
            ff[i][k] = NONFLUID;
        }
    }
    
    for (int i = cylinder[0][x_]; i < cylinder[1][x_]; i ++) {
        for (int k = cylinder[0][y_]; k < cylinder[1][y_]; k ++) {
            bf[i][k][x_] = NONE;
            bf[i][k][y_] = NONE;
        }
    }
    for (int k = SBOUND; k < NBOUND; k ++) {
        // the in-flow boundary
        bf[WBOUND - 1][k][x_] = IFLOW;
        // the out-flow boundary
        bf[EBOUND - 1][k][x_] = OFLOW;
    }
    for (int i = WBOUND; i < EBOUND; i ++) {
        // the wall boundaries
        bf[i][SBOUND - 1][y_] = NSLIP;
        bf[i][NBOUND - 1][y_] = NSLIP;
    }
    // boundaries of the square cylinder
    for (int i = cylinder[0][x_]; i < cylinder[1][x_]; i ++) {
        bf[i][cylinder[0][y_] - 1][y_] = NSLIP;
        bf[i][cylinder[1][y_] - 1][y_] = NSLIP;
    }
    for (int k = cylinder[0][y_]; k < cylinder[1][y_]; k ++) {
        bf[cylinder[0][x_] - 1][k][x_] = NSLIP;
        bf[cylinder[1][x_] - 1][k][x_] = NSLIP;
    }
}

void bE1(
    int i, int k, 
    int bf[NNX][NNY][2], 
    int &ee1, int &ew1, int &en1, int &es1
) {
    ee1 = bf[i    ][k    ][x_];
    en1 = bf[i    ][k    ][y_];
    ew1 = bf[i - 1][k    ][x_];
    es1 = bf[i    ][k - 1][y_];
}

void bE2(
    int i, int k,
    int bf[NNX][NNY][2], 
    int &ee1, int &ew1, int &en1, int &es1, 
    int &ee2, int &ew2, int &en2, int &es2
) {
    ee1 = bf[i    ][k    ][x_];
    en1 = bf[i    ][k    ][y_];
    ew1 = bf[i - 1][k    ][x_];
    es1 = bf[i    ][k - 1][y_];

    ee2 = bf[i + 1][k    ][x_];
    en2 = bf[i    ][k + 1][y_];
    ew2 = bf[i - 2][k    ][x_];
    es2 = bf[i    ][k - 2][y_];
}

// to extrapolate cell center values
double extrapolCC(int type, double base, double dist, double bcvalue) {
    if (type == DIRICHLET) {
        return (2 * bcvalue - base);
    }
    else if (type == NEUMANN) {
        return (base + dist * bcvalue);
    }
    return 0.0;
}

// to extrapolate cell face values
double extrapolCF(int type, double base, double dist, double bcvalue) {
    if (type == DIRICHLET) {
        return bcvalue;
    }
    else if (type == NEUMANN) {
        return (base + dist * bcvalue);
    }
    return 0.0;
}

// calculate ||div(U)||
double calcdiv(
    double uF[NNX][NNY][2], double u[NNX][NNY][2], 
    int    ff[NNX][NNY]   , int    bf[NNX][NNY][2]
) {
    double tdiver = 0;
    #pragma acc kernels loop independent reduction(+:tdiver) collapse(2) present(uF, u, ff, bf) copy(tdiver)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double ufe, ufw, vfn, vfs; // velocity at each direction
                double uc0, vc0;           // velocity at cell center
                double dudx, dvdy;         // velocity gradient at cell center
                double diver;
                int    ee1, ew1, en1, es1;     // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                ufe = (ee1 == NONE)? uF[i    ][k    ][u_] : extrapolCF(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                ufw = (ew1 == NONE)? uF[i - 1][k    ][u_] : extrapolCF(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                vfn = (en1 == NONE)? uF[i    ][k    ][v_] : extrapolCF(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                vfs = (es1 == NONE)? uF[i    ][k - 1][v_] : extrapolCF(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);

                dudx   = (ufe - ufw) / (dx);
                dvdy   = (vfn - vfs) / (dy);
                diver  = dudx + dvdy;
                tdiver = tdiver + diver * diver;
            }
        }
    }
    return sqrt(tdiver);
}

// interpolate cell face values from cell center values
double calcUF(
    double u[NNX][NNY][2], double uF[NNX][NNY][2], 
    int    ff[NNX][NNY]  , int    bf[NNX][NNY][2]
) {
    #pragma acc kernels loop independent collapse(2) present(u, uF, ff, bf)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double ue1, uc0, uw1, vn1, vc0, vs1; // velocity at each direction
                double ufe, ufw, vfn, vfs;           // velocity at cell faces
                int    ee1, ew1, en1, es1;           // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0 = u[i][k][u_];
                vc0 = u[i][k][v_];

                ue1 = (ee1 == NONE)? u[i + 1][k    ][u_] : extrapolCC(BTU[ee1], uc0,   dx, BCU[ee1][u_]);
                uw1 = (ew1 == NONE)? u[i - 1][k    ][u_] : extrapolCC(BTU[ew1], uc0, - dx, BCU[ew1][u_]);
                vn1 = (en1 == NONE)? u[i    ][k + 1][v_] : extrapolCC(BTU[en1], vc0,   dy, BCU[en1][v_]);
                vs1 = (es1 == NONE)? u[i    ][k - 1][v_] : extrapolCC(BTU[es1], vc0, - dy, BCU[es1][v_]);

                ufe = (ee1 == NONE)? 0.5 * (uc0 + ue1) : extrapolCF(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                ufw = (ew1 == NONE)? 0.5 * (uw1 + uc0) : extrapolCF(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                vfn = (en1 == NONE)? 0.5 * (vc0 + vn1) : extrapolCF(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                vfs = (es1 == NONE)? 0.5 * (vs1 + vc0) : extrapolCF(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);

                uF[i    ][k    ][u_] = ufe;
                uF[i - 1][k    ][u_] = ufw;
                uF[i    ][k    ][v_] = vfn;
                uF[i    ][k - 1][v_] = vfs;
            }
        }
    }
}

// initialize pressure potential π
void initPi(double pi[NNX][NNY]) {
    #pragma acc kernels loop collapse(2) present(pi)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            pi[i][k] = 0;
        }
    }
}

// initialize velocity and pressure
void init(
    double u[NNX][NNY][2], double uF[NNX][NNY][2], 
    double p[NNX][NNY]   , double pi[NNX][NNY]   , 
    int    ff[NNX][NNY]  , int    bf[NNX][NNY][2]
) {
    #pragma acc kernels loop collapse(2) present(u)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            u[i][k][u_] = 0;
            u[i][k][v_] = 0;
        }
    }
    calcUF(u, uF);
    initPi(pi);
    #pragma acc kernels loop collapse(2) present(p)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            p[i][k] = 0;
        }
    }
}