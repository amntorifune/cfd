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

const int NCX      = 128;
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
const double dx    = 1.0 / NCX;
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

// cell types
const int FLUID        = 1;
const int NONFLUID     = 0;

/**********************************************************************************/
/******************************* boundary parameters ******************************/
// boundary indices
// 2 boundaries: the walls and the lid
const int NONE         = 0;
const int WALL         = 1;
const int LID          = 2;

// boundary condition types
const int DIRICHLET    = 1;
const int NEUMANN      = 2;

// dirichlet boundary values for verlocity
const double UBCWALL   = 0.0;
const double VBCWALL   = 0.0;
const double UBCLID    = 1.0;
const double VBCLID    = 0.0;

// neumann boundary values pressure
const double DPBC      = 0;

// boundary condition type for velocity of each boundary
const int    BTU[]     = {
    NONE     ,          // non-boundary
    DIRICHLET,          // the walls
    DIRICHLET,          // the lid
};
// boundary condition type for pressure of each boundary
const int    BTP[]     = {
    NONE   ,            // non-boundary
    NEUMANN,            // the walls
    NEUMANN,            // the lid
};
// boundary condition value for velocity of each boundary
const double BCU[][2]  = {
    {NONE   , NONE   }, // non-boundary
    {UBCWALL, VBCWALL}, // the walls
    {UBCLID , VBCLID }, // the lid
};
// boundary condition value for pressure of each boundary
const double BCP[]     = {
    NONE,               // non-boundary
    DPBC,               // the walls
    DPBC                // the lid
};
/**********************************************************************************/
/**********************************************************************************/
void setupGeom(int ff[NNX][NNY], int bf[NNX][NNY][2]) {
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            ff[i][k] = NONFLUID;
        }
    }
    // inside the cavity
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            ff[i][k] = FLUID;
        }
    }
    
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            bf[i][k][x_] = NONE;
            bf[i][k][y_] = NONE;
        }
    }
    for (int k = SBOUND; k < NBOUND; k ++) {
        // the left wall
        bf[WBOUND - 1][k][x_] = WALL;
        // the right wall
        bf[EBOUND - 1][k][x_] = WALL;
    }
    for (int i = WBOUND; i < EBOUND; i ++) {
        // the bottom wall
        bf[i][SBOUND - 1][y_] = WALL;
        // the lid
        bf[i][NBOUND - 1][y_] = LID;
    }
}
/**********************************************************************************/

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
double bcc(int type, double base, double dist, double bcvalue) {
    if (type == DIRICHLET) {
        return (2 * bcvalue - base);
    }
    else if (type == NEUMANN) {
        return (base + dist * bcvalue);
    }
    return 0.0;
}

// to extrapolate cell face values
double bcf(int type, double base, double dist, double bcvalue) {
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
    #pragma acc kernels loop independent reduction(+:tdiver) collapse(2) present(uF, u, ff, bf, BTU, BCU) copy(tdiver)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double ufe, ufw, vfn, vfs; // velocity at each direction
                double uc0, vc0;           // velocity at cell center
                double dudx, dvdy;         // velocity gradient at cell center
                double diver;
                int    ee1, ew1, en1, es1; // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0 =  u[i    ][k    ][u_];
                vc0 =  u[i    ][k    ][v_];

                ufe = uF[i    ][k    ][u_];
                ufw = uF[i - 1][k    ][u_];
                vfn = uF[i    ][k    ][v_];
                vfs = uF[i    ][k - 1][v_];

                ufe = (ee1 == NONE)? ufe : bcf(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                ufw = (ew1 == NONE)? ufw : bcf(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                vfn = (en1 == NONE)? vfn : bcf(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                vfs = (es1 == NONE)? vfs : bcf(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);

                dudx   = (ufe - ufw) / dx;
                dvdy   = (vfn - vfs) / dy;
                diver  = dudx + dvdy;
                tdiver = tdiver + diver * diver;
            }
        }
    }
    return sqrt(tdiver);
}

// interpolate cell face values from cell center values
void calcUF(
    double u[NNX][NNY][2], double uF[NNX][NNY][2], 
    int    ff[NNX][NNY]  , int    bf[NNX][NNY][2]
) {
    #pragma acc kernels loop independent collapse(2) present(u, uF, ff, bf, BTU, BCU)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double ue1, uc0, uw1, vn1, vc0, vs1; // velocity at each direction
                double ufe, ufw, vfn, vfs;           // velocity at cell faces
                int    ee1, ew1, en1, es1;           // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0 = u[i    ][k    ][u_];
                vc0 = u[i    ][k    ][v_];

                ue1 = u[i + 1][k    ][u_];
                uw1 = u[i - 1][k    ][u_];
                vn1 = u[i    ][k + 1][v_];
                vs1 = u[i    ][k - 1][v_];

                ue1 = (ee1 == NONE)? ue1 : bcc(BTU[ee1], uc0,   dx, BCU[ee1][u_]);
                uw1 = (ew1 == NONE)? uw1 : bcc(BTU[ew1], uc0, - dx, BCU[ew1][u_]);
                vn1 = (en1 == NONE)? vn1 : bcc(BTU[en1], vc0,   dy, BCU[en1][v_]);
                vs1 = (es1 == NONE)? vs1 : bcc(BTU[es1], vc0, - dy, BCU[es1][v_]);

                ufe = 0.5 * (uc0 + ue1);
                ufw = 0.5 * (uw1 + uc0);
                vfn = 0.5 * (vc0 + vn1);
                vfs = 0.5 * (vs1 + vc0);

                ufe = (ee1 == NONE)? ufe : bcf(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                ufw = (ew1 == NONE)? ufw : bcf(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                vfn = (en1 == NONE)? vfn : bcf(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                vfs = (es1 == NONE)? vfs : bcf(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);

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
    calcUF(u, uF, ff, bf);
    initPi(pi);
    #pragma acc kernels loop collapse(2) present(p)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            p[i][k] = 0;
        }
    }
}

// copy U to UN or UF to UFN
void cpUt(double ut[NNX][NNY][2], double utN[NNX][NNY][2]) {
    #pragma acc kernels loop independent collapse(2) present(ut, utN)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            utN[i][k][u_] = ut[i][k][u_];
            utN[i][k][v_] = ut[i][k][v_];
        }
    }
}

// copy Pi to PiN
void cpPi(double pi[NNX][NNY], double piN[NNX][NNY]) {
    #pragma acc kernels loop independent collapse(2) present(pi, piN)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            piN[i][k] = pi[i][k];
        }
    }
}

int sgn(double x) {
    if (x == 0) {
        return 0;
    }
    else if (x > 0) {
        return 1;
    }
    else {
        return - 1;
    }
}

double minmod(double x, double y) {
    return sgn(x) * max(0.0, min(abs(x), sgn(x) * y));
}

// interpolated value φ+ at cell face: φA φB | φC φD
// φ+ = φC - ε/4 * ((1 - κ) * Δ+ + (1 + κ) * Δ-)
// Δ+ = minmod(d+, β * d-), Δ- = minmod(d-, β * d+)
// d+ = φD - φC, d- = φC - φB
double interpolP(double b, double c, double d) {
    double dP = d - c;
    double dM = c - b;
    double DP = minmod(dP, BETA * dM);
    double DM = minmod(dM, BETA * dP);
    return (c - EPSILON * 0.25 * ((1 - KAPPA) * DP + (1 + KAPPA) * DM));
}
// interpolated value φ- at cell face: φA φB | φC φD
// φ- = φB + ε/4 * ((1 - κ) * Δ- + (1 + κ) * Δ+)
// Δ+ = minmod(d+, β * d-), Δ- = minmod(d-, β * d+)
// d+ = φC - φB, d- = φB - φA
double interpolM(double a, double b, double c) {
    double dP = c - b;
    double dM = b - a;
    double DP = minmod(dP, BETA * dM);
    double DM = minmod(dM, BETA * dP);
    return (b + EPSILON * 0.25 * ((1 - KAPPA) * DM + (1 + KAPPA) * DP));
}
// numerical flux at cell face: φA φB | φC φD
// f~ = 1/2 * ((f+ + f-) - |u@face| * (φ+ - φ-))
// f+- = u@face * φ+-
double flux(double a, double b, double c, double d, double uf) {
    double phiP = interpolP(b, c, d);
    double phiM = interpolM(a, b, c);
    double flxP = uf * phiP;
    double flxM = uf * phiM;
    return 0.5 * (flxP + flxM - abs(uf) * (phiP - phiM));
}

// u* = u + Δt * (- Convection + Viscoucity - grad(P))
void prediction(
    double u[NNX][NNY][2], double uN[NNX][NNY][2], double uF[NNX][NNY][2],
    double p[NNX][NNY],
    int    ff[NNX][NNY],   int    bf[NNX][NNY][2]
) {
    cpUt(u, uN);
    #pragma acc kernels loop independent collapse(2) present(u, uN, uF, p, ff, bf, BTU, BCU, BTP, BCP)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double fue, fuw, fun, fus, fve, fvw, fvn, fvs;      // flux variables
                double uw2, uw1, us2, us1, uc0, un1, un2, ue1, ue2; // velocity at each direction
                double vw2, vw1, vs2, vs1, vc0, vn1, vn2, ve1, ve2; // velocity at each direction
                double ufe, ufw, vfn, vfs;                          // velocity at cell faces
                double duudx, dvudy, duvdx, dvvdy;                  // 1st derivatives
                double ddudxx, ddudyy, ddvdxx, ddvdyy;              // 2nd derivatives
                double pe1, pw1, pc0, pn1, ps1;                     // pressure at each direction
                double dprdxfe, dprdxfw, dprdyfn, dprdyfs;          // pressure derivatives at cell faces
                double dprdxcc, dprdycc;                            // pressure derivatives at cell center
                int    ee1, ew1, en1, es1;                          // first level boundary indicator
                int    ee2, ew2, en2, es2;                          // second level boundary indicator
                bE2(i, k, bf, ee1, ew1, en1, es1, ee2, ew2, en2, es2);

                uc0 = uN[i    ][k    ][u_];
                vc0 = uN[i    ][k    ][v_];

                ue1 = uN[i + 1][k    ][u_];
                uw1 = uN[i - 1][k    ][u_];
                un1 = uN[i    ][k + 1][u_];
                us1 = uN[i    ][k - 1][u_];
                ve1 = uN[i + 1][k    ][v_];
                vw1 = uN[i - 1][k    ][v_];
                vn1 = uN[i    ][k + 1][v_];
                vs1 = uN[i    ][k - 1][v_];

                ue2 = uN[i + 2][k    ][u_];
                uw2 = uN[i - 2][k    ][u_];
                un2 = uN[i    ][k + 2][u_];
                us2 = uN[i    ][k - 2][u_];
                ve2 = uN[i + 2][k    ][v_];
                vw2 = uN[i - 2][k    ][v_];
                vn2 = uN[i    ][k + 2][v_];
                vs2 = uN[i    ][k - 2][v_];

                ufe = uF[i    ][k    ][u_];
                ufw = uF[i - 1][k    ][u_];
                vfn = uF[i    ][k    ][v_];
                vfs = uF[i    ][k - 1][v_];

                pc0 = p[i    ][k    ];

                pe1 = p[i + 1][k    ];
                pw1 = p[i - 1][k    ];
                pn1 = p[i    ][k + 1];
                ps1 = p[i    ][k - 1];

                ue1 = (ee1 == NONE)? ue1 : bcc(BTU[ee1], uc0,   dx, BCU[ee1][u_]);
                uw1 = (ew1 == NONE)? uw1 : bcc(BTU[ew1], uc0, - dx, BCU[ew1][u_]);
                un1 = (en1 == NONE)? un1 : bcc(BTU[en1], uc0,   dy, BCU[en1][u_]);
                us1 = (es1 == NONE)? us1 : bcc(BTU[es1], uc0, - dy, BCU[es1][u_]);
                ve1 = (ee1 == NONE)? ve1 : bcc(BTU[ee1], vc0,   dx, BCU[ee1][v_]);
                vw1 = (ew1 == NONE)? vw1 : bcc(BTU[ew1], vc0, - dx, BCU[ew1][v_]);
                vn1 = (en1 == NONE)? vn1 : bcc(BTU[en1], vc0,   dy, BCU[en1][v_]);
                vs1 = (es1 == NONE)? vs1 : bcc(BTU[es1], vc0, - dy, BCU[es1][v_]);

                ue2 = (ee2 == NONE)? ue2 : bcc(BTU[ee2], ue1,   dx, BCU[ee2][u_]);
                uw2 = (ew2 == NONE)? uw2 : bcc(BTU[ew2], uw1, - dx, BCU[ew2][u_]);
                un2 = (en2 == NONE)? un2 : bcc(BTU[en2], un1,   dy, BCU[en2][u_]);
                us2 = (es2 == NONE)? us2 : bcc(BTU[es2], us1, - dy, BCU[es2][u_]);
                ve2 = (ee2 == NONE)? ve2 : bcc(BTU[ee2], ve1,   dx, BCU[ee2][v_]);
                vw2 = (ew2 == NONE)? vw2 : bcc(BTU[ew2], vw1, - dx, BCU[ew2][v_]);
                vn2 = (en2 == NONE)? vn2 : bcc(BTU[en2], vn1,   dy, BCU[en2][v_]);
                vs2 = (es2 == NONE)? vs2 : bcc(BTU[es2], vs1, - dy, BCU[es2][v_]);

                ue2 = (ee1 == NONE)? ue2 : bcc(BTU[ee1], uw1,   3 * dx, BCU[ee1][u_]);
                uw2 = (ew1 == NONE)? uw2 : bcc(BTU[ew1], ue1, - 3 * dx, BCU[ew1][u_]);
                un2 = (en1 == NONE)? un2 : bcc(BTU[en1], us1,   3 * dy, BCU[en1][u_]);
                us2 = (es1 == NONE)? us2 : bcc(BTU[es1], un1, - 3 * dy, BCU[es1][u_]);
                ve2 = (ee1 == NONE)? ve2 : bcc(BTU[ee1], vw1,   3 * dx, BCU[ee1][v_]);
                vw2 = (ew1 == NONE)? vw2 : bcc(BTU[ew1], ve1, - 3 * dx, BCU[ew1][v_]);
                vn2 = (en1 == NONE)? vn2 : bcc(BTU[en1], vs1,   3 * dy, BCU[en1][v_]);
                vs2 = (es1 == NONE)? vs2 : bcc(BTU[es1], vn1, - 3 * dy, BCU[es1][v_]);

                ufe = (ee1 == NONE)? ufe : bcf(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                ufw = (ew1 == NONE)? ufw : bcf(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                vfn = (en1 == NONE)? vfn : bcf(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                vfs = (es1 == NONE)? vfs : bcf(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);

                
                pe1 = (ee1 == NONE)? pe1 : bcc(BTP[ee1], pc0,   dx, BCP[ee1]);
                pw1 = (ew1 == NONE)? pw1 : bcc(BTP[ew1], pc0, - dx, BCP[ew1]);
                pn1 = (en1 == NONE)? pn1 : bcc(BTP[en1], pc0,   dy, BCP[en1]);
                ps1 = (es1 == NONE)? ps1 : bcc(BTP[es1], pc0, - dy, BCP[es1]);

                fue = flux(uw1, uc0, ue1, ue2, ufe);
                fuw = flux(uw2, uw1, uc0, ue1, ufw);
                fun = flux(us1, uc0, un1, un2, vfn);
                fus = flux(us2, us1, uc0, un1, vfs);
                fve = flux(vw1, vc0, ve1, ve2, ufe);
                fvw = flux(vw2, vw1, vc0, ve1, ufw);
                fvn = flux(vs1, vc0, vn1, vn2, vfn);
                fvs = flux(vs2, vs1, vc0, vn1, vfs);

                duudx   = (fue - fuw) / dx;
                dvudy   = (fun - fus) / dy;
                duvdx   = (fve - fvw) / dx;
                dvvdy   = (fvn - fvs) / dy;
                ddudxx  = (uw1 - 2 * uc0 + ue1) / (dx * dx);
                ddudyy  = (us1 - 2 * uc0 + un1) / (dy * dy);
                ddvdxx  = (vw1 - 2 * vc0 + ve1) / (dx * dx);
                ddvdyy  = (vs1 - 2 * vc0 + vn1) / (dy * dy);

                dprdxfe = (pe1 - pc0) / dx;
                dprdxfw = (pc0 - pw1) / dx;
                dprdyfn = (pn1 - pc0) / dy;
                dprdyfs = (pc0 - ps1) / dy;
                dprdxcc = (dprdxfe + dprdxfw) * 0.5;
                dprdycc = (dprdyfn + dprdyfs) * 0.5;

                uc0     = uc0 + dt * (- duudx - dvudy - dprdxcc + (ddudxx + ddudyy) / Re);
                vc0     = vc0 + dt * (- duvdx - dvvdy - dprdycc + (ddvdxx + ddvdyy) / Re);

                u[i][k][u_] = uc0;
                u[i][k][v_] = vc0;
            }
        }
    }
    calcUF(u, uF, ff, bf);
}

// div(grad(π)) = div(u*) / Δt
// π is the pressure potential, which is different from pressure P
int poisson(
    double pi[NNX][NNY]   , double piN[NNX][NNY]  , 
    double uF[NNX][NNY][2], double u[NNX][NNY][2] , 
    int    ff[NNX][NNY]   , int    bf[NNX][NNY][2]
) {
    double E     = 1E-5;
    int    MAXIT = 100;
    bool   con   = false;
    int    it    = 0;
    double R     = 0;
    while (con == false) {
        R   = 0;
        con = true;
        cpPi(pi, piN);
        #pragma acc kernels loop independent reduction(+:R) collapse(2) present(pi, piN, uF, u, ff, bf, BTU, BCU, BTP, BCP) copy(R)
        for (int i = WBOUND; i < EBOUND; i ++) {
            for (int k = SBOUND; k < NBOUND; k ++) {
                if (ff[i][k] == FLUID) {
                    double ufe, ufw, vfn, vfs;             // velocity at cell faces
                    double uc0, vc0;                       // velocity at cell center
                    double pie1, piw1, pic0, pin1, pis1;   // pressure potential at each direction
                    double dudx, dvdy;                     // velocity gradient at cell center
                    double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // 1st derivatives at cell faces
                    double ddpdxx, ddpdyy;                 // 2nd derivatives at cell center
                    double psi, res;
                    int    ee1, ew1, en1, es1;             // first level boundary indicator
                    bE1(i, k, bf, ee1, ew1, en1, es1);

                    uc0  =  u[i    ][k    ][u_];
                    vc0  =  u[i    ][k    ][v_];

                    ufe  = uF[i    ][k    ][u_];
                    ufw  = uF[i - 1][k    ][u_];
                    vfn  = uF[i    ][k    ][v_];
                    vfs  = uF[i    ][k - 1][v_];

                    pic0 = pi[i    ][k    ];

                    pie1 = pi[i + 1][k    ];
                    piw1 = pi[i - 1][k    ];
                    pin1 = pi[i    ][k + 1];
                    pis1 = pi[i    ][k - 1];

                    ufe  = (ee1 == NONE)? ufe : bcf(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                    ufw  = (ew1 == NONE)? ufw : bcf(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                    vfn  = (en1 == NONE)? vfn : bcf(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                    vfs  = (es1 == NONE)? vfs : bcf(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);

                    pie1 = (ee1 == NONE)? pie1 : bcc(BTP[ee1], pic0,   dx, BCP[ee1]);
                    piw1 = (ew1 == NONE)? piw1 : bcc(BTP[ew1], pic0, - dx, BCP[ew1]);
                    pin1 = (en1 == NONE)? pin1 : bcc(BTP[en1], pic0,   dy, BCP[en1]);
                    pis1 = (es1 == NONE)? pis1 : bcc(BTP[es1], pic0, - dy, BCP[es1]);

                    dudx     = (ufe - ufw) / dx;
                    dvdy     = (vfn - vfs) / dy;

                    dpdxfe   = (pie1 - pic0) / dx;
                    dpdxfw   = (pic0 - piw1) / dx;
                    dpdyfn   = (pin1 - pic0) / dy;
                    dpdyfs   = (pic0 - pis1) / dy;
                    ddpdxx   = (dpdxfe - dpdxfw) / dx;
                    ddpdyy   = (dpdyfn - dpdyfs) / dy;

                    psi      = (dudx + dvdy) / dt;
                    res      = dtau * (ddpdxx + ddpdyy - psi);

                    pic0     = pic0 + res;
                    R        = R + res * res;
                    pi[i][k] = pic0;
                }
            }
        }
        if (sqrt(R) > E) {
            con = false;
        }
        it += 1;
        if (it >= MAXIT) {
            con = true;
        }
    }
    return it;
}

// P = P + π
void updateP(double p[NNX][NNY], double pi[NNX][NNY], int ff[NNX][NNY]) {
    #pragma acc kernels loop independent collapse(2) present(p, pi, ff)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                p[i][k] = p[i][k] + pi[i][k];
            }
        }
    }
}

// u = u* - Δt * grad(π)
double projection(
    double u[NNX][NNY][2], double uF[NNX][NNY][2], double uFN[NNX][NNY][2],
    double pi[NNX][NNY]  , double p[NNX][NNY], 
    int    ff[NNX][NNY]  , int    bf[NNX][NNY][2]
) {
    cpUt(uF, uFN);
    #pragma acc kernels loop independent collapse(2) present(u, uF, uFN, pi, ff, bf, BTU, BCU, BTP, BCP)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double pie1, piw1, pic0, pin1, pis1;   // pressure potential at each direction
                double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // pressure gradient at cell faces
                double dpdxcc, dpdycc;                 // pressure gradient at cell center
                double uc0, vc0;                       // velocity at cell center
                double ufe, ufw, vfn, vfs;             // velocity at cell faces
                int    ee1, ew1, en1, es1;             // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0  =   u[i    ][k    ][u_];
                vc0  =   u[i    ][k    ][v_];

                ufe  = uFN[i    ][k    ][u_];
                ufw  = uFN[i - 1][k    ][u_];
                vfn  = uFN[i    ][k    ][v_];
                vfs  = uFN[i    ][k - 1][v_];

                pic0 =  pi[i    ][k    ];

                pie1 =  pi[i + 1][k    ];
                piw1 =  pi[i - 1][k    ];
                pin1 =  pi[i    ][k + 1];
                pis1 =  pi[i    ][k - 1];

                ufe  = (ee1 == NONE)? ufe : bcf(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                ufw  = (ew1 == NONE)? ufw : bcf(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                vfn  = (en1 == NONE)? vfn : bcf(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                vfs  = (es1 == NONE)? vfs : bcf(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);
                
                pie1 = (ee1 == NONE)? pie1 : bcc(BTP[ee1], pic0,   dx, BCP[ee1]);
                piw1 = (ew1 == NONE)? piw1 : bcc(BTP[ew1], pic0, - dx, BCP[ew1]);
                pin1 = (en1 == NONE)? pin1 : bcc(BTP[en1], pic0,   dy, BCP[en1]);
                pis1 = (es1 == NONE)? pis1 : bcc(BTP[es1], pic0, - dy, BCP[es1]);

                dpdxfe = (pie1 - pic0) / dx;
                dpdxfw = (pic0 - piw1) / dx;
                dpdyfn = (pin1 - pic0) / dy;
                dpdyfs = (pic0 - pis1) / dy;
                dpdxcc = 0.5 * (dpdxfe + dpdxfw);
                dpdycc = 0.5 * (dpdyfn + dpdyfs);

                uc0 = uc0 - dt * dpdxcc;
                vc0 = vc0 - dt * dpdycc;
                
                ufe = ufe - dt * dpdxfe;
                ufw = ufw - dt * dpdxfw;
                vfn = vfn - dt * dpdyfn;
                vfs = vfs - dt * dpdyfs;

                ufe = (ee1 == NONE)? ufe : bcf(BTU[ee1], uc0,   0.5 * dx, BCU[ee1][u_]);
                ufw = (ew1 == NONE)? ufw : bcf(BTU[ew1], uc0, - 0.5 * dx, BCU[ew1][u_]);
                vfn = (en1 == NONE)? vfn : bcf(BTU[en1], vc0,   0.5 * dy, BCU[en1][v_]);
                vfs = (es1 == NONE)? vfs : bcf(BTU[es1], vc0, - 0.5 * dy, BCU[es1][v_]);
                
                uF[i    ][k    ][u_] = ufe;
                uF[i - 1][k    ][u_] = ufw;
                uF[i    ][k    ][v_] = vfn;
                uF[i    ][k - 1][v_] = vfs;

                u[i][k][u_]          = uc0;
                u[i][k][v_]          = vc0;
            }
        }
    }
    updateP(p, pi, ff);
    return calcdiv(uF, u, ff, bf);
}

void ns2d(
    double u[NNX][NNY][2] , double uN[NNX][NNY][2],
    double uF[NNX][NNY][2], double uFN[NNX][NNY][2],
    double pi[NNX][NNY]   , double piN[NNX][NNY], double p[NNX][NNY],
    int    ff[NNX][NNY]   , int    bf[NNX][NNY][2]
) {
    double t        = 0;
    double E        = 2.5E-3;
    double diver    = calcdiv(uF, u, ff, bf);
    double maxdiver = 0;
    int    it;
    printf("\rt = %5.8lf, poi = %5d, max div = %5.8lf", t, 0, diver);
    fflush(stdout);
    while(t < Tr) {
        diver    = DBL_MAX;
        maxdiver = 0;
        prediction(u, uN, uF, p, ff, bf);
        while(diver > E) {
            initPi(pi);
            it    = poisson(pi, piN, uF, u, ff, bf);
            diver = projection(u, uF, uFN, pi, p, ff, bf);
            if (diver > maxdiver) {
                maxdiver = diver;
            }
            printf("\rt = %5.8lf, poi = %5d, div = %5.8lf, max div = %5.8lf", t, it, diver, maxdiver);
            fflush(stdout);
        }
        t += dt;
    }
    printf("\n");
}

void calcgrid(
    double u[NNX][NNY][2] , double uF[NNX][NNY][2], double uG[NGX][NGY][2],
    double p[NNX][NNY]    , double pG[NGX][NGY]   ,
    int    ff[NNX][NNY]   , int    bf[NNX][NNY][2]
) {
    int    OFST = GHST - 1;
    double ufn, ufs, vfe, vfw, vfn, vfs, ufe, ufw; // cell face velocities at each direction
    double une, unw, use, usw, vne, vnw, vse, vsw; // cell center velocities at each direction
    double pne, pnw, pse, psw;                     // pressure at each direction
    double ugc, vgc;                               // velocity at grid point
    double pgc;                                    // pressure at grid point
    int    fne, fnw, fse, fsw;
    int    bfn, bfs, bfe, bfw;
    // int    cfnu, cfsu, cfeu, cfwu;
    // int    cfnv, cfsv, cfev, cfwv;
    for (int i = 0; i < NGX; i ++) {
        for (int k = 0; k < NGY; k ++) {
            une = u[i + OFST + 1][k + OFST + 1][u_];
            unw = u[i + OFST    ][k + OFST + 1][u_];
            use = u[i + OFST + 1][k + OFST    ][u_];
            usw = u[i + OFST    ][k + OFST    ][u_];
            vne = u[i + OFST + 1][k + OFST + 1][v_];
            vnw = u[i + OFST    ][k + OFST + 1][v_];
            vse = u[i + OFST + 1][k + OFST    ][v_];
            vsw = u[i + OFST    ][k + OFST    ][v_];

            ufn = uF[i + OFST    ][k + OFST + 1][u_];
            ufs = uF[i + OFST    ][k + OFST    ][u_];
            vfe = uF[i + OFST + 1][k + OFST    ][v_];
            vfw = uF[i + OFST    ][k + OFST    ][v_];

            pne = p[i + OFST + 1][k + OFST + 1];
            pnw = p[i + OFST    ][k + OFST + 1];
            pse = p[i + OFST + 1][k + OFST    ];
            psw = p[i + OFST    ][k + OFST    ];

            fne = ff[i + OFST + 1][k + OFST + 1];
            fnw = ff[i + OFST    ][k + OFST + 1];
            fse = ff[i + OFST + 1][k + OFST    ];
            fsw = ff[i + OFST    ][k + OFST    ];

            bfn = bf[i + OFST    ][k + OFST + 1][x_];
            bfs = bf[i + OFST    ][k + OFST    ][x_];
            bfe = bf[i + OFST + 1][k + OFST    ][y_];
            bfw = bf[i + OFST    ][k + OFST    ][y_];


            if (fne + fnw + fse + fsw == 4 * FLUID) {
                pgc =(fne * pne + fnw * pnw + fse * pse + fsw * psw) / (fne + fnw + fse + fsw);
                ugc = 0.5 * (ufn + ufs);
                vgc = 0.5 * (vfe + vfw);
            }
            else if (fne + fnw + fse + fsw == 3 * FLUID) {
                pgc =(fne * pne + fnw * pnw + fse * pse + fsw * psw) / (fne + fnw + fse + fsw);
                if (fne == NONFLUID) {
                    ufe = bcf(BTU[bfe], use,   0.5 * dy, BCU[bfe][u_]);
                    vfn = bcf(BTU[bfn], vnw,   0.5 * dx, BCU[bfn][v_]);
                    ugc = 0.5 * (ufn + ufe);
                    vgc = 0.5 * (vfn + vfe);
                }
                else if (fnw == NONFLUID) {
                    ufw = bcf(BTU[bfw], usw,   0.5 * dy, BCU[bfw][u_]);
                    vfn = bcf(BTU[bfn], vne, - 0.5 * dx, BCU[bfn][v_]);
                    ugc = 0.5 * (ufn + ufw);
                    vgc = 0.5 * (vfn + vfw);
                }
                else if (fsw == NONFLUID) {
                    ufw = bcf(BTU[bfw], unw, - 0.5 * dy, BCU[bfw][u_]);
                    vfs = bcf(BTU[bfs], vse, - 0.5 * dx, BCU[bfs][v_]);
                    ugc = 0.5 * (ufw + ufs);
                    vgc = 0.5 * (vfw + vfs);
                }
                else if (fse == NONFLUID) {
                    ufe = bcf(BTU[bfe], une, - 0.5 * dy, BCU[bfe][u_]);
                    vfs = bcf(BTU[bfs], vsw,   0.5 * dx, BCU[bfs][v_]);
                    ugc = 0.5 * (ufe + ufs);
                    vgc = 0.5 * (vfe + vfs);
                }
            }
            else if (fne + fnw + fse + fsw == 2 * FLUID) {
                pgc =(fne * pne + fnw * pnw + fse * pse + fsw * psw) / (fne + fnw + fse + fsw);
                if (fne + fse == NONFLUID) {
                    ugc = 0.5 * (ufn + ufs);
                    vgc = bcf(BTU[bfn], vfw,   0.5 * dx, BCU[bfn][v_]);
                }
                else if (fnw + fsw == NONFLUID) {
                    ugc = 0.5 * (ufn + ufs);
                    vgc = bcf(BTU[bfs], vfe, - 0.5 * dx, BCU[bfs][v_]);
                }
                else if (fnw + fne == NONFLUID) {
                    ugc = bcf(BTU[bfw], ufs,   0.5 * dy, BCU[bfw][u_]);
                    vgc = 0.5 * (vfe + vfw);
                }
                else if (fsw + fse == NONFLUID) {
                    ugc = bcf(BTU[bfe], ufn, - 0.5 * dy, BCU[bfe][u_]);
                    vgc = 0.5 * (vfe + vfw);
                }
                else {
                    ugc = 0.5 * (ufn + ufs);
                    vgc = 0.5 * (vfe + vfw);
                }
            }
            else if (fne + fnw + fse + fsw == FLUID) {
                pgc =(fne * pne + fnw * pnw + fse * pse + fsw * psw) / (fne + fnw + fse + fsw);
                if (fne == FLUID) {
                    ufe = bcf(BTU[bfe], une, - 0.5 * dy, BCU[bfe][u_]);
                    vfn = bcf(BTU[bfn], vne, - 0.5 * dx, BCU[bfn][v_]);
                    ugc = 0.5 * (ufn + ufe);
                    vgc = 0.5 * (vfn + vfe);
                }
                else if (fnw == FLUID) {
                    ufw = bcf(BTU[bfw], unw, - 0.5 * dy, BCU[bfw][u_]);
                    vfn = bcf(BTU[bfn], vnw,   0.5 * dx, BCU[bfn][v_]);
                    ugc = 0.5 * (ufn + ufw);
                    vgc = 0.5 * (vfn + vfw);
                }
                else if (fsw == FLUID) {
                    ufw = bcf(BTU[bfw], usw,   0.5 * dy, BCU[bfw][u_]);
                    vfs = bcf(BTU[bfs], vsw,   0.5 * dx, BCU[bfs][v_]);
                    ugc = 0.5 * (ufw + ufs);
                    vgc = 0.5 * (vfw + vfs);
                }
                else if (fse == FLUID) {
                    ufe = bcf(BTU[bfe], use,   0.5 * dy, BCU[bfe][u_]);
                    vfs = bcf(BTU[bfs], vse, - 0.5 * dx, BCU[bfs][v_]);
                    ugc = 0.5 * (ufe + ufs);
                    vgc = 0.5 * (vfe + vfs);
                }
            }
            else if (fne + fnw + fse + fsw == NONFLUID) {
                pgc = 0.25 * (pne + pnw + pse + psw);
                ugc = 0.5 * (ufn + ufs);
                vgc = 0.5 * (vfe + vfw);
            }

            uG[i][k][u_] = ugc;
            uG[i][k][v_] = vgc;
            pG[i][k]     = pgc;
        }
    }
}

void o2f(double uG[NGX][NGY][2], double pG[NGX][NGY]) {
    FILE  *fo;
    char   fname[128];
    double xpos, ypos;
    sprintf(fname, "UVP_Re%d.ns2d.csv", int(Re), int(Tr));
    fo = fopen(fname, "w+t");

    if ( fo == NULL ) {
        printf("\nERROR when opening file\n");
    }
    else {
        fprintf(fo, "x,y,z,u,v,p\n");
        for (int k = 0; k < NGY; k ++) {
            for (int i = 0; i < NGX; i ++) {
                xpos = i * dx;
                ypos = k * dy;
                fprintf(fo, "%5.8lf,%5.8lf,%5.8lf,%5.8lf,%5.8lf,%5.8lf\n", xpos, ypos, 0.0, uG[i][k][u_], uG[i][k][v_], pG[i][k]);
            }
        }
        fclose(fo);
    }
}

int main(int argc, char ** argv) {
    if (argc < 3) {
        printf("ns2d Re time\n");
        return 0;
    }
    Re = strtod(argv[1], NULL);
    Tr = strtod(argv[2], NULL);
    dt = min(CFL * min(dx, dy) * 0.5, 0.5 * Re * dx * dx * dy * dy / (2 * (dx * dx + dy * dy)));
    printf("Re = %lf, T = %lf, dt = %lf, dtau = %lf\n", Re, Tr, dt, dtau);

    setupGeom(FF, BF);
    #pragma acc enter data copyin(U, UN, UF, UFN, Pi, PiN, P, FF, BF, BTU, BTP, BCU, BCP)
    init(U, UF, P, Pi, FF, BF);
    ns2d(U, UN, UF, UFN, Pi, PiN, P, FF, BF);
    #pragma acc exit data copyout(U, UF, P)
    calcgrid(U, UF, UG, P, PG, FF, BF);
    o2f(UG, PG);
    

    return 0;
}
