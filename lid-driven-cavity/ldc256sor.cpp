#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <algorithm>

using namespace std;

#define _u         0
#define _v         1
#define _p         2
#define _pi        3

#define _x         0
#define _y         1

#define _dirichlet 0
#define _neumann   1

const int NCX      = 256;
const int NCY      = 256;
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
const double CFL   = 0.5;
const double DIFFC = 0.5;

// solver parameters
const int EPSILON3 = 1;
const int EPSILON1 = 0;
const double KAPPA = 1.0 / 3.0;
const double BETA  = (3 - KAPPA) / (1 - KAPPA);
const double OMEGA = 1.8;

double Re;  // Reynolds number
double dt;  // timestep length
double Tr;  // simulation time span

// data vectors
// note that Pi is the potential and P is the pressure
double U[NNX][NNY][2], UA[NNX][NNY][2], UF[NNX][NNY][2], UFN[NNX][NNY][2], UG[NGX][NGY][2]; // velocity
double Pi[NNX][NNY], P[NNX][NNY], PG[NGX][NGY];                                             // pressure
double DIV[NNX][NNY];
/**********************************************************************************/
/********************************** geometry setup ********************************/
// cell types
const int FLUID        = 1;
const int NONFLUID     = 0;
const int BLACK        = 0;
const int RED          = 1;

int FF[NNX][NNY];    // fluid flag
int BF[NNX][NNY][2]; // boundary flag

// boundary indices
const int    NONE      = 0;
const int    WALL      = 1;
const int    LID       = 2;
// boundary types
const int    DIRICHLET = 1;
const int    NEUMANN   = 2;
// boundary condition values
const double UWALL     = 0;
const double VWALL     = 0;
const double PWALL     = 0;
const double PIWALL    = 0;
const double ULID      = 1;
const double VLID      = 0;
const double PLID      = 0;
const double PILID     = 0;

const int BT[3][4][2]  = {
//    u          v          p          π
//    D  N       D  N       D  N       D  N
    {{0, 0}   , {0, 0}   , {0, 0}   , {0, 0}   }, // no boundary
    {{1, 0}   , {1, 0}   , {0, 1}   , {0, 1}   }, // wall
    {{1, 0}   , {1, 0}   , {0, 1}   , {0, 1}   }, // lid

//  {NONE     , NONE     , NONE     , NONE     }, // no boundary
//  {DIRICHLET, DIRICHLET, NEUMANN  , NEUMANN  }, // wall
//  {DIRICHLET, DIRICHLET, NEUMANN  , NEUMANN  }, // lid
};

double BC(int index, int variable, double i, double k, double t) {
    double bc[3][4] = {
    //   u      v      p      π
        {0.0  , 0.0  , 0.0  , 0.0   }, // no boundary
        {UWALL, VWALL, PWALL, PIWALL}, // wall
        {ULID , VLID , PLID , PILID }, // lid

    };

    return bc[index][variable];
}

void setupGeom(int ff[NNX][NNY], int bf[NNX][NNY][2]) {
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            ff[i][k] = NONFLUID;
        }
    }
    // the cavity
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            ff[i][k] = FLUID;
        }
    }
    
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            bf[i][k][_x] = NONE;
            bf[i][k][_y] = NONE;
        }
    }
    // bottom wall
    for (int i = WBOUND; i < EBOUND; i ++) {
        bf[i][SBOUND - 1][_y] = WALL;
    }
    // left wall
    for (int k = SBOUND; k < NBOUND; k ++) {
        bf[WBOUND - 1][k][_x] = WALL;
    }
    // right wall
    for (int k = SBOUND; k < NBOUND; k ++) {
        bf[EBOUND - 1][k][_x] = WALL;
    }
    // lid
    for (int i = WBOUND; i < EBOUND; i ++) {
        bf[i][NBOUND - 1][_y] = LID;
    }
}
/**********************************************************************************/
/******************************* fetch boundary flags *****************************/
void bE1(
    int i, int k, 
    int bf[NNX][NNY][2], 
    int &ee1, int &ew1, int &en1, int &es1
) {
    ee1 = bf[i    ][k    ][_x];
    en1 = bf[i    ][k    ][_y];
    ew1 = bf[i - 1][k    ][_x];
    es1 = bf[i    ][k - 1][_y];
}

void bE2(
    int i, int k,
    int bf[NNX][NNY][2], 
    int &ee1, int &ew1, int &en1, int &es1, 
    int &ee2, int &ew2, int &en2, int &es2
) {
    ee1 = bf[i    ][k    ][_x];
    en1 = bf[i    ][k    ][_y];
    ew1 = bf[i - 1][k    ][_x];
    es1 = bf[i    ][k - 1][_y];

    ee2 = bf[i + 1][k    ][_x];
    en2 = bf[i    ][k + 1][_y];
    ew2 = bf[i - 2][k    ][_x];
    es2 = bf[i    ][k - 2][_y];
}
/**********************************************************************************/
/**************** extrapolation according to boundary conditions ******************/
double bcc(const int type[2], double base, double direction, double bcvalue) {
    if (type[_dirichlet]) {
        return (2 * bcvalue - base);
    }
    if (type[_neumann  ]) {
        return (base + direction * bcvalue);
    }
    return 0.0;
}

double bcf(const int type[2], double base, double direction, double bcvalue) {
    if (type[_dirichlet]) {
        return bcvalue;
    }
    if (type[_neumann  ]) {
        return (base + direction * 0.5 * bcvalue);
    }
    return 0.0;
}
/**********************************************************************************/
/****************************** calculate ||div(U)|| ******************************/
double calcdiv(
    double  uF[NNX][NNY][2], double  u[NNX][NNY][2], 
    int     ff[NNX][NNY]   , int    bf[NNX][NNY][2],
    double div[NNX][NNY]   , double  t
) {
    double tdiver = 0;
    #pragma acc kernels loop independent reduction(+:tdiver) collapse(2) present(uF, u, ff, bf, BT, div) copy(tdiver)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double ufe, ufw, vfn, vfs; // velocity at each direction
                double uc0, vc0;           // velocity at cell center
                double dudx, dvdy;         // velocity gradient at cell center
                double diver;
                int    ee1, ew1, en1, es1; // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0 =  u[i    ][k    ][_u];
                vc0 =  u[i    ][k    ][_v];

                ufe = uF[i    ][k    ][_u];
                ufw = uF[i - 1][k    ][_u];
                vfn = uF[i    ][k    ][_v];
                vfs = uF[i    ][k - 1][_v];

                ufe = (ee1 == NONE)? ufe : bcf(BT[ee1][_u], uc0,   dx, BC(ee1, _u, i + 0.5, k      , t));
                ufw = (ew1 == NONE)? ufw : bcf(BT[ew1][_u], uc0, - dx, BC(ew1, _u, i - 0.5, k      , t));
                vfn = (en1 == NONE)? vfn : bcf(BT[en1][_v], vc0,   dy, BC(en1, _v, i      , k + 0.5, t));
                vfs = (es1 == NONE)? vfs : bcf(BT[es1][_v], vc0, - dy, BC(es1, _v, i      , k - 0.5, t));

                dudx      = (ufe - ufw) / dx;
                dvdy      = (vfn - vfs) / dy;
                diver     = dudx + dvdy;
                tdiver    = tdiver + diver * diver;
                div[i][k] = diver;
            }
        }
    }
    return sqrt(tdiver);
}
/**********************************************************************************/
/************** interpolate cell face values from cell center values **************/
void uc2uf(
    double uC[NNX][NNY][2], double uF[NNX][NNY][2], 
    int    ff[NNX][NNY]   , int    bf[NNX][NNY][2],
    int    t              , bool   ifbc
) {
    #pragma acc kernels loop independent collapse(2) present(uC, uF, ff, bf, BT)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double ue1, uc0, uw1, vn1, vc0, vs1; // velocity at each direction
                double ufe, ufw, vfn, vfs;           // velocity at cell faces
                int    ee1, ew1, en1, es1;           // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0 = uC[i    ][k    ][_u];
                vc0 = uC[i    ][k    ][_v];

                ue1 = uC[i + 1][k    ][_u];
                uw1 = uC[i - 1][k    ][_u];
                vn1 = uC[i    ][k + 1][_v];
                vs1 = uC[i    ][k - 1][_v];

                if(ifbc) {
                    ue1 = (ee1 == NONE)? ue1 : bcc(BT[ee1][_u], uc0,   dx, BC(ee1, _u, i + 0.5, k      , t));
                    uw1 = (ew1 == NONE)? uw1 : bcc(BT[ew1][_u], uc0, - dx, BC(ew1, _u, i - 0.5, k      , t));
                    vn1 = (en1 == NONE)? vn1 : bcc(BT[en1][_v], vc0,   dy, BC(en1, _v, i      , k + 0.5, t));
                    vs1 = (es1 == NONE)? vs1 : bcc(BT[es1][_v], vc0, - dy, BC(es1, _v, i      , k - 0.5, t));
                }
                
                ufe = 0.5 * (uc0 + ue1);
                ufw = 0.5 * (uw1 + uc0);
                vfn = 0.5 * (vc0 + vn1);
                vfs = 0.5 * (vs1 + vc0);

                if(ifbc) {
                    ufe = (ee1 == NONE)? ufe : bcf(BT[ee1][_u], uc0,   dx, BC(ee1, _u, i + 0.5, k      , t));
                    ufw = (ew1 == NONE)? ufw : bcf(BT[ew1][_u], uc0, - dx, BC(ew1, _u, i - 0.5, k      , t));
                    vfn = (en1 == NONE)? vfn : bcf(BT[en1][_v], vc0,   dy, BC(en1, _v, i      , k + 0.5, t));
                    vfs = (es1 == NONE)? vfs : bcf(BT[es1][_v], vc0, - dy, BC(es1, _v, i      , k - 0.5, t));

                }
                
                uF[i    ][k    ][_u] = ufe;
                uF[i - 1][k    ][_u] = ufw;
                uF[i    ][k    ][_v] = vfn;
                uF[i    ][k - 1][_v] = vfs;
            }
        }
    }
}
/**********************************************************************************/
/******************************** initialization **********************************/
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
    double  u[NNX][NNY][2], double uF[NNX][NNY][2], 
    double  p[NNX][NNY]   , double pi[NNX][NNY]   , 
    int    ff[NNX][NNY]   , int    bf[NNX][NNY][2]
) {
    #pragma acc kernels loop collapse(2) present(u)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            u[i][k][_u] = 0;
            u[i][k][_v] = 0;
        }
    }
    // #pragma acc kernels loop independent collapse(2) present(u, step)
    // for (int i = WBOUND; i < EBOUND; i ++) {
    //     for (int k = step[1][_y]; k < NBOUND; k ++) {
    //         u[i][k][_u] = BC(INFLOW, _u, i, k, 0.0);
    //         // u[i][k][_u] = UINFLOW;
    //         u[i][k][_v] = 0;
    //     }
    // }
    uc2uf(u, uF, ff, bf, 0.0, true);
    initPi(pi);
    #pragma acc kernels loop collapse(2) present(p)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            p[i][k] = 0;
        }
    }
}
/**********************************************************************************/
/************************************* copy ***************************************/
// copy U to UN or UF to UFN
void cpUt(double ut[NNX][NNY][2], double utN[NNX][NNY][2]) {
    #pragma acc kernels loop independent collapse(2) present(ut, utN)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            utN[i][k][_u] = ut[i][k][_u];
            utN[i][k][_v] = ut[i][k][_v];
        }
    }
}

// copy Pi to PiN or P to PN
void cpPt(double pt[NNX][NNY], double ptN[NNX][NNY]) {
    #pragma acc kernels loop independent collapse(2) present(pt, ptN)
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            ptN[i][k] = pt[i][k];
        }
    }
}
/**********************************************************************************/
/************************************* flux ***************************************/
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
double R(double dm, double dp) {
    double numer = dm;
    double denom = dp;
    double e     = 1E-9;
    numer = (numer < 0)? numer - e : numer + e;
    denom = (denom < 0)? denom - e : denom + e;
    return numer / denom;
}
// albada1 limiter
double albada1(double r) {
    return (r + r * r) / (1 + r * r);
}
// albada2 limiter
double albada2(double r) {
    return 2 * r / (1 + r * r);
}
// superbee limiter
double superbee(double r) {
    double tmp = max(min(2.0 * r, 1.0), min(r, 2.0));
    return max(0.0, tmp);
}
// minmod limiter
double minmod(double r) {
    return max(0.0, min(1.0, r));
}
// interpolated value φR at cell face: φl2 φl1 | φr1 φr2
// φR = φr1 - ε/4 * limiter(r) * ((1 - κ) * Δ+ + (1 + κ) * Δ-)
// d+ = φr2 - φr1, d- = φr1 - φl1
// r = d- / d+
double interpolR(double l1, double r1, double r2, int eps) {
    double dm = r1 - l1;
    double dp = r2 - r1;
    return r1 - eps * 0.25 * minmod(R(dm, dp)) * ((1 - KAPPA) * dp + (1 + KAPPA) * dm);
}
// interpolated value φL at cell face: φl2 φl1 | φr1 φr2
// φL = φl1 + ε/4 *  limiter(r) * ((1 - κ) * Δ- + (1 + κ) * Δ+)
// d+ = φr1 - φl1, d- = φl1 - φl2
// r = d- / d+
double interpolL(double l2, double l1, double r1, int eps) {
    double dm = l1 - l2;
    double dp = r1 - l1;
    return l1 + eps * 0.25 * minmod(R(dm, dp)) * ((1 - KAPPA) * dm + (1 + KAPPA) * dp);
}
// numerical flux at cell face: φl2 φl1 | φr1 φr2
// f~ = 1/2 * ((fR + fL) - |u@face| * (φR - φL))
// fRL = u@face * φRL
double flux(double l2, double l1, double r1, double r2, double uf, int epsL, int epsR) {
    double phfR = interpolR(l1, r1, r2, epsR);
    double phfL = interpolL(l2, l1, r1, epsL);
    double flxR = uf * phfR;
    double flxL = uf * phfL;
    return 0.5 * (flxR + flxL - abs(uf) * (phfR - phfL));
}
/**********************************************************************************/
/************************************* SMAC ***************************************/
// u* = u + Δt * (- Convection + Viscoucity)
// uP@center = u*@center - Δt * grad(p)@center
// uP@face = average(u*@center) - Δt * grad(p)@face
void prediction(
    double  u[NNX][NNY][2], double  uA[NNX][NNY][2] , 
    double uF[NNX][NNY][2], double uFN[NNX][NNY][2],
    double  p[NNX][NNY]   , int     ff[NNX][NNY]    , 
    int    bf[NNX][NNY][2], double   t
) {
    // calculate u* (intermediate value without pressure gradient) at cell center
    #pragma acc kernels loop independent collapse(2) present(u, uA, uF, ff, bf, BT)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double fue, fuw, fun, fus, fve, fvw, fvn, fvs;      // flux variables
                double uw2, uw1, us2, us1, uc0, un1, un2, ue1, ue2; // velocity at each direction
                double vw2, vw1, vs2, vs1, vc0, vn1, vn2, ve1, ve2; // velocity at each direction
                double ufe, ufw, vfn, vfs;                          // velocity at cell faces
                double duudx, dvudy, duvdx, dvvdy;                  // 1st derivatives
                double ddudxx, ddudyy, ddvdxx, ddvdyy;              // 2nd derivatives
                int    ee1, ew1, en1, es1;                          // first level boundary indicator
                int    ee2, ew2, en2, es2;                          // second level boundary indicator
                int    epseR, epswR, epsnR, epssR;
                int    epseL, epswL, epsnL, epssL;
                bE2(i, k, bf, ee1, ew1, en1, es1, ee2, ew2, en2, es2);

                uc0 =  u[i    ][k    ][_u];
                vc0 =  u[i    ][k    ][_v];

                ue1 =  u[i + 1][k    ][_u];
                uw1 =  u[i - 1][k    ][_u];
                un1 =  u[i    ][k + 1][_u];
                us1 =  u[i    ][k - 1][_u];
                ve1 =  u[i + 1][k    ][_v];
                vw1 =  u[i - 1][k    ][_v];
                vn1 =  u[i    ][k + 1][_v];
                vs1 =  u[i    ][k - 1][_v];

                ue2 =  u[i + 2][k    ][_u];
                uw2 =  u[i - 2][k    ][_u];
                un2 =  u[i    ][k + 2][_u];
                us2 =  u[i    ][k - 2][_u];
                ve2 =  u[i + 2][k    ][_v];
                vw2 =  u[i - 2][k    ][_v];
                vn2 =  u[i    ][k + 2][_v];
                vs2 =  u[i    ][k - 2][_v];

                ufe = uF[i    ][k    ][_u];
                ufw = uF[i - 1][k    ][_u];
                vfn = uF[i    ][k    ][_v];
                vfs = uF[i    ][k - 1][_v];

                ue1 = (ee1 == NONE)? ue1 : bcc(BT[ee1][_u], uc0,   dx, BC(ee1, _u, i + 0.5, k      , t));
                uw1 = (ew1 == NONE)? uw1 : bcc(BT[ew1][_u], uc0, - dx, BC(ew1, _u, i - 0.5, k      , t));
                un1 = (en1 == NONE)? un1 : bcc(BT[en1][_u], uc0,   dy, BC(en1, _u, i      , k + 0.5, t));
                us1 = (es1 == NONE)? us1 : bcc(BT[es1][_u], uc0, - dy, BC(es1, _u, i      , k - 0.5, t));
                ve1 = (ee1 == NONE)? ve1 : bcc(BT[ee1][_v], vc0,   dx, BC(ee1, _v, i + 0.5, k      , t));
                vw1 = (ew1 == NONE)? vw1 : bcc(BT[ew1][_v], vc0, - dx, BC(ew1, _v, i - 0.5, k      , t));
                vn1 = (en1 == NONE)? vn1 : bcc(BT[en1][_v], vc0,   dy, BC(en1, _v, i      , k + 0.5, t));
                vs1 = (es1 == NONE)? vs1 : bcc(BT[es1][_v], vc0, - dy, BC(es1, _v, i      , k - 0.5, t));

                ue2 = (ee2 == NONE)? ue2 : bcc(BT[ee2][_u], ue1,   dx, BC(ee2, _u, i + 1.5, k      , t));
                uw2 = (ew2 == NONE)? uw2 : bcc(BT[ew2][_u], uw1, - dx, BC(ew2, _u, i - 1.5, k      , t));
                un2 = (en2 == NONE)? un2 : bcc(BT[en2][_u], un1,   dy, BC(en2, _u, i      , k + 1.5, t));
                us2 = (es2 == NONE)? us2 : bcc(BT[es2][_u], us1, - dy, BC(es2, _u, i      , k - 1.5, t));
                ve2 = (ee2 == NONE)? ve2 : bcc(BT[ee2][_v], ve1,   dx, BC(ee2, _v, i + 1.5, k      , t));
                vw2 = (ew2 == NONE)? vw2 : bcc(BT[ew2][_v], vw1, - dx, BC(ew2, _v, i - 1.5, k      , t));
                vn2 = (en2 == NONE)? vn2 : bcc(BT[en2][_v], vn1,   dy, BC(en2, _v, i      , k + 1.5, t));
                vs2 = (es2 == NONE)? vs2 : bcc(BT[es2][_v], vs1, - dy, BC(es2, _v, i      , k - 1.5, t));

                ufe = (ee1 == NONE)? ufe : bcf(BT[ee1][_u], uc0,   dx, BC(ee1, _u, i + 0.5, k      , t));
                ufw = (ew1 == NONE)? ufw : bcf(BT[ew1][_u], uc0, - dx, BC(ew1, _u, i - 0.5, k      , t));
                vfn = (en1 == NONE)? vfn : bcf(BT[en1][_v], vc0,   dy, BC(en1, _v, i      , k + 0.5, t));
                vfs = (es1 == NONE)? vfs : bcf(BT[es1][_v], vc0, - dy, BC(es1, _v, i      , k - 0.5, t));

                epseL = EPSILON3;
                epseR = (ee1 == NONE)? EPSILON3 : EPSILON1;
                epswL = (ew1 == NONE)? EPSILON3 : EPSILON1;
                epswR = EPSILON3;
                epsnL = EPSILON3;
                epsnR = (en1 == NONE)? EPSILON3 : EPSILON1;
                epssL = (es1 == NONE)? EPSILON3 : EPSILON1;
                epssR = EPSILON3;

                // epseL = EPSILON1;
                // epseR = EPSILON1;
                // epswL = EPSILON1;
                // epswR = EPSILON1;
                // epsnL = EPSILON1;
                // epsnR = EPSILON1;
                // epssL = EPSILON1;
                // epssR = EPSILON1;

                fue = flux(uw1, uc0, ue1, ue2, ufe, epseL, epseR);
                fuw = flux(uw2, uw1, uc0, ue1, ufw, epswL, epswR);
                fun = flux(us1, uc0, un1, un2, vfn, epsnL, epsnR);
                fus = flux(us2, us1, uc0, un1, vfs, epssL, epssR);
                fve = flux(vw1, vc0, ve1, ve2, ufe, epseL, epseR);
                fvw = flux(vw2, vw1, vc0, ve1, ufw, epswL, epswR);
                fvn = flux(vs1, vc0, vn1, vn2, vfn, epsnL, epsnR);
                fvs = flux(vs2, vs1, vc0, vn1, vfs, epssL, epssR);

                duudx  = (fue - fuw) / dx;
                dvudy  = (fun - fus) / dy;
                duvdx  = (fve - fvw) / dx;
                dvvdy  = (fvn - fvs) / dy;
                ddudxx = (uw1 - 2 * uc0 + ue1) / (dx * dx);
                ddudyy = (us1 - 2 * uc0 + un1) / (dy * dy);
                ddvdxx = (vw1 - 2 * vc0 + ve1) / (dx * dx);
                ddvdyy = (vs1 - 2 * vc0 + vn1) / (dy * dy);

                uc0    = uc0 + dt * (- duudx - dvudy + (ddudxx + ddudyy) / Re);
                vc0    = vc0 + dt * (- duvdx - dvvdy + (ddvdxx + ddvdyy) / Re);

                uA[i][k][_u] = uc0;
                uA[i][k][_v] = vc0;
            }
        }
    }
    uc2uf(uA, uF, ff, bf, t, false);
    cpUt(uF, uFN);
    // calculate uP (prediction value with pressure gradient) at cell faces and cell center
    #pragma acc kernels loop independent collapse(2) present(u, uF, uA, uFN, p, ff, bf, BT)
    for(int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if(ff[i][k] == FLUID) {
                double pe1, pw1, pc0, pn1, ps1;        // pressure at each direction
                double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // pressure derivatives at cell faces
                double dpdxc0, dpdyc0;                 // pressure derivatives at cell center
                double ufe, ufw, vfn, vfs;             // velocity at cell faces
                double uc0, vc0;                       // velocity at cell center
                int    ee1, ew1, en1, es1;             // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0 =  uA[i    ][k    ][_u];
                vc0 =  uA[i    ][k    ][_v];

                ufe = uFN[i    ][k    ][_u];
                ufw = uFN[i - 1][k    ][_u];
                vfn = uFN[i    ][k    ][_v];
                vfs = uFN[i    ][k - 1][_v];

                pc0 =   p[i    ][k    ];
                pe1 =   p[i + 1][k    ];
                pw1 =   p[i - 1][k    ];
                pn1 =   p[i    ][k + 1];
                ps1 =   p[i    ][k - 1];

                pe1 = (ee1 == NONE)? pe1 : bcc(BT[ee1][_p], pc0,   dx, BC(ee1, _p, i + 0.5, k      , t));
                pw1 = (ew1 == NONE)? pw1 : bcc(BT[ew1][_p], pc0, - dx, BC(ew1, _p, i - 0.5, k      , t));
                pn1 = (en1 == NONE)? pn1 : bcc(BT[en1][_p], pc0,   dy, BC(en1, _p, i      , k + 0.5, t));
                ps1 = (es1 == NONE)? ps1 : bcc(BT[es1][_p], pc0, - dy, BC(es1, _p, i      , k - 0.5, t));

                dpdxfe = (pe1 - pc0) / dx;
                dpdxfw = (pc0 - pw1) / dx;
                dpdyfn = (pn1 - pc0) / dy;
                dpdyfs = (pc0 - ps1) / dy;
                dpdxc0 = 0.5 * (dpdxfe + dpdxfw);
                dpdyc0 = 0.5 * (dpdyfn + dpdyfs);

                uc0 = uc0 - dt * dpdxc0;
                vc0 = vc0 - dt * dpdyc0;
                ufe = ufe - dt * dpdxfe;
                ufw = ufw - dt * dpdxfw;
                vfn = vfn - dt * dpdyfn;
                vfs = vfs - dt * dpdyfs;

                ufe  = (ee1 == NONE)? ufe  : bcf(BT[ee1][_u ], uc0 ,   dx, BC(ee1, _u , i + 0.5, k      , t));
                ufw  = (ew1 == NONE)? ufw  : bcf(BT[ew1][_u ], uc0 , - dx, BC(ew1, _u , i - 0.5, k      , t));
                vfn  = (en1 == NONE)? vfn  : bcf(BT[en1][_v ], vc0 ,   dy, BC(en1, _v , i      , k + 0.5, t));
                vfs  = (es1 == NONE)? vfs  : bcf(BT[es1][_v ], vc0 , - dy, BC(es1, _v , i      , k - 0.5, t));

                uF[i    ][k    ][_u] = ufe;
                uF[i - 1][k    ][_u] = ufw;
                uF[i    ][k    ][_v] = vfn;
                uF[i    ][k - 1][_v] = vfs;

                u[i][k][_u]          = uc0;
                u[i][k][_v]          = vc0;
            }
        }
    }
}
// div(grad(π)) = div(uP) / Δt
int poisson(
    double pi[NNX][NNY]   , double uF[NNX][NNY][2], 
    double  u[NNX][NNY][2], int    ff[NNX][NNY]   , 
    int    bf[NNX][NNY][2], double &r             , 
    double  t
) {
    double E     = 1E-6;
    int    MAXIT = 1000;
    bool   con   = false;
    int    it    = 0;
    double R     = 0;
    while (con == false) {
        R   = 0;
        con = true;
        #pragma acc kernels loop independent reduction(+:R) collapse(2) present(pi, uF, u, ff, bf, BT) copy(R)
        for (int i = WBOUND; i < EBOUND; i ++) {
            for (int k = SBOUND; k < NBOUND; k ++) {
                if (ff[i][k] == FLUID && (i + k) % 2 == BLACK) {
                    double ufe, ufw, vfn, vfs;           // velocity at cell faces
                    double uc0, vc0;                     // velocity at cell center
                    double pie1, piw1, pic0, pin1, pis1; // pressure potential at each direction
                    double pica;                         // πc*
                    double dudx, dvdy;                   // velocity gradient at cell center
                    double ae1, aw1, ac0, an1, as1;      // coefficient
                    double psi, res;                     //
                    int    ee1, ew1, en1, es1;           // first level boundary indicator
                    bE1(i, k, bf, ee1, ew1, en1, es1);

                    uc0  =  u[i    ][k    ][_u];
                    vc0  =  u[i    ][k    ][_v];
                    ufe  = uF[i    ][k    ][_u];
                    ufw  = uF[i - 1][k    ][_u];
                    vfn  = uF[i    ][k    ][_v];
                    vfs  = uF[i    ][k - 1][_v];

                    pic0 = pi[i    ][k    ];
                    pie1 = pi[i + 1][k    ];
                    piw1 = pi[i - 1][k    ];
                    pin1 = pi[i    ][k + 1];
                    pis1 = pi[i    ][k - 1];

                    ufe  = (ee1 == NONE)? ufe  : bcf(BT[ee1][_u ], uc0 ,   dx, BC(ee1, _u , i + 0.5, k      , t));
                    ufw  = (ew1 == NONE)? ufw  : bcf(BT[ew1][_u ], uc0 , - dx, BC(ew1, _u , i - 0.5, k      , t));
                    vfn  = (en1 == NONE)? vfn  : bcf(BT[en1][_v ], vc0 ,   dy, BC(en1, _v , i      , k + 0.5, t));
                    vfs  = (es1 == NONE)? vfs  : bcf(BT[es1][_v ], vc0 , - dy, BC(es1, _v , i      , k - 0.5, t));

                    pie1 = (ee1 == NONE)? pie1 : bcc(BT[ee1][_pi], pic0,   dx, BC(ee1, _pi, i + 0.5, k      , t));
                    piw1 = (ew1 == NONE)? piw1 : bcc(BT[ew1][_pi], pic0, - dx, BC(ew1, _pi, i - 0.5, k      , t));
                    pin1 = (en1 == NONE)? pin1 : bcc(BT[en1][_pi], pic0,   dy, BC(en1, _pi, i      , k + 0.5, t));
                    pis1 = (es1 == NONE)? pis1 : bcc(BT[es1][_pi], pic0, - dy, BC(es1, _pi, i      , k - 0.5, t));

                    ac0  = - 2.0 / (dx * dx) - 2.0 / (dy * dy);
                    ae1  =   1.0 / (dx * dx);
                    aw1  =   1.0 / (dx * dx);
                    an1  =   1.0 / (dy * dy);
                    as1  =   1.0 / (dy * dy);

                    dudx = (ufe - ufw) / dx;
                    dvdy = (vfn - vfs) / dy;
                    psi  = (dudx + dvdy) / dt;

                    pica = (psi - ae1 * pie1 - aw1 * piw1 - an1 * pin1 - as1 * pis1) / ac0;
                    res  = OMEGA * (pica - pic0);
                    pic0 = pic0 + res;
                    R    = R + res * res;

                    pi[i][k] = pic0;
                }
            }
        }
        #pragma acc kernels loop independent reduction(+:R) collapse(2) present(pi, uF, u, ff, bf, BT) copy(R)
        for (int i = WBOUND; i < EBOUND; i ++) {
            for (int k = SBOUND; k < NBOUND; k ++) {
                if (ff[i][k] == FLUID && (i + k) % 2 == RED) {
                    double ufe, ufw, vfn, vfs;           // velocity at cell faces
                    double uc0, vc0;                     // velocity at cell center
                    double pie1, piw1, pic0, pin1, pis1; // pressure potential at each direction
                    double pica;                         // πc*
                    double dudx, dvdy;                   // velocity gradient at cell center
                    double ae1, aw1, ac0, an1, as1;      // coefficient
                    double psi, res;                     //
                    int    ee1, ew1, en1, es1;           // first level boundary indicator
                    bE1(i, k, bf, ee1, ew1, en1, es1);

                    uc0  =  u[i    ][k    ][_u];
                    vc0  =  u[i    ][k    ][_v];
                    ufe  = uF[i    ][k    ][_u];
                    ufw  = uF[i - 1][k    ][_u];
                    vfn  = uF[i    ][k    ][_v];
                    vfs  = uF[i    ][k - 1][_v];

                    pic0 = pi[i    ][k    ];
                    pie1 = pi[i + 1][k    ];
                    piw1 = pi[i - 1][k    ];
                    pin1 = pi[i    ][k + 1];
                    pis1 = pi[i    ][k - 1];

                    ufe  = (ee1 == NONE)? ufe  : bcf(BT[ee1][_u ], uc0 ,   dx, BC(ee1, _u , i + 0.5, k      , t));
                    ufw  = (ew1 == NONE)? ufw  : bcf(BT[ew1][_u ], uc0 , - dx, BC(ew1, _u , i - 0.5, k      , t));
                    vfn  = (en1 == NONE)? vfn  : bcf(BT[en1][_v ], vc0 ,   dy, BC(en1, _v , i      , k + 0.5, t));
                    vfs  = (es1 == NONE)? vfs  : bcf(BT[es1][_v ], vc0 , - dy, BC(es1, _v , i      , k - 0.5, t));

                    pie1 = (ee1 == NONE)? pie1 : bcc(BT[ee1][_pi], pic0,   dx, BC(ee1, _pi, i + 0.5, k      , t));
                    piw1 = (ew1 == NONE)? piw1 : bcc(BT[ew1][_pi], pic0, - dx, BC(ew1, _pi, i - 0.5, k      , t));
                    pin1 = (en1 == NONE)? pin1 : bcc(BT[en1][_pi], pic0,   dy, BC(en1, _pi, i      , k + 0.5, t));
                    pis1 = (es1 == NONE)? pis1 : bcc(BT[es1][_pi], pic0, - dy, BC(es1, _pi, i      , k - 0.5, t));

                    ac0  = - 2.0 / (dx * dx) - 2.0 / (dy * dy);
                    ae1  =   1.0 / (dx * dx);
                    aw1  =   1.0 / (dx * dx);
                    an1  =   1.0 / (dy * dy);
                    as1  =   1.0 / (dy * dy);

                    dudx = (ufe - ufw) / dx;
                    dvdy = (vfn - vfs) / dy;
                    psi  = (dudx + dvdy) / dt;

                    pica = (psi - ae1 * pie1 - aw1 * piw1 - an1 * pin1 - as1 * pis1) / ac0;
                    res  = OMEGA * (pica - pic0);
                    pic0 = pic0 + res;
                    R    = R + res * res;

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
    r = sqrt(R);
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
// u = uP - Δt * grad(π)
double projection(
    double   u[NNX][NNY][2], double  uF[NNX][NNY][2],
    double uFN[NNX][NNY][2], double  pi[NNX][NNY]   , 
    double   p[NNX][NNY]   , int     ff[NNX][NNY]   , 
    int     bf[NNX][NNY][2], double div[NNX][NNY]   , 
    double   t
) {
    cpUt(uF, uFN);
    #pragma acc kernels loop independent collapse(2) present(u, uF, uFN, pi, ff, bf, BT)
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            if (ff[i][k] == FLUID) {
                double pie1, piw1, pic0, pin1, pis1;       // pressure potential at each direction
                double dpidxfe, dpidxfw, dpidyfn, dpidyfs; // pressure gradient at cell faces
                double dpidxcc, dpidycc;                   // pressure gradient at cell center
                double uc0, vc0;                           // velocity at cell center
                double ufe, ufw, vfn, vfs;                 // velocity at cell faces
                int    ee1, ew1, en1, es1;                 // first level boundary indicator
                bE1(i, k, bf, ee1, ew1, en1, es1);

                uc0  =   u[i    ][k    ][_u];
                vc0  =   u[i    ][k    ][_v];
                ufe  = uFN[i    ][k    ][_u];
                ufw  = uFN[i - 1][k    ][_u];
                vfn  = uFN[i    ][k    ][_v];
                vfs  = uFN[i    ][k - 1][_v];

                pic0 =  pi[i    ][k    ];
                pie1 =  pi[i + 1][k    ];
                piw1 =  pi[i - 1][k    ];
                pin1 =  pi[i    ][k + 1];
                pis1 =  pi[i    ][k - 1];

                ufe  = (ee1 == NONE)? ufe  : bcf(BT[ee1][_u ], uc0 ,   dx, BC(ee1, _u , i + 0.5, k      , t));
                ufw  = (ew1 == NONE)? ufw  : bcf(BT[ew1][_u ], uc0 , - dx, BC(ew1, _u , i - 0.5, k      , t));
                vfn  = (en1 == NONE)? vfn  : bcf(BT[en1][_v ], vc0 ,   dy, BC(en1, _v , i      , k + 0.5, t));
                vfs  = (es1 == NONE)? vfs  : bcf(BT[es1][_v ], vc0 , - dy, BC(es1, _v , i      , k - 0.5, t));
                
                pie1 = (ee1 == NONE)? pie1 : bcc(BT[ee1][_pi], pic0,   dx, BC(ee1, _pi, i + 0.5, k      , t));
                piw1 = (ew1 == NONE)? piw1 : bcc(BT[ew1][_pi], pic0, - dx, BC(ew1, _pi, i - 0.5, k      , t));
                pin1 = (en1 == NONE)? pin1 : bcc(BT[en1][_pi], pic0,   dy, BC(en1, _pi, i      , k + 0.5, t));
                pis1 = (es1 == NONE)? pis1 : bcc(BT[es1][_pi], pic0, - dy, BC(es1, _pi, i      , k - 0.5, t));

                dpidxfe = (pie1 - pic0) / dx;
                dpidxfw = (pic0 - piw1) / dx;
                dpidyfn = (pin1 - pic0) / dy;
                dpidyfs = (pic0 - pis1) / dy;
                dpidxcc = 0.5 * (dpidxfe + dpidxfw);
                dpidycc = 0.5 * (dpidyfn + dpidyfs);

                uc0 = uc0 - dt * dpidxcc;
                vc0 = vc0 - dt * dpidycc;
                ufe = ufe - dt * dpidxfe;
                ufw = ufw - dt * dpidxfw;
                vfn = vfn - dt * dpidyfn;
                vfs = vfs - dt * dpidyfs;

                ufe  = (ee1 == NONE)? ufe  : bcf(BT[ee1][_u ], uc0 ,   dx, BC(ee1, _u , i + 0.5, k      , t));
                ufw  = (ew1 == NONE)? ufw  : bcf(BT[ew1][_u ], uc0 , - dx, BC(ew1, _u , i - 0.5, k      , t));
                vfn  = (en1 == NONE)? vfn  : bcf(BT[en1][_v ], vc0 ,   dy, BC(en1, _v , i      , k + 0.5, t));
                vfs  = (es1 == NONE)? vfs  : bcf(BT[es1][_v ], vc0 , - dy, BC(es1, _v , i      , k - 0.5, t));
                
                uF[i    ][k    ][_u] = ufe;
                uF[i - 1][k    ][_u] = ufw;
                uF[i    ][k    ][_v] = vfn;
                uF[i    ][k - 1][_v] = vfs;

                u[i][k][_u]          = uc0;
                u[i][k][_v]          = vc0;
            }
        }
    }
    updateP(p, pi, ff);
    return calcdiv(uF, u, ff, bf, div, t);
}

void ns2d(
    double   u[NNX][NNY][2], double  uA[NNX][NNY][2],
    double  uF[NNX][NNY][2], double uFN[NNX][NNY][2],
    double  pi[NNX][NNY]   , double   p[NNX][NNY]   , 
    int     ff[NNX][NNY]   , int     bf[NNX][NNY][2], 
    double div[NNX][NNY]
) {
    double t        = 0;
    double E        = 2E-4;
    double diver    = calcdiv(uF, u, ff, bf, div, t);
    double maxdiver = 0;
    int    it, it2;
    double r;
    printf("\rt = %5.5lf, poi = %5d, max div = %5.7lf\n", t, 0, diver);
    fflush(stdout);
    while(t < Tr) {
        diver    = DBL_MAX;
        maxdiver = 0;
        it2      = 0;
        prediction(u, uA, uF, uFN, p, ff, bf, t);
        while(diver > E) {
            initPi(pi);
            it    = poisson(pi, uF, u, ff, bf, r, t);
            diver = projection(u, uF, uFN, pi, p, ff, bf, div, t);
            if (diver > maxdiver) {
                maxdiver = diver;
            }
            printf("\rt = %5.5lf, it = %5d, poi = %5d, div = %5.7lf, res = %5.7lf, max div = %5.8lf", t, it2, it, diver, r, maxdiver);
            fflush(stdout);
            it2 += 1;
        }
        t += dt;
    }
    printf("\n");
}

void o2fo(double u[NNX][NNY][2], double p[NNX][NNY], double uF[NNX][NNY][2], double div[NNX][NNY]) {
    FILE  *fo;
    char   fname[128];
    double xpos, ypos;
    sprintf(fname, "Re%d@c.sor256.csv", int(Re));
    fo = fopen(fname, "w+t");

    if ( fo == NULL ) {
        printf("\nERROR when opening file\n");
    }
    else {
        fprintf(fo, "x,y,z,u,v,p,div\n");
        for (int k = 0; k < NNY; k ++) {
            for (int i = 0; i < NNX; i ++) {
                xpos = (i - WBOUND) * dx + 0.5 * dx;
                ypos = (k - SBOUND) * dy + 0.5 * dy;
                fprintf(fo, "%5.8lf,%5.8lf,%5.8lf,%5.8lf,%5.8lf,%5.8lf,%5.8lf\n", xpos, ypos, 0.0, u[i][k][_u], u[i][k][_v], p[i][k], div[i][k]);
            }
        }
        fclose(fo);
    }

    sprintf(fname, "Re%d@fx.sor256.csv", int(Re));
    fo = fopen(fname, "w+t");
    if ( fo == NULL ) {
        printf("\nERROR when opening file\n");
    }
    else {
        fprintf(fo, "x,y,z,uf\n");
        for (int k = 0; k < NNY; k ++) {
            for (int i = 0; i < NNX; i ++) {
                xpos = (i - WBOUND) * dx + dx;
                ypos = (k - SBOUND) * dy + 0.5 * dy;
                fprintf(fo, "%5.8lf,%5.8lf,%5.8lf,%5.8lf\n", xpos, ypos, 0.0, uF[i][k][_u]);
            }
        }
        fclose(fo);
    }

    sprintf(fname, "Re%d@fy.sor256.csv", int(Re));
    fo = fopen(fname, "w+t");
    if ( fo == NULL ) {
        printf("\nERROR when opening file\n");
    }
    else {
        fprintf(fo, "x,y,z,vf\n");
        for (int k = 0; k < NNY; k ++) {
            for (int i = 0; i < NNX; i ++) {
                xpos = (i - WBOUND) * dx + 0.5 * dx;
                ypos = (k - SBOUND) * dy + dy;
                fprintf(fo, "%5.8lf,%5.8lf,%5.8lf,%5.8lf\n", xpos, ypos, 0.0, uF[i][k][_v]);
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
    double U0 = ULID;
    double V0 = VLID;
    Re = strtod(argv[1], NULL);
    Tr = strtod(argv[2], NULL);
    dt = min(CFL / (U0 / dx + V0 / dy), DIFFC * Re * dx * dx * dy * dy / (2 * (dx * dx + dy * dy)));
    printf("Re = %lf, T = %lf, dt = %lf, dtau = %lf\n", Re, Tr, dt, dtau);

    clock_t start, finish;

    setupGeom(FF, BF);
    #pragma acc enter data copyin(U, UA, UF, UFN, Pi, P, FF, BF, BT, DIV)
    init(U, UF, P, Pi, FF, BF);

    start = clock();

    ns2d(U, UA, UF, UFN, Pi, P, FF, BF, DIV);

    finish = clock();
    #pragma acc exit data copyout(U, UF, P, DIV)
    o2fo(U, P, UF, DIV);

    printf("CPU time: %f sec\n", float(finish - start) / CLOCKS_PER_SEC);
    
    return 0;
}
