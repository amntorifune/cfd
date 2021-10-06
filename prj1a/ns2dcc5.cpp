#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <algorithm>

using namespace std;

#define u_ 0
#define v_ 1

// cell numbers
const int NC = 128;
const int NN = NC + 4;
const int LLBOUND = 2;
const int RUBOUND = LLBOUND + NC;
const int NG = NC + 1;

// bouncary values
const double UBCe  = 0;
const double UBCw  = 0;
const double UBCn  = 1;
const double UBCs  = 0;
const double VBCe  = 0;
const double VBCw  = 0;
const double VBCn  = 0;
const double VBCs  = 0;
const double DPBCe = 0;
const double DPBCw = 0;
const double DPBCn = 0;
const double DPBCs = 0;

// other constants
const double dd = 1.0 / NC;
const double B0 = 2.0 / (dd * dd) + 2.0 / (dd * dd);
const double dtau = 1.0 / B0;

// solver parameters
const int EPSILON = 1;
const double KAPPA = 1.0 / 3.0;
const double BETA = (3 - KAPPA) / (1 - KAPPA);

// Reynolds number
double Re;
// timestep length
double dt;

// data vectors
double U[NN][NN][2], UN[NN][NN][2], UF[NN][NN][2], UG[NG][NG][2];
double P[NN][NN], PN[NN][NN], PG[NG][NG];

void cpU(double u[NN][NN][2], double uN[NN][NN][2]) {
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            uN[i][k][u_] = u[i][k][u_];
            uN[i][k][v_] = u[i][k][v_];
        }
    }
}
void cpP(double p[NN][NN], double pN[NN][NN]) {
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            pN[i][k] = p[i][k];
        }
    }
}

void be(int i, int k, int &ee, int &ew, int &en, int &es) {
    ee = (i < RUBOUND - 1)? 1 : 0;
    ew = (i > LLBOUND    )? 1 : 0;
    en = (k < RUBOUND - 1)? 1 : 0;
    es = (k > LLBOUND    )? 1 : 0;
}

double flux(double phiL, double phiR, double uF) {
    return uF * (phiL + phiR) * 0.5;
}

void fs1(double u[NN][NN][2], double uN[NN][NN][2], double uF[NN][NN][2]) {
    cpU(u, uN);
    double fe, fw, fn, fs;                         // flux variables
    double ue, uw, uc, un, us, ve, vw, vc, vn, vs; // velocity at each direction
    double ufe, ufw, vfn, vfs;                     // velocity at cell faces
    int ee, ew, en, es;                            // boundary indicator
    double duudx, dvudy, duvdx, dvvdy;             // 1st derivatives
    double ddudxx, ddudyy, ddvdxx, ddvdyy;         // 2nd derivatives
    for (int i = LLBOUND; i < RUBOUND; i ++) {
        for (int k = LLBOUND; k < RUBOUND; k ++) {
            be(i, k, ee, ew, en, es);
            uc = uN[i][k][u_];
            vc = uN[i][k][v_];
            ue =  ee * uN[i + 1][k    ][u_] + (1 - ee) * (2 * UBCe - uc);
            uw =  ew * uN[i - 1][k    ][u_] + (1 - ew) * (2 * UBCw - uc);
            un =  en * uN[i    ][k + 1][u_] + (1 - en) * (2 * UBCn - uc);
            us =  es * uN[i    ][k - 1][u_] + (1 - es) * (2 * UBCs - uc);
            ve =  ee * uN[i + 1][k    ][v_] + (1 - ee) * (2 * VBCe - vc);
            vw =  ew * uN[i - 1][k    ][v_] + (1 - ew) * (2 * VBCw - vc);
            vn =  en * uN[i    ][k + 1][v_] + (1 - en) * (2 * VBCn - vc);
            vs =  es * uN[i    ][k - 1][v_] + (1 - es) * (2 * VBCs - vc);
            ufe = ee * uF[i    ][k    ][u_] + (1 - ee) * UBCe;
            ufw = ew * uF[i - 1][k    ][u_] + (1 - ew) * UBCw;
            vfn = en * uF[i    ][k    ][v_] + (1 - en) * VBCn;
            vfs = es * uF[i    ][k - 1][v_] + (1 - es) * VBCs;

            // duu/dx = (fe~ - fw~)/dx
            fe = flux(uc, ue, ufe);
            fw = flux(uw, uc, ufw);
            duudx = (fe - fw) / dd;
            // dvu/dx = (fn~ = fs~)/dy
            fn = flux(uc, un, vfn);
            fs = flux(us, uc, vfs);
            dvudy = (fn - fs) / dd;
            // viscoucity
            ddudxx = (uw - 2 * uc + ue) / (dd * dd);
            ddudyy = (un - 2 * uc + us) / (dd * dd);
            u[i][k][u_] = uc + dt * (- duudx - dvudy + (ddudxx + ddudyy) / Re);

            // duv/dx = (fe~ - fw~)/dx
            fe = flux(vc, ve, ufe);
            fw = flux(vw, vc, ufw);
            duvdx = (fe - fw) / dd;
            // dvv/dy = (fn~ - fs~)/dy
            fn = flux(vc, vn, vfn);
            fs = flux(vs, vc, vfs);
            dvvdy = (fn - fs) / dd;
            // viscoucity
            ddvdxx = (vw - 2 * vc + ve) / (dd * dd);
            ddvdyy = (vn - 2 * vc + vs) / (dd * dd);
            u[i][k][v_] = vc + dt * (- duvdx - dvvdy + (ddvdxx + ddvdyy) / Re);
        }
    }
    // update velocity at cell faces: uf = (u + u_neighbor)/2
    for (int i = LLBOUND; i < RUBOUND; i ++) {
        for (int k = LLBOUND; k < RUBOUND; k ++) {
            be(i, k, ee, ew, en, es);
            ufe = ee * 0.5 * (u[i    ][k    ][u_] + u[i + 1][k    ][u_]) + (1 - ee) * UBCe;
            ufw = ew * 0.5 * (u[i - 1][k    ][u_] + u[i    ][k    ][u_]) + (1 - ew) * UBCw;
            vfn = en * 0.5 * (u[i    ][k    ][v_] + v[i    ][k + 1][v_]) + (1 - en) * VBCn;
            vfs = es * 0.5 * (u[i    ][k - 1][v_] + v[i    ][k    ][v_]) + (1 - es) * VBCs;
            uF[i    ][k    ][u_] = ufe;
            uF[i - 1][k    ][u_] = ufw;
            uF[i    ][k    ][v_] = vfn;
            uF[i    ][k - 1][v_] = vfs;
        }
    }
}

int poisson(double p[NN][NN], double pN[NN][NN], double uF[NN][NN][2]) {
    double E = 0.02;
    int MAXIT = 10000;
    bool con = false;
    int it = 0;
    int ee, ew, en, es;                    // boundary indicator
    double ufe, ufw, vfn, vfs;             // velocity at cell faces
    double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // dp/dx at cell faces
    double ddpdxx, ddpdyy;                 // 2nd derivatives at cell center
    double dudx, dvdy;                     // velocity gradient at cell center
    double res;
    while (con == false) {
        con = true;
        cpP(p, pN);
        for (int i = LLBOUND; i < RUBOUND; i ++) {
            for (int k = LLBOUND; k < RUBOUND; k ++) {
                be(i, k, ee, ew, en, es);
                dpdxfe = ee * (pN[i + 1][k    ] - pN[i    ][k    ]) / dd + (1 - ee) * DPBCe;
                dpdxfw = ew * (pN[i    ][k    ] - pN[i - 1][k    ]) / dd + (1 - ew) * DPBCw;
                dpdyfn = en * (pN[i    ][k + 1] - pN[i    ][k    ]) / dd + (1 - en) * DPBCn;
                dpdyfs = es * (pN[i    ][k    ] - pN[i    ][k - 1]) / dd + (1 - es) * DPBCs;
                ddpdxx = (dpdxfe - dpdxfw) / dd;
                ddpdyy = (dpdyfn - dpdyfs) / dd;
                ufe = ee * uF[i    ][k    ][u_] + (1 - ee) * UBCe;
                ufw = ew * uF[i - 1][k    ][u_] + (1 - ew) * UBCw;
                vfn = en * uF[i    ][k    ][v_] + (1 - en) * VBCn;
                vfs = es * uF[i    ][k - 1][v_] + (1 - es) * VBCs;
                dudx = (ufe - ufw) / dd;
                dvdy = (vfn - vfs) / dd;
                res = dtau * (ddpdxx + ddpdyy - (dudx + dvdy) / dt);
                p[i][k] = pN[i][k] + res;
                if (abs(res) > abs(p[i][k]) * E) {
                    con = false;
                }
            }
        }
        it += 1;
        if (it >= MAXIT) {
            break;
        }
    }
    return it;
}

void fs2(double u[NN][NN][2], double uF[NN][NN][2], double p[NN][NN]) {
    int ee, ew, en, es;                    // boundary indicator
    double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // pressure gradient at cell faces
    double dpdxcc, dpdycc;                 // pressure gradient at cell center
    double ufe, ufw, vfn, vfs;             // velocity at cell faces
    double ucc, vcc;                       // velocity at cell center
    for (int i = LLBOUND; i < RUBOUND; i ++) {
        for (int k = LLBOUND; k < RUBOUND; k ++) {
            be(i, k, ee, ew, en, es);
            ucc = u[i][k][u_];
            vcc = u[i][k][v_];
            ufe = uF[i    ][k    ][u_];
            ufw = uF[i - 1][k    ][u_];
            vfn = uF[i    ][k    ][v_];
            vfs = uF[i    ][k - 1][v_];
            dpdxfe = ee * (p[i + 1][k    ] - p[i    ][k    ]) / dd + (1 - ee) * DPBCe;
            dpdxfw = ew * (p[i    ][k    ] - p[i - 1][k    ]) / dd + (1 - ew) * DPBCw;
            dpdyfn = en * (p[i    ][k + 1] - p[i    ][k    ]) / dd + (1 - en) * DPBCn;
            dpdyfs = es * (p[i    ][k    ] - p[i    ][k - 1]) / dd + (1 - es) * DPBCs;
            dpdxcc = 0.5 * (dpdxfe + dpdxfw);
            dpdycc = 0.5 * (dpdyfn + dpdyfs);

            // update velocity at cell center
            u[i][k][u_] = ucc - dt * dpdxcc;
            u[i][k][v_] = vcc - dt * dpdycc;
            
            // update velocity at cell faces
        }
    }
}