#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <algorithm>

using namespace std;

#define u_ 0
#define v_ 1

// cell numbers
// the data domain is φ[0:NN, 0:NN]
// the computational domain is φ[SWBOUND:NEBOUND, SWBOUND:NEBOUND]
const int NC = 256;
const int GHST = 2;
const int NN = NC + 2 * GHST;
const int SWBOUND = GHST;
const int NEBOUND = SWBOUND + NC;
const int NG = NC + 1;

// boundary values
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
const double FUBCe = 0;
const double FUBCw = 0;
const double FUBCn = 0;
const double FUBCs = 0;
const double FVBCe = 0;
const double FVBCw = 0;
const double FVBCn = 0;
const double FVBCs = 0;

// other constants
const double dd = 1.0 / NC;
const double B0 = 2.0 / (dd * dd) + 2.0 / (dd * dd);
const double dtau = 1.0 / B0;
const double CFL = 1.0;

// solver parameters
const int EPSILON = 1;
const double KAPPA = 1.0 / 3.0;
const double BETA = (3 - KAPPA) / (1 - KAPPA);

double Re;  // Reynolds number
double dt;  // timestep length
double ter; // simulation time span

// data vectors
double U[NN][NN][2], UN[NN][NN][2], UF[NN][NN][2], UFN[NN][NN][2], UG[NG][NG][2];
double P[NN][NN], PN[NN][NN], Pr[NN][NN], PrG[NG][NG];

// fetch boundary information
void bE(int i, int k, int &ee1, int &ew1, int &en1, int &es1, int &ee2, int &ew2, int &en2, int &es2) {
    ee1 = (i < NEBOUND - 1)? 1 : 0;
    ew1 = (i > SWBOUND    )? 1 : 0;
    en1 = (k < NEBOUND - 1)? 1 : 0;
    es1 = (k > SWBOUND    )? 1 : 0;

    ee2 = (i < NEBOUND - 2)? 1 : 0;
    ew2 = (i > SWBOUND + 1)? 1 : 0;
    en2 = (k < NEBOUND - 2)? 1 : 0;
    es2 = (k > SWBOUND + 1)? 1 : 0;
}

double calcdiv(double uF[NN][NN][2]) {
    double tdiver = 0;
    #pragma acc kernels loop independent reduction(+:tdiver) collapse(2) present(uF) copy(tdiver)
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double ufe, ufw, vfn, vfs; // velocity at each direction
            double dudx, dvdy;         // velocity gradient at cell center
            double diver;
            int    ee1, ew1, en1, es1;     // first level boundary indicator
            int    ee2, ew2, en2, es2;     // second level boundary indicator
            bE(i, k, ee1, ew1, en1, es1, ee2, ew2, en2, es2);
            ufe    = (ee1)? uF[i    ][k    ][u_] : UBCe;
            ufw    = (ew1)? uF[i - 1][k    ][u_] : UBCw;
            vfn    = (en1)? uF[i    ][k    ][v_] : VBCn;
            vfs    = (es1)? uF[i    ][k - 1][v_] : VBCs;
            dudx   = (ufe - ufw) / (dd);
            dvdy   = (vfn - vfs) / (dd);
            diver  = dudx + dvdy;
            tdiver = tdiver + diver * diver;
        }
    }
    return sqrt(tdiver);
}

// interpolate cell face values from cell center values
void calcUF(double u[NN][NN][2], double uF[NN][NN][2]) {
    #pragma acc kernels loop independent collapse(2) present(u, uF)
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double ue, uc, uw, vn, vc, vs; // velocity at each direction
            double ufe, ufw, vfn, vfs;     // velocity at cell faces
            int    ee1, ew1, en1, es1;     // first level boundary indicator
            int    ee2, ew2, en2, es2;     // second level boundary indicator
            bE(i, k, ee1, ew1, en1, es1, ee2, ew2, en2, es2);

            uc  = u[i    ][k    ][u_];
            vc  = u[i    ][k    ][v_];

            ue  = (ee1)? u[i + 1][k    ][u_] : (2 * UBCe - uc);
            uw  = (ew1)? u[i - 1][k    ][u_] : (2 * UBCw - uc);
            vn  = (en1)? u[i    ][k + 1][v_] : (2 * VBCn - vc);
            vs  = (es1)? u[i    ][k - 1][v_] : (2 * VBCs - vc);

            ufe = (ee1)? 0.5 * (uc + ue) : UBCe;
            ufw = (ew1)? 0.5 * (uw + uc) : UBCw;
            vfn = (en1)? 0.5 * (vc + vn) : VBCn;
            vfs = (es1)? 0.5 * (vs + vc) : VBCs;
            uF[i    ][k    ][u_] = ufe;
            uF[i - 1][k    ][u_] = ufw;
            uF[i    ][k    ][v_] = vfn;
            uF[i    ][k - 1][v_] = vfs;
        }
    }
}

// initialize pressure potentia
void initP(double p[NN][NN]) {
    #pragma acc kernels loop collapse(2) present(p)
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            p[i][k] = 0;
        }
    }
}

// initialize velocity and pressure
void init(double u[NN][NN][2], double uF[NN][NN][2], double p[NN][NN], double pr[NN][NN]) {
    #pragma acc kernels loop collapse(2) present(u)
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            u[i][k][u_] = 0;
            u[i][k][v_] = 0;
        }
    }
    #pragma acc kernels loop present(u)
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        u[i][NEBOUND - 1][u_] = UBCn;
    }
    calcUF(u, uF);

    initP(p);

    #pragma acc kernels loop collapse(2) present(pr)
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            pr[i][k] = 0;
        }
    }
}

// copy U to UN or UF to UFN
void cpUUF(double ut[NN][NN][2], double utN[NN][NN][2]) {
    #pragma acc kernels loop independent collapse(2) present(ut, utN)
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            utN[i][k][u_] = ut[i][k][u_];
            utN[i][k][v_] = ut[i][k][v_];
        }
    }
}
// copy P to PN
void cpP(double p[NN][NN], double pN[NN][NN]) {
    #pragma acc kernels loop independent collapse(2) present(p, pN)
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            pN[i][k] = p[i][k];
        }
    }
}

int sgn(double x) {
    if (x >= 0) {
        return 1;
    }
    else {
        return -1;
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
    return c - EPSILON * 0.25 * ((1 - KAPPA) * DP + (1 + KAPPA) * DM);
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
    return b + EPSILON * 0.25 * ((1 - KAPPA) * DM + (1 + KAPPA) * DP);
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

// u* = u + Δt * (- Convection + Viscoucity - grad(Pressure))
void fs1(double u[NN][NN][2], double uN[NN][NN][2], double uF[NN][NN][2], double pr[NN][NN]) {
    cpUUF(u, uN);
    #pragma acc kernels loop independent collapse(2) present(u, uN, uF, pr)
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double fue, fuw, fun, fus, fve, fvw, fvn, fvs;      // flux variables
            double uw2, uw1, us2, us1, uc0, un1, un2, ue1, ue2; // velocity at each direction
            double vw2, vw1, vs2, vs1, vc0, vn1, vn2, ve1, ve2; // velocity at each direction
            double ufe, ufw, vfn, vfs;                          // velocity at cell faces
            double duudx, dvudy, duvdx, dvvdy;                  // 1st derivatives
            double ddudxx, ddudyy, ddvdxx, ddvdyy;              // 2nd derivatives
            double pre, prw, prn, prs, prc;                     // pressure at each direction
            double dprdxfe, dprdxfw, dprdyfn, dprdyfs;          // pressure derivatives at cell faces
            double dprdxcc, dprdycc;                            // pressure derivatives at cell center
            int    ee1, ew1, en1, es1;                          // first level boundary indicator
            int    ee2, ew2, en2, es2;                          // second level boundary indicator
            bE(i, k, ee1, ew1, en1, es1, ee2, ew2, en2, es2);

            uc0     = uN[i    ][k    ][u_];
            vc0     = uN[i    ][k    ][v_];
            ue1     = (ee1)? uN[i + 1][k    ][u_] : (2 * UBCe - uc0);
            uw1     = (ew1)? uN[i - 1][k    ][u_] : (2 * UBCw - uc0);
            un1     = (en1)? uN[i    ][k + 1][u_] : (2 * UBCn - uc0);
            us1     = (es1)? uN[i    ][k - 1][u_] : (2 * UBCs - uc0);
            ve1     = (ee1)? uN[i + 1][k    ][v_] : (2 * VBCe - vc0);
            vw1     = (ew1)? uN[i - 1][k    ][v_] : (2 * VBCw - vc0);
            vn1     = (en1)? uN[i    ][k + 1][v_] : (2 * VBCn - vc0);
            vs1     = (es1)? uN[i    ][k - 1][v_] : (2 * VBCs - vc0);
            ue2     = (ee2)? uN[i + 2][k    ][u_] : (2 * UBCe - ue1);
            uw2     = (ew2)? uN[i - 2][k    ][u_] : (2 * UBCw - uw1);
            un2     = (en2)? uN[i    ][k + 2][u_] : (2 * UBCn - un1);
            us2     = (es2)? uN[i    ][k - 2][u_] : (2 * UBCs - us1);
            ve2     = (ee2)? uN[i + 2][k    ][v_] : (2 * VBCe - ve1);
            vw2     = (ew2)? uN[i - 2][k    ][v_] : (2 * VBCw - vw1);
            vn2     = (en2)? uN[i    ][k + 2][v_] : (2 * VBCn - vn1);
            vs2     = (es2)? uN[i    ][k - 2][v_] : (2 * VBCs - vs1);

            ufe     = (ee1)? uF[i    ][k    ][u_] : UBCe;
            ufw     = (ew1)? uF[i - 1][k    ][u_] : UBCw;
            vfn     = (en1)? uF[i    ][k    ][v_] : VBCn;
            vfs     = (es1)? uF[i    ][k - 1][v_] : VBCs;

            prc     = pr[i    ][k    ];
            pre     = pr[i + 1][k    ];
            prw     = pr[i - 1][k    ];
            prn     = pr[i    ][k + 1];
            prs     = pr[i    ][k - 1];

            fue     = (ee1)? flux(uw1, uc0, ue1, ue2, ufe) : FUBCe;
            fuw     = (ew1)? flux(uw2, uw1, uc0, ue1, ufw) : FUBCw;
            fun     = (en1)? flux(us1, uc0, un1, un2, vfn) : FUBCn;
            fus     = (es1)? flux(us2, us1, uc0, un1, vfs) : FUBCs;
            fve     = (ee1)? flux(vw1, vc0, ve1, ve2, ufe) : FVBCe;
            fvw     = (ew1)? flux(vw2, vw1, vc0, ve1, ufw) : FVBCw;
            fvn     = (en1)? flux(vs1, vc0, vn1, vn2, vfn) : FVBCn;
            fvs     = (es1)? flux(vs2, vs1, vc0, vn1, vfs) : FVBCs;

            duudx   = (fue - fuw) / dd;
            dvudy   = (fun - fus) / dd;
            duvdx   = (fve - fvw) / dd;
            dvvdy   = (fvn - fvs) / dd;
            ddudxx  = (uw1 - 2 * uc0 + ue1) / (dd * dd);
            ddudyy  = (us1 - 2 * uc0 + un1) / (dd * dd);
            ddvdxx  = (vw1 - 2 * vc0 + ve1) / (dd * dd);
            ddvdyy  = (vs1 - 2 * vc0 + vn1) / (dd * dd);

            dprdxfe = (ee1)? (pre - prc) / dd : DPBCe;
            dprdxfw = (ew1)? (prc - prw) / dd : DPBCw;
            dprdyfn = (en1)? (prn - prc) / dd : DPBCn;
            dprdyfs = (es1)? (prc - prs) / dd : DPBCs;
            dprdxcc = (dprdxfe + dprdxfw) * 0.5;
            dprdycc = (dprdyfn + dprdyfs) * 0.5;

            uc0     = uc0 + dt * (- duudx - dvudy - dprdxcc + (ddudxx + ddudyy) / Re);
            vc0     = vc0 + dt * (- duvdx - dvvdy - dprdycc + (ddvdxx + ddvdyy) / Re);

            u[i][k][u_] = uc0;
            u[i][k][v_] = vc0;
        }
    }
    calcUF(u, uF);
}

// div(grad(p)) = div(u*) / Δt
// p is the pressure potential, which is different from pressure Pr
int poisson(double p[NN][NN], double pN[NN][NN], double uF[NN][NN][2]) {
    double E  = 1E-5;
    int MAXIT = 100;
    bool con  = false;
    int it    = 0;
    double R  = 0;
    while (con == false) {
        R   = 0;
        con = true;
        cpP(p, pN);
        #pragma acc kernels loop independent reduction(+:R) collapse(2) present(p, pN, uF) copy(R)
        for (int i = SWBOUND; i < NEBOUND; i ++) {
            for (int k = SWBOUND; k < NEBOUND; k ++) {
                double ufe, ufw, vfn, vfs;             // velocity at cell faces
                double pe, pw, pc, pn, ps;             // pressure potential at each direction
                double dudx, dvdy;                     // velocity gradient at cell center
                double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // 1st derivatives at cell faces
                double ddpdxx, ddpdyy;                 // 2nd derivatives at cell center
                double psi, res;
                int    ee1, ew1, en1, es1;             // first level boundary indicator
                int    ee2, ew2, en2, es2;             // second level boundary indicator
                bE(i, k, ee1, ew1, en1, es1, ee2, ew2, en2, es2);

                ufe     = (ee1)? uF[i    ][k    ][u_] : UBCe;
                ufw     = (ew1)? uF[i - 1][k    ][u_] : UBCw;
                vfn     = (en1)? uF[i    ][k    ][v_] : VBCn;
                vfs     = (es1)? uF[i    ][k - 1][v_] : VBCs;

                pc      = pN[i    ][k    ];
                pe      = pN[i + 1][k    ];
                pw      = pN[i - 1][k    ];
                pn      = pN[i    ][k + 1];
                ps      = pN[i    ][k - 1];

                dudx    = (ufe - ufw) / (dd);
                dvdy    = (vfn - vfs) / (dd);

                dpdxfe  = (ee1)? (pe - pc) / dd : DPBCe;
                dpdxfw  = (ew1)? (pc - pw) / dd : DPBCw;
                dpdyfn  = (en1)? (pn - pc) / dd : DPBCn;
                dpdyfs  = (es1)? (pc - ps) / dd : DPBCs;
                ddpdxx  = (dpdxfe - dpdxfw) / dd;
                ddpdyy  = (dpdyfn - dpdyfs) / dd;
                psi     = (dudx + dvdy) / dt;
                res     = dtau * (ddpdxx + ddpdyy - psi);

                pc      = pc + res;
                R       = R + res * res;

                p[i][k] = pc;
            }
        }
        if (sqrt(R) > E) {
            con = false;
        }
        it += 1;
        if (it >= MAXIT) {
            break;
        }
    }
    return it;
}

// u = u* - Δt * grad(p)
double fs2(double u[NN][NN][2], double uF[NN][NN][2], double uFN[NN][NN][2], double p[NN][NN]) {
    cpUUF(uF, uFN);
    #pragma acc kernels loop independent collapse(2) present(u, uF, uFN, p)
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double pe, pw, pc, pn, ps;             // pressure potential at each direction
            double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // pressure gradient at cell faces
            double dpdxcc, dpdycc;                 // pressure gradient at cell center
            double ucc, vcc;                       // velocity at cell center
            double ufe, ufw, vfn, vfs;             // velocity at cell faces
            int    ee1, ew1, en1, es1;             // first level boundary indicator
            int    ee2, ew2, en2, es2;             // second level boundary indicator
            bE(i, k, ee1, ew1, en1, es1, ee2, ew2, en2, es2);

            ucc    = u[i    ][k    ][u_];
            vcc    = u[i    ][k    ][v_];
            ufe    = (ee1)? uFN[i    ][k    ][u_] : UBCe;
            ufw    = (ew1)? uFN[i - 1][k    ][u_] : UBCw;
            vfn    = (en1)? uFN[i    ][k    ][v_] : VBCn;
            vfs    = (es1)? uFN[i    ][k - 1][v_] : VBCs;

            pc     = p[i    ][k    ];
            pe     = p[i + 1][k    ];
            pw     = p[i - 1][k    ];
            pn     = p[i    ][k + 1];
            ps     = p[i    ][k - 1];

            dpdxfe  = (ee1)? (pe - pc) / dd : DPBCe;
            dpdxfw  = (ew1)? (pc - pw) / dd : DPBCw;
            dpdyfn  = (en1)? (pn - pc) / dd : DPBCn;
            dpdyfs  = (es1)? (pc - ps) / dd : DPBCs;
            dpdxcc = 0.5 * (dpdxfe + dpdxfw);
            dpdycc = 0.5 * (dpdyfn + dpdyfs);

            ufe    = (ee1)? (ufe - dt * dpdxfe) : UBCe;
            ufw    = (ew1)? (ufw - dt * dpdxfw) : UBCw;
            vfn    = (en1)? (vfn - dt * dpdyfn) : VBCn;
            vfs    = (es1)? (vfs - dt * dpdyfs) : VBCs;
            ucc    = ucc - dt * dpdxcc;
            vcc    = vcc - dt * dpdycc;

            uF[i    ][k    ][u_] = ufe;
            uF[i - 1][k    ][u_] = ufw;
            uF[i    ][k    ][v_] = vfn;
            uF[i    ][k - 1][v_] = vfs;
            u[i][k][u_]          = ucc;
            u[i][k][v_]          = vcc;
        }
    }
    return calcdiv(uF);
}

// Pr = Pr + p
void updatePr(double p[NN][NN], double pr[NN][NN]) {
    #pragma acc kernels loop independent collapse(2) present(p, pr)
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            pr[i][k] = pr[i][k] + p[i][k];
        }
    }
}

void ns2d(double u[NN][NN][2], double uN[NN][NN][2], double uF[NN][NN][2], double uFN[NN][NN][2], double p[NN][NN], double pN[NN][NN], double pr[NN][NN]) {
    double t        = 0;
    double E        = 5E-3;
    double diver    = calcdiv(uF);
    double maxdiver = 0;
    int it;
    printf("\rt = %5.8lf, poi = %5d, max div = %5.8lf", t, 0, diver);
    fflush(stdout);
    while (t < ter) {
        diver    = DBL_MAX;
        maxdiver = 0;
        fs1(u, uN, uF, pr);
        while (diver > E) {
            initP(p);
            it = poisson(p, pN, uF);
            diver = fs2(u, uF, uFN, p);
            if (diver > maxdiver) {
                maxdiver = diver;
            }
            printf("\rt = %5.8lf, poi = %5d, div = %5.8lf, max div = %5.8lf", t, it, diver, maxdiver);
            fflush(stdout);
            updatePr(p, pr);
        }
        t += dt;
    }
    printf("\n");
}

void calcgrid(double u[NN][NN][2], double uG[NG][NG][2], double pr[NN][NN], double prG[NG][NG]) {
    int OFST = GHST - 1;
    double une, unw, use, usw, vne, vnw, vse, vsw; // velocity at each direction
    double pne, pnw, pse, psw;                     // pressure at each direction
    double ugc, vgc;                               // velocity at grid point
    double pgc;                                    // pressure at grid point
    for (int i = 0; i < NG; i ++) {
        for (int k = 0; k < NG; k ++) {
            une = u[i + OFST + 1][k + OFST + 1][u_];
            unw = u[i + OFST    ][k + OFST + 1][u_];
            use = u[i + OFST + 1][k + OFST    ][u_];
            usw = u[i + OFST    ][k + OFST    ][u_];
            vne = u[i + OFST + 1][k + OFST + 1][v_];
            vnw = u[i + OFST    ][k + OFST + 1][v_];
            vse = u[i + OFST + 1][k + OFST    ][v_];
            vsw = u[i + OFST    ][k + OFST    ][v_];
            ugc = 0.25 * (une + unw + use + usw);
            vgc = 0.25 * (vne + vnw + vse + vsw);

            uG[i][k][u_] = ugc;
            uG[i][k][v_] = vgc;
        }
    }
    for (int i = 0; i < NG; i ++) {
        uG[i][0     ][u_] = UBCs;
        uG[i][0     ][v_] = VBCs;
        uG[i][NG - 1][u_] = UBCn;
        uG[i][NG - 1][v_] = VBCn;
    }
    for (int k = 0; k < NG; k ++) {
        uG[0     ][k][u_] = UBCw;
        uG[0     ][k][v_] = VBCw;
        uG[NG - 1][k][u_] = UBCe;
        uG[NG - 1][k][v_] = VBCe;
    }
    uG[0     ][0     ][u_] = 0.5 * (UBCw + UBCs);
    uG[0     ][NG - 1][u_] = 0.5 * (UBCw + UBCn);
    uG[NG - 1][0     ][u_] = 0.5 * (UBCe + UBCs);
    uG[NG - 1][NG - 1][u_] = 0.5 * (UBCe + UBCn);
    uG[0     ][0     ][v_] = 0.5 * (VBCw + VBCs);
    uG[0     ][NG - 1][v_] = 0.5 * (VBCw + VBCn);
    uG[NG - 1][0     ][v_] = 0.5 * (VBCe + VBCs);
    uG[NG - 1][NG - 1][v_] = 0.5 * (VBCe + VBCn);
    for (int i = 0; i < NG; i ++) {
        for (int k = 0; k < NG; k ++) {
            pne = pr[i + OFST + 1][k + OFST + 1];
            pnw = pr[i + OFST    ][k + OFST + 1];
            pse = pr[i + OFST + 1][k + OFST    ];
            psw = pr[i + OFST    ][k + OFST    ];
            pgc = 0.25 * (pne + pnw + pse + psw);

            prG[i][k] = pgc;
        }
    }
    for (int i = 0; i < NG; i ++) {
        int k = 0;
        pne = pr[i + OFST + 1][k + OFST + 1];
        pnw = pr[i + OFST    ][k + OFST + 1];
        prG[i][k] = 0.5 * (pne + pnw);
        k = NG - 1;
        pse = pr[i + OFST + 1][k + OFST    ];
        psw = pr[i + OFST    ][k + OFST    ];
        prG[i][k] = 0.5 * (pse + psw);
    }
    for (int k = 0; k < NG; k ++) {
        int i = 0;
        pne = pr[i + OFST + 1][k + OFST + 1];
        pse = pr[i + OFST + 1][k + OFST    ];
        prG[i][k] = 0.5 * (pne + pse);
        i = NG - 1;
        pnw = pr[i + OFST    ][k + OFST + 1];
        psw = pr[i + OFST    ][k + OFST    ];
        prG[i][k] = 0.5 * (pnw + psw);
    }
    prG[0     ][0     ] = pr[OFST + 1     ][OFST + 1     ];
    prG[0     ][NG - 1] = pr[OFST + 1     ][NG - 1 + OFST];
    prG[NG - 1][0     ] = pr[NG - 1 + OFST][OFST + 1     ];
    prG[NG - 1][NG - 1] = pr[NG - 1 + OFST][NG - 1 + OFST];
}

void o2f(double uG[NG][NG][2], double pG[NG][NG]) {
    FILE *fo;
    char fname[128];
    double xpos, ypos;
    sprintf(fname, "UVP_Re%d.ns2dccg256.csv", int(Re), int(ter));
    fo = fopen(fname, "w+t");

    if ( fo == NULL ) {
        printf("\nERROR when opening file\n");
    }
    else {
        fprintf(fo, "x,y,z,u,v,p\n");
        for (int k = 0; k < NG; k ++) {
            for (int i = 0; i < NG; i ++) {
                xpos = i * dd;
                ypos = k * dd;
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
    ter = strtod(argv[2], NULL);
    dt = min(CFL * dd * 0.5, 0.5 * Re * dd * dd * dd * dd / (2 * (dd * dd + dd * dd)));
    printf("Re = %lf, T = %lf, dt = %lf, dtau = %lf\n", Re, ter, dt, dtau);

    #pragma acc enter data copyin(U, UN, UF, UFN, P, PN, Pr)
    init(U, UF, P, Pr);
    ns2d(U, UN, UF, UFN, P, PN, Pr);
    #pragma acc exit data copyout(U, UF, Pr)
    calcgrid(U, UG, Pr, PrG);
    o2f(UG, PrG);

    return 0;
}
