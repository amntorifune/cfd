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
const int NC = 128;
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
double P[NN][NN], PN[NN][NN], PG[NG][NG];

double calcdiv(double uF[NN][NN][2]) {
    double tdiver = 0;
    #pragma omp parallel for reduction(+:tdiver)
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double ufe, ufw, vfn, vfs; // velocity at each direction
            double dudx, dvdy;         // velocity gradient at cell center
            double diver;
            ufe    = uF[i    ][k    ][u_];
            ufw    = uF[i - 1][k    ][u_];
            vfn    = uF[i    ][k    ][v_];
            vfs    = uF[i    ][k - 1][v_];
            dudx   = (ufe - ufw) / (dd);
            dvdy   = (vfn - vfs) / (dd);
            diver  = dudx + dvdy;
            tdiver = tdiver + diver * diver;
        }
    }
    return sqrt(tdiver);
}

// apply boundary conditions to velocity at cell center
void applyBCU(double u[NN][NN][2]) {
    #pragma omp parallel for
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        u[i][SWBOUND - 1][u_] = 2 * UBCs - u[i][SWBOUND    ][u_];
        u[i][NEBOUND    ][u_] = 2 * UBCn - u[i][NEBOUND - 1][u_];
        u[i][SWBOUND - 1][v_] = 2 * VBCs - u[i][SWBOUND    ][v_];
        u[i][NEBOUND    ][v_] = 2 * VBCn - u[i][NEBOUND - 1][v_];
    }
    #pragma omp parallel for
    for (int k = SWBOUND; k < NEBOUND; k ++) {
        u[SWBOUND - 1][k][u_] = 2 * UBCw - u[SWBOUND    ][k][u_];
        u[NEBOUND    ][k][u_] = 2 * UBCe - u[NEBOUND - 1][k][u_];
        u[SWBOUND - 1][k][v_] = 2 * VBCw - u[SWBOUND    ][k][v_];
        u[NEBOUND    ][k][v_] = 2 * VBCe - u[NEBOUND - 1][k][v_];
    }
}
// apply boundary conditions to velocity at cell faces
void applyBCUF(double uF[NN][NN][2]) {
    #pragma omp parallel for
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        uF[i][SWBOUND - 1][v_] = VBCs;
        uF[i][NEBOUND - 1][v_] = VBCn;
    }
    #pragma omp parallel for
    for (int k = SWBOUND; k < NEBOUND; k ++) {
        uF[SWBOUND - 1][k][u_] = UBCw;
        uF[NEBOUND - 1][k][u_] = UBCe;
    }
}
// apply boundary conditions to pressure potential
void applyBCP(double p[NN][NN]) {
    #pragma omp parallel for
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        p[i][NEBOUND    ] = p[i][NEBOUND - 1] + DPBCn * dd;
        p[i][SWBOUND - 1] = p[i][SWBOUND    ] - DPBCs * dd;
    }
    #pragma omp parallel for
    for (int k = SWBOUND; k < NEBOUND; k ++) {
        p[NEBOUND    ][k] = p[NEBOUND - 1][k] + DPBCe * dd;
        p[SWBOUND - 1][k] = p[SWBOUND    ][k] - DPBCw * dd;
    }
}

// interpolate cell face values from cell center values
void calcUF(double u[NN][NN][2], double uF[NN][NN][2]) {
    #pragma omp parallel for
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double ue, uc, uw, vn, vc, vs; // velocity at each direction
            double ufe, ufw, vfn, vfs;     // velocity at cell faces
            ue  = u[i + 1][k    ][u_];
            uc  = u[i    ][k    ][u_];
            uw  = u[i - 1][k    ][u_];
            vn  = u[i    ][k + 1][v_];
            vc  = u[i    ][k    ][v_];
            vs  = u[i    ][k - 1][v_];
            ufe = 0.5 * (uc + ue);
            ufw = 0.5 * (uw + uc);
            vfn = 0.5 * (vc + vn);
            vfs = 0.5 * (vs + vc);
            uF[i    ][k    ][u_] = ufe;
            uF[i - 1][k    ][u_] = ufw;
            uF[i    ][k    ][v_] = vfn;
            uF[i    ][k - 1][v_] = vfs;
        }
    }
}

// fetch boundary information
void bE(int i, int k, int &ee, int &ew, int &en, int &es) {
    ee = (i < NEBOUND - 1)? 1 : 0;
    ew = (i > SWBOUND    )? 1 : 0;
    en = (k < NEBOUND - 1)? 1 : 0;
    es = (k > SWBOUND    )? 1 : 0;
}

// initialize velocity and pressure
void init(double u[NN][NN][2], double uF[NN][NN][2], double p[NN][NN]) {
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            u[i][k][u_] = 0;
            u[i][k][v_] = 0;
        }
    }
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        u[i][NEBOUND - 1][u_] = UBCn;
    }
    applyBCU(u);
    calcUF(u, uF);
    applyBCUF(uF);

    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            p[i][k] = 0;
        }
    }
    applyBCP(p);
}

// copy U to UN or UF to UFN
void cpUUF(double u[NN][NN][2], double uN[NN][NN][2]) {
    #pragma omp parallel for
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            uN[i][k][u_] = u[i][k][u_];
            uN[i][k][v_] = u[i][k][v_];
        }
    }
}
// copy P to PN
void cpP(double p[NN][NN], double pN[NN][NN]) {
    #pragma omp parallel for
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

void fs1(double u[NN][NN][2], double uN[NN][NN][2], double uF[NN][NN][2]) {
    cpUUF(u, uN);
    #pragma omp parallel for
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double fue, fuw, fun, fus, fve, fvw, fvn, fvs;      // flux variables
            double uw2, uw1, us2, us1, uc0, un1, un2, ue1, ue2; // velocity at each direction
            double vw2, vw1, vs2, vs1, vc0, vn1, vn2, ve1, ve2; // velocity at each direction
            double ufe, ufw, vfn, vfs;                          // velocity at cell faces
            double duudx, dvudy, duvdx, dvvdy;                  // 1st derivatives
            double ddudxx, ddudyy, ddvdxx, ddvdyy;              // 2nd derivatives
            int ee, ew, en, es;                              // boundary indicators
            bE(i, k, ee, ew, en, es);
            uw2    = uN[i - 2][k    ][u_];
            uw1    = uN[i - 1][k    ][u_];
            us2    = uN[i    ][k - 2][u_];
            us1    = uN[i    ][k - 1][u_];
            uc0    = uN[i    ][k    ][u_];
            un1    = uN[i    ][k + 1][u_];
            un2    = uN[i    ][k + 2][u_];
            ue1    = uN[i + 1][k    ][u_];
            ue2    = uN[i + 2][k    ][u_];
            vw2    = uN[i - 2][k    ][v_];
            vw1    = uN[i - 1][k    ][v_];
            vs2    = uN[i    ][k - 2][v_];
            vs1    = uN[i    ][k - 1][v_];
            vc0    = uN[i    ][k    ][v_];
            vn1    = uN[i    ][k + 1][v_];
            vn2    = uN[i    ][k + 2][v_];
            ve1    = uN[i + 1][k    ][v_];
            ve2    = uN[i + 2][k    ][v_];

            ufe    = uF[i    ][k    ][u_];
            ufw    = uF[i - 1][k    ][u_];
            vfn    = uF[i    ][k    ][v_];
            vfs    = uF[i    ][k - 1][v_];

            fue    = (ee)? flux(uw1, uc0, ue1, ue2, ufe) : FUBCe;
            fuw    = (ew)? flux(uw2, uw1, uc0, ue1, ufw) : FUBCw;
            fun    = (en)? flux(us1, uc0, un1, un2, vfn) : FUBCn;
            fus    = (es)? flux(us2, us1, uc0, un1, vfs) : FUBCs;
            fve    = (ee)? flux(vw1, vc0, ve1, ve2, ufe) : FVBCe;
            fvw    = (ew)? flux(vw2, vw1, vc0, ve1, ufw) : FVBCw;
            fvn    = (en)? flux(vs1, vc0, vn1, vn2, vfn) : FVBCn;
            fvs    = (es)? flux(vs2, vs1, vc0, vn1, vfs) : FVBCs;

            duudx  = (fue - fuw) / dd;
            dvudy  = (fun - fus) / dd;
            duvdx  = (fve - fvw) / dd;
            dvvdy  = (fvn - fvs) / dd;
            ddudxx = (uw1 - 2 * uc0 + ue1) / (dd * dd);
            ddudyy = (us1 - 2 * uc0 + un1) / (dd * dd);
            ddvdxx = (vw1 - 2 * vc0 + ve1) / (dd * dd);
            ddvdyy = (vs1 - 2 * vc0 + vn1) / (dd * dd);

            uc0    = uc0 + dt * (- duudx - dvudy + (ddudxx + ddudyy) / Re);
            vc0    = vc0 + dt * (- duvdx - dvvdy + (ddvdxx + ddvdyy) / Re);

            u[i][k][u_] = uc0;
            u[i][k][v_] = vc0;
        }
    }
    calcUF(u, uF);
}

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
        #pragma omp parallel for reduction(+:R)
        for (int i = SWBOUND; i < NEBOUND; i ++) {
            for (int k = SWBOUND; k < NEBOUND; k ++) {
                double ufe, ufw, vfn, vfs; // velocity at cell faces
                double pe, pw, pc, pn, ps; // pressure potential at each direction
                double dudx, dvdy;         // velocity gradient at cell center
                double ddpdxx, ddpdyy;     // 2nd derivatives at cell center
                double psi, res;
                ufe     = uF[i    ][k    ][u_];
                ufw     = uF[i - 1][k    ][u_];
                vfn     = uF[i    ][k    ][v_];
                vfs     = uF[i    ][k - 1][v_];
                pe      = pN[i + 1][k    ];
                pw      = pN[i - 1][k    ];
                pc      = pN[i    ][k    ];
                pn      = pN[i    ][k + 1];
                ps      = pN[i    ][k - 1];

                dudx    = (ufe - ufw) / (dd);
                dvdy    = (vfn - vfs) / (dd);
                ddpdxx  = (pw - 2 * pc + pe) / (dd * dd);
                ddpdyy  = (ps - 2 * pc + pn) / (dd * dd);
                psi     = (dudx + dvdy) / dt;
                res     = dtau * (ddpdxx + ddpdyy - psi);

                pc      = pc + res;
                R       = R + res * res;

                p[i][k] = pc;
            }
        }
        applyBCP(p);
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

double fs2(double u[NN][NN][2], double uF[NN][NN][2], double uFN[NN][NN][2], double p[NN][NN]) {
    cpUUF(uF, uFN);
    #pragma omp parallel for
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            double pe, pw, pc, pn, ps;             // pressure potential at each direction
            double dpdxfe, dpdxfw, dpdyfn, dpdyfs; // pressure gradient at cell faces
            double dpdxcc, dpdycc;                 // pressure gradient at cell center
            double ucc, vcc;                       // velocity at cell center
            double ufe, ufw, vfn, vfs;             // velocity at cell faces
            ucc    =   u[i    ][k    ][u_];
            vcc    =   u[i    ][k    ][v_];
            ufe    = uFN[i    ][k    ][u_];
            ufw    = uFN[i - 1][k    ][u_];
            vfn    = uFN[i    ][k    ][v_];
            vfs    = uFN[i    ][k - 1][v_];
            pe     =   p[i + 1][k    ];
            pw     =   p[i - 1][k    ];
            pc     =   p[i    ][k    ];
            pn     =   p[i    ][k + 1];
            ps     =   p[i    ][k - 1];

            dpdxfe = (pe - pc) / dd;
            dpdxfw = (pc - pw) / dd;
            dpdyfn = (pn - pc) / dd;
            dpdyfs = (pc - ps) / dd;
            dpdxcc = 0.5 * (dpdxfe + dpdxfw);
            dpdycc = 0.5 * (dpdyfn + dpdyfs);

            ufe    = ufe - dt * dpdxfe;
            ufw    = ufw - dt * dpdxfw;
            vfn    = vfn - dt * dpdyfn;
            vfs    = vfs - dt * dpdyfs;
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
    applyBCU(u);
    applyBCUF(uF);
    return calcdiv(uF);
}

void ns2d(double u[NN][NN][2], double uN[NN][NN][2], double uF[NN][NN][2], double uFN[NN][NN][2], double p[NN][NN], double pN[NN][NN]) {
    double t        = 0;
    double E        = 1E-2;
    double diver    = calcdiv(uF);
    double maxdiver = 0;
    int it;
    printf("\rt = %5.8lf, poi = %5d, max div = %5.8lf", t, 0, diver);
    fflush(stdout);
    while (t < ter) {
        diver    = DBL_MAX;
        maxdiver = 0;
        fs1(u, uN, uF);
        while (diver > E) {
            it = poisson(p, pN, uF);
            diver = fs2(u, uF, uFN, p);
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

void calcgrid(double u[NN][NN][2], double uG[NG][NG][2], double p[NN][NN], double pG[NG][NG]) {
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
            pne = p[i + OFST + 1][k + OFST + 1];
            pnw = p[i + OFST    ][k + OFST + 1];
            pse = p[i + OFST + 1][k + OFST    ];
            psw = p[i + OFST    ][k + OFST    ];
            pgc = 0.25 * (pne + pnw + pse + psw);

            pG[i][k] = pgc;
        }
    }
    pG[0     ][0     ] = p[OFST + 1     ][OFST + 1     ];
    pG[0     ][NG - 1] = p[OFST + 1     ][NG - 1 + OFST];
    pG[NG - 1][0     ] = p[NG - 1 + OFST][OFST + 1     ];
    pG[NG - 1][NG - 1] = p[NG - 1 + OFST][NG - 1 + OFST];
}

void o2f(double uG[NG][NG][2], double pG[NG][NG]) {
    FILE *fo;
    char fname[128];
    double xpos, ypos;
    sprintf(fname, "UVP_Re%d_t%d.ns2dcc8.4.csv", int(Re), int(ter));
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

    init(U, UF, P);
    ns2d(U, UN, UF, UFN, P, PN);
    calcgrid(U, UG, P, PG);
    o2f(UG, PG);

    return 0;
}
