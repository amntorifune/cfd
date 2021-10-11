#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <algorithm>

using namespace std;

#define u_ 0
#define v_ 1

// cell numbers
// the data domain is φ[0:NN, 0:NN]
// the computational domain is φ[SWBOUND:NEBOUND, SWBOUND:NEBOUND]
const int NC = 128;
const int GHST = 1;
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
double U[NN][NN][2], UN[NN][NN][2], UF[NN][NN][2], UG[NG][NG][2];
double P[NN][NN], PN[NN][NN], PG[NG][NG];

// apply boundary conditions to velocity at cell center
void applyBCU(double u[NN][NN][2]) {
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        u[i][SWBOUND - 1][u_] = 2 * UBCs - u[i][SWBOUND    ][u_];
        u[i][NEBOUND    ][u_] = 2 * UBCn - u[i][NEBOUND - 1][u_];
        u[i][SWBOUND - 1][v_] = 2 * VBCs - u[i][SWBOUND    ][v_];
        u[i][NEBOUND    ][v_] = 2 * VBCn - u[i][NEBOUND - 1][v_];
    }
    for (int k = SWBOUND; k < NEBOUND; k ++) {
        u[SWBOUND - 1][k][u_] = 2 * UBCw - u[SWBOUND    ][k][u_];
        u[NEBOUND    ][k][u_] = 2 * UBCe - u[NEBOUND - 1][k][u_];
        u[SWBOUND - 1][k][v_] = 2 * VBCw - u[SWBOUND    ][k][v_];
        u[NEBOUND    ][k][v_] = 2 * VBCe - u[NEBOUND - 1][k][v_];
    }
}
// apply boundary conditions to pressure potential
void applyBCP(double p[NN][NN]) {
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        p[i][NEBOUND    ] = p[i][NEBOUND - 1] + DPBCn * dd;
        p[i][SWBOUND - 1] = p[i][SWBOUND    ] - DPBCs * dd;
    }
    for (int k = SWBOUND; k < NEBOUND; k ++) {
        p[NEBOUND    ][k] = p[NEBOUND - 1][k] + DPBCe * dd;
        p[SWBOUND - 1][k] = p[SWBOUND    ][k] - DPBCw * dd;
    }
}

// initialize velocity and pressure
void init(double u[NN][NN][2], double p[NN][NN]) {
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

    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            p[i][k] = 0;
        }
    }
    applyBCP(p);
}

// copy U to UN
void cpU(double u[NN][NN][2], double uN[NN][NN][2]) {
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            uN[i][k][u_] = u[i][k][u_];
            uN[i][k][v_] = u[i][k][v_];
        }
    }
}
// copy P to PN
void cpP(double p[NN][NN], double pN[NN][NN]) {
    for (int i = 0; i < NN; i ++) {
        for (int k = 0; k < NN; k ++) {
            pN[i][k] = p[i][k];
        }
    }
}

// 2nd-order central flux
double flux(double L, double R, double uF) {
    return uF * (L + R) * 0.5;
}

void fs1(double u[NN][NN][2], double uN[NN][NN][2]) {
    double fue, fuw, fun, fus, fve, fvw, fvn, fvs; // flux variables
    double ue, uw, uc, un, us, ve, vw, vc, vn, vs; // velocity at each direction
    double ufe, ufw, vfn, vfs;                     // velocity at cell faces
    double duudx, dvudy, duvdx, dvvdy;             // 1st derivatives
    double ddudxx, ddudyy, ddvdxx, ddvdyy;         // 2nd derivatives
    cpU(u, uN);
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            ue     = uN[i + 1][k    ][u_];
            uw     = uN[i - 1][k    ][u_];
            uc     = uN[i    ][k    ][u_];
            un     = uN[i    ][k + 1][u_];
            us     = uN[i    ][k - 1][u_];
            ve     = uN[i + 1][k    ][v_];
            vw     = uN[i - 1][k    ][v_];
            vc     = uN[i    ][k    ][v_];
            vn     = uN[i    ][k + 1][v_];
            vs     = uN[i    ][k - 1][v_];

            ufe    = 0.5 * (uc + ue);
            ufw    = 0.5 * (uw + uc);
            vfn    = 0.5 * (vc + vn);
            vfs    = 0.5 * (vs + vc);

            fue    = flux(uc, ue, ufe);
            fuw    = flux(uw, uc, ufw);
            fun    = flux(uc, un, vfn);
            fus    = flux(us, uc, vfs);
            fve    = flux(vc, ve, ufe);
            fvw    = flux(vw, vc, ufw);
            fvn    = flux(vc, vn, vfn);
            fvs    = flux(vs, vc, vfs);

            duudx  = (fue - fuw) / dd;
            dvudy  = (fun - fus) / dd;
            duvdx  = (fve - fvw) / dd;
            dvvdy  = (fvn - fvs) / dd;
            ddudxx = (uw - 2 * uc + ue) / (dd * dd);
            ddudyy = (us - 2 * uc + un) / (dd * dd);
            ddvdxx = (vw - 2 * vc + ve) / (dd * dd);
            ddvdyy = (vs - 2 * vc + vn) / (dd * dd);

            uc     = uc + dt * (- duudx - dvudy + (ddudxx + ddudyy) / Re);
            vc     = vc + dt * (- duvdx - dvvdy + (ddvdxx + ddvdyy) / Re);

            u[i][k][u_] = uc;
            u[i][k][v_] = vc;
        }
    }
}

int poisson(double p[NN][NN], double pN[NN][NN], double u[NN][NN][2]) {
    double E  = 0.02;
    int MAXIT = 10000;
    bool con  = false;
    int it    = 0;
    double ue, uw, vn, vs;     // velocity at each direction
    double pe, pw, pc, pn, ps; // pressure potential at each direction
    double dudx, dvdy;         // velocity gradient at cell center
    double ddpdxx, ddpdyy;     // 2nd derivatives at cell center
    double psi, res;
    while (con == false) {
        con = true;
        cpP(p, pN);
        for (int i = SWBOUND; i < NEBOUND; i ++) {
            for (int k = SWBOUND; k < NEBOUND; k ++) {
                ue      =  u[i + 1][k    ][u_];
                uw      =  u[i - 1][k    ][u_];
                vn      =  u[i    ][k + 1][v_];
                vs      =  u[i    ][k - 1][v_];
                pe      = pN[i + 1][k    ];
                pw      = pN[i - 1][k    ];
                pc      = pN[i    ][k    ];
                pn      = pN[i    ][k + 1];
                ps      = pN[i    ][k - 1];

                dudx    = (ue - uw) / (2 * dd);
                dvdy    = (vn - vs) / (2 * dd);
                ddpdxx  = (pw - 2 * pc + pe) / (dd * dd);
                ddpdyy  = (ps - 2 * pc + pn) / (dd * dd);
                psi     = (dudx + dvdy) / dt;
                res     = dtau * (ddpdxx + ddpdyy - psi);

                pc      = pc + res;
                if (abs(res) > abs(pc) * E) {
                    con = false;
                }

                p[i][k] = pc;
            }
        }
        applyBCP(p);
        it += 1;
        if (it >= MAXIT) {
            break;
        }
    }
    return it;
}

void fs2(double u[NN][NN][2], double p[NN][NN]) {
    double pe, pw, pn, ps; // pressure potential at each direction
    double dpdx, dpdy;     // pressure gradient at cell center
    double uc, vc;         // velocity at cell center
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            uc   = u[i    ][k    ][u_];
            vc   = u[i    ][k    ][v_];
            pe   = p[i + 1][k    ];
            pw   = p[i - 1][k    ];
            pn   = p[i    ][k + 1];
            ps   = p[i    ][k - 1];

            dpdx = (pe - pw) / (2 * dd);
            dpdy = (pn - ps) / (2 * dd);

            uc   = uc - dt * dpdx;
            vc   = vc - dt * dpdy;

            u[i][k][u_] = uc;
            u[i][k][v_] = vc;
        }
    }
    applyBCU(u);
}

double calcdiv(double u[NN][NN][2]) {
    double ue, uw, vn, vs;     // velocity at each direction
    double dudx, dvdy;         // velocity gradient at cell center
    double tdiv, div;
    tdiv = 0;
    for (int i = SWBOUND; i < NEBOUND; i ++) {
        for (int k = SWBOUND; k < NEBOUND; k ++) {
            ue   = u[i + 1][k    ][u_];
            uw   = u[i - 1][k    ][u_];
            vn   = u[i    ][k + 1][v_];
            vs   = u[i    ][k - 1][v_];
            dudx = (ue - uw) / (2 * dd);
            dvdy = (vn - vs) / (2 * dd);
            div  = abs(dudx + dvdy);
            if (div > tdiv) {
                tdiv = div;
            }
        }
    }
    return tdiv;
}

void ns2d(double u[NN][NN][2], double uN[NN][NN][2], double p[NN][NN], double pN[NN][NN]) {
    double t   = 0;
    double div = calcdiv(u);
    int it;
    printf("\rt = %5.8lf, poi = %5d, max div = %5.8lf", t, 0, div);
    fflush(stdout);
    while (t < ter) {
        fs1(u, uN);
        it = poisson(p, pN, u);
        fs2(u, p);
        div = calcdiv(u);
        t += dt;
        printf("\rt = %5.8lf, poi = %5d, max div = %5.8lf", t, it, div);
        fflush(stdout);
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
    uG[0     ][NG - 1][u_] = 0.5 * (UBCe + UBCs);
    uG[NG - 1][0     ][u_] = 0.5 * (UBCw + UBCn);
    uG[NG - 1][NG - 1][u_] = 0.5 * (UBCe + UBCn);
    uG[0     ][0     ][v_] = 0.5 * (VBCw + VBCs);
    uG[0     ][NG - 1][v_] = 0.5 * (VBCe + VBCs);
    uG[NG - 1][0     ][v_] = 0.5 * (VBCw + VBCn);
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
}

void o2f(double uG[NG][NG][2], double pG[NG][NG]) {
    FILE *fo;
    char fname[128];
    double xpos, ypos;
    sprintf(fname, "UVP_Re%d_t%d.ns2dcc8.csv", int(Re), int(ter));
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
    dt = min(CFL * dd * 0.5, Re * dd * dd * dd * dd / (2 * (dd * dd + dd * dd)));
    printf("Re = %lf, T = %lf, dt = %lf, dtau = %lf\n", Re, ter, dt, dtau);

    init(U, P);
    ns2d(U, UN, P, PN);
    calcgrid(U, UG, P, PG);
    o2f(UG, PG);

    return 0;
}
