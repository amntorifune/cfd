#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

using namespace std;

const int NC = 128;
const int NN = NC + 2;
const int NG = NC + 1;
const double cfl = 1.0;
const double dd = 1.0 / NC;
const double B0 = 4.0 / (dd * dd);
const double dtau = 1.0 / B0;
const double eps = 1.0;
const double ka = 1.0 / 3.0;
const double beta = (3 - ka) / (1 - ka);
double dt;
double Re;
double ter;

double U1[NN][NN], U2[NN][NN], UG[NG][NG];
double V1[NN][NN], V2[NN][NN], VG[NG][NG];
double P1[NN][NN], P2[NN][NN], PG[NG][NG];
double (*u)[NN] = U1;
double (*v)[NN] = V1;
double (*p)[NN] = P1;
double (*ut)[NN] = U2;
double (*vt)[NN] = V2;
double (*pt)[NN] = P2;
double (*ug)[NG] = UG;
double (*vg)[NG] = VG;
double (*pg)[NG] = PG;

int maxdiv;
int iter;

void bfu(void) {
    for (int i = 0; i < NN; i ++) {
        ut[i][0] = u[i][0];
        ut[i][NN - 1] = u[i][NN - 1];
    }
    for (int j = 0; j < NN; j ++) {
        ut[0][j] = u[0][j];
        ut[NN - 1][j] = u[NN - 1][j];
    }
}

void bfv(void) {
    for (int i = 0; i < NN; i ++) {
        vt[i][0] = v[i][0];
        vt[i][NN - 1] = v[i][NN - 1];
    }
    for (int j = 0; j < NN; j ++) {
        vt[0][j] = v[0][j];
        vt[NN - 1][j] = v[NN - 1][j];
    }
}

void bfp(void) {
    for (int i = 0; i < NN; i ++) {
        pt[i][0] = p[i][0];
        pt[i][NN - 1] = p[i][NN - 1];
    }
    for (int j = 0; j < NN; j ++) {
        pt[0][j] = p[0][j];
        pt[NN - 1][j] = p[NN - 1][j];
    }
}

int sgn(double x) {
    if (x > 0) {
        return 1;
    }
    else if (x < 0) {
        return -1;
    }
    else {
        return 0;
    }
}

// minmod(x, y) = sgn(x) * max{0, min{|x|, sgn(x) * y}}
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
    double DP = minmod(dP, beta * dM);
    double DM = minmod(dM, beta * dP);
    return c - eps * ((1 - ka) * DP + (1 + ka) * DM) / 4.0;
}

// interpolated value φ- at cell face: φA φB | φC φD
// φ- = φB + ε/4 * ((1 - κ) * Δ- + (1 + κ) * Δ+)
// Δ+ = minmod(d+, β * d-), Δ- = minmod(d-, β * d+)
// d+ = φC - φB, d- = φB - φA
double interpolM(double a, double b, double c) {
    double dP = c - b;
    double dM = b - a;
    double DP = minmod(dP, beta * dM);
    double DM = minmod(dM, beta * dP);
    return b + eps * ((1 - ka) * DM + (1 + ka) * DP) / 4.0;
}

// numerical flux at cell face: φA φB | φC φD
// f~ = 1/2 * ((f+ + f-) - |u@face| * (φ+ - φ-))
// f+- = u@face * φ+-
double flux(double a, double b, double c, double d, double uf) {
    double phiP = interpolP(b, c, d);
    double phiM = interpolM(a, b, c);
    double flxP = uf * phiP;
    double flxM = uf * phiM;
    return (flxP + flxM - abs(uf) * (phiP - phiM)) / 2.0;
}

// du/dt + duu/dx + dvu/dy = 1/Re * (ddu/dxx + ddu/dyy)
// dv/dt + duv/dx + dvv/dy = 1/Re * (ddv/dxx + ddv/dyy)
// equations above can be written as: dφ/dt = - (dF/dx + dG/dy) + viscousity
// F = u * φ, G = v * φ
// dF/dx = (fe~ - fw~) / Δx
// dG/dy = (fn~ - fs~) / Δy
//     fn
// fw cell fe
//     fs
void fs1(void) {
    double a, b, c, d, Uf;
    for (int i = 1; i < NN - 1; i ++) {
        for (int j = 1; j < NN - 1; j ++) {
            // fw~ = (uu)w~, φ = u
            //    = 1/2 * ((uu)w+ + (uu)w- - |uw| * (uw+ - uw-))
            //    = f~(ui-2, ui-1, ui, ui+1, uw)
            // (uu)w+- = fw+- = uw * uw+-
            // uw = u@wface = (ui-1 + ui) / 2
            // uw+ = u+(ui-1, ui, ui+1), uw- = u-(ui-2, ui-1, ui)
            double fw = 0.0;
            if (i > 1) {
                a = u[i - 2][j];
                b = u[i - 1][j];
                c = u[i][j];
                d = u[i + 1][j];
                Uf = (u[i - 1][j] + u[i][j]) / 2.0;
                fw = flux(a, b, c, d, Uf);
            }
            // fe~ = (uu)e~, φ = u
            //    = 1/2 * ((uu)e+ + (uu)e- - |ue| * (ue+ - ue-))
            //    = f~(ui-1, ui, ui+1, ui+2, ue)
            // (uu)e+- = fe+- = ue * ue+-
            // ue = u@eface = (ui + ui+1) / 2
            // ue+ = u+(ui, ui+1, ui+2), ue- = u-(ui-1, ui, ui+1)
            double fe = 0.0;
            if (i < NN - 2) {
                a = u[i - 1][j];
                b = u[i][j];
                c = u[i + 1][j];
                d = u[i + 2][j];
                Uf = (u[i][j] + u[i + 1][j]) / 2.0;
                fe = flux(a, b, c, d, Uf);
            }
            // fs~ = (vu)s~, φ = u
            //    = 1/2 * ((vu)s+ + (vu)s- - |vs| * (us+ - us-))
            //    = f~(uj-2, uj-1, uj, uj+1, vs)
            // (vu)s+- = fs+- = vs * us+-
            // vs = v@sface = (vj-1 + vj) / 2
            // us+ = u+(uj-1, uj, uj+1), us- = u-(uj-2, uj-1, uj)
            double fs = 0.0;
            if (j > 1) {
                a = u[i][j - 2];
                b = u[i][j - 1];
                c = u[i][j];
                d = u[i][j + 1];
                Uf = (v[i][j - 1] + v[i][j]) / 2.0;
                fs = flux(a, b, c, d, Uf);
            }
            // fn~ = (vu)n~, φ = u
            //    = 1/2 * ((vu)n+ + (vu)n- - |vn| * (un+ - un-))
            //    = f~(uj-1, uj, uj+1, uj+2, vn)
            // (vu)n+- = fn+- = vn * un+-
            // vn = v@nface = (vj + vj+1) / 2
            // un+ = u+(uj, uj+1, uj+2), un- = u-(uj-1, uj, uj+1)
            double fn = 0.0;
            if (j < NN - 2) {
                a = u[i][j - 1];
                b = u[i][j];
                c = u[i][j + 1];
                d = u[i][j + 2];
                Uf = (v[i][j] + v[i][j + 1]) / 2.0;
                fn = flux(a, b, c, d, Uf);
            }
            double duudx = (fe - fw) / dd;
            double dvudy = (fn - fs) / dd;
            double ddudxx = (u[i - 1][j] - 2 * u[i][j] + u[i + 1][j]) / (dd * dd);
            double ddudyy = (u[i][j - 1] - 2 * u[i][j] + u[i][j + 1]) / (dd * dd);
            ut[i][j] = u[i][j] + dt * (- duudx - dvudy + (ddudxx + ddudyy) / Re);
        }
    }
    for (int i = 1; i < NN - 1; i ++) {
        for (int j = 1; j < NN - 1; j ++) {
            // fw~ = (uv)w~, φ = v
            //    = 1/2 * ((uv)w+ + (uv)w- - |uw| * (vw+ - vw-))
            //    = f~(vi-2, vi-2, vi, vi+1, uw)
            // (uv)w+- = fw+- = uw * vw+-
            // uw = u@wface = (ui-1 + ui) / 2
            // vw+ = v+(vi-1, vi, vi+1), vw- = v-(vi-2, vi-1, vi)
            double fw = 0.0;
            if (i > 1) {
                a = v[i - 2][j];
                b = v[i - 1][j];
                c = v[i][j];
                d = v[i + 1][j];
                Uf = (u[i - 1][j] + u[i][j]) / 2.0;
                fw = flux(a, b, c, d, Uf);
            }
            // fe~ = (uv)e~, φ = v
            //    = 1/2 * ((uv)e+ + (uv)e- - |ue| * (ve+ - ve-))
            //    = f~(vi-1, vi, vi+1, vi+2)
            // (uv)e+- = fe+- = ue * ve+-
            // ue = u@eface = (ui + ui+1) / 2
            // ve+ = v+(vi, vi+1, vi+2), ve- = v-(vi-1, vi, vi+1)
            double fe = 0.0;
            if (i < NN - 2) {
                a = v[i - 1][j];
                b = v[i][j];
                c = v[i + 1][j];
                d = v[i + 2][j];
                Uf = (u[i][j] + u[i + 1][j]) / 2.0;
                fe = flux(a, b, c, d, Uf);
            }
            // fs~ = (vv)s~, φ = v
            //    = 1/2 * ((vv)s+ + (vv)s- - |vs| * (vs+ - vs-))
            //    = f~(vj-2, vj-1, vj, vj+1)
            // (vv)s+- = fs+- = vs * vs+-
            // vs = v@sface = (vj-1 + vj) / 2
            // vs+ = v+(vj-1, vj, vj+1), vs- = v-(vj-2, vj-1, vj)
            double fs = 0.0;
            if (j > 1) {
                a = v[i][j - 2];
                b = v[i][j - 1];
                c = v[i][j];
                d = v[i][j + 1];
                Uf = (v[i][j - 1] + v[i][j]) / 2.0;
                fs = flux(a, b, c, d, Uf);
            }
            // fn~ = (vv)n~, φ = v
            //    = 1/2 * ((vv)n+ + (vv)n- - |vn| * (vn+ - vn-))
            //    = f~(vj-1, vj, vj+1, vj+2)
            // (vv)n+- = fn+- = vn * vn+-
            // vn = v@nface = (vj + vj+1) / 2
            // vn+ = v+(vj, vj+1, vj+2), vn- = v-(vj-1, vj, vj+1)
            double fn = 0.0;
            if (j < NN - 2) {
                a = v[i][j - 1];
                b = v[i][j];
                c = v[i][j + 1];
                d = v[i][j + 2];
                Uf = (v[i][j] + v[i][j + 1]) / 2.0;
                fn = flux(a, b, c, d, Uf);
            }
            double duvdx = (fe - fw) / dd;
            double dvvdy = (fn - fs) / dd;
            double ddvdxx = (v[i - 1][j] - 2 * v[i][j] + v[i + 1][j]) / (dd * dd);
            double ddvdyy = (v[i][j - 1] - 2 * v[i][j] + v[i][j + 1]) / (dd * dd);
            vt[i][j] = u[i][j] + dt * (- duvdx - dvvdy + (ddvdxx + ddvdyy) / Re);
        }
    }
}

int main(int argc, char ** argv) {
    return 0;
}