include <stdio.h>
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

// U[X, Y, u or v]
double U[NN][NN][2], Un[NN][NN][2], Ug[NG][NG];
double P[NN][NN], Pn[NN][NN], Pg[NG][NG];

// boundary condition for U
void bcu(void) {

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