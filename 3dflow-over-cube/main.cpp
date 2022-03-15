#include "var.h"
#include "func.h"
#include <stdio.h>
#include <stdlib.h>

void paraout(void) {
    FILE *fo;
    fo = fopen("para.csv", "w+t");
    if (fo == NULL) {
        return;
    }
    else {
        fprintf(fo, "x,y,z,f,fxu,fyu,fzu,fxp,fyp,fzp,bxu,bxv,bxw,byu,byv,byw,bzu,bzv,bzw,bxp,byp,bzp,k1x1,k2x2,k3x3,j,c1,c2,c3,c7,c8,c9,g11,g22,g33\n");
        for (int k = K1 - 2; k <= K2 + 2; k ++) {
            for (int j = J1 - 2; j <= J2 + 2; j ++) {
                for (int i = I1 - 2; i <= I2 + 2; i ++) {
                    double XX, YY, ZZ;
                    int    FF;
                    int    FXU;
                    int    FYU;
                    int    FZU;
                    int    FXP, FYP, FZP;
                    double BXU, BXV, BXW;
                    double BYU, BYV, BYW;
                    double BZU, BZV, BZW;
                    double BXP, BYP, BZP;
                    double KX1, KX2, KX3;
                    double DET;
                    double C1, C2, C3, C7, C8, C9;
                    double G11, G22, G33;

                    XX  =   X[i][j][k][x];
                    YY  =   X[i][j][k][y];
                    ZZ  =   X[i][j][k][z];
                    FF  =   F[i][j][k];
                    FXU = FFU[i][j][k][x];
                    FYU = FFU[i][j][k][y];
                    FZU = FFU[i][j][k][z];
                    FXP = FFP[i][j][k][x];
                    FYP = FFP[i][j][k][y];
                    FZP = FFP[i][j][k][z];
                    BXU = BBU[i][j][k][x][u];
                    BXV = BBU[i][j][k][x][v];
                    BXW = BBU[i][j][k][x][w];
                    BYU = BBU[i][j][k][y][u];
                    BYV = BBU[i][j][k][y][v];
                    BYW = BBU[i][j][k][y][w];
                    BZU = BBU[i][j][k][z][u];
                    BZV = BBU[i][j][k][z][v];
                    BZW = BBU[i][j][k][z][w];
                    BXP = BBP[i][j][k][x];
                    BYP = BBP[i][j][k][y];
                    BZP = BBP[i][j][k][z];
                    KX1 =  KX[i][j][k][x];
                    KX2 =  KX[i][j][k][y];
                    KX3 =  KX[i][j][k][z];
                    DET =   J[i][j][k];
                    C1  =   C[i][j][k][0];
                    C2  =   C[i][j][k][1];
                    C3  =   C[i][j][k][2];
                    C7  =   C[i][j][k][3];
                    C8  =   C[i][j][k][4];
                    C9  =   C[i][j][k][5];
                    G11 =   G[i][j][k][x];
                    G22 =   G[i][j][k][y];
                    G33 =   G[i][j][k][z];
                    fprintf(fo, "%lf,%lf,%lf,%d,%d,%d,%d,%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", XX, YY, ZZ, FF, FXU, FYU, FZU, FXP, FYP, FZP, BXU, BXV, BXW, BYU, BYV, BYW, BZU, BZV, BZW, BXP, BYP, BZP, KX1, KX2, KX3, DET, C1, C2, C3, C7, C8, C9, G11, G22, G33);
                }
            }
        }
    }
}

int main(void) {
    init(U, UA, UC, UU, UUA, FFU, BBU, P, FFP, BBP, DIV, X, KX, J, G, C, SGS, F);
    paraout();
}
