#include "topo.h"
#include "flag.h"
#include "var.h"
#include "func.h"
#include <stdio.h>
#include <stdlib.h>

void setupf(void) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            F[i][j] = NONFLUID;
        }
    }

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            F[i][j] = FLUID;
        }
    }

    for (int i = I3; i <= I4; i ++) {
        for (int j = J3; j <= J4; j ++) {
            F[i][j] = NONFLUID;
        }
    }
}

void setupbb(void) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            FFU[i][j][0]    = NOTHING;
            FFU[i][j][1]    = NOTHING;
            BBU[i][j][0][0] = 0.0;
            BBU[i][j][0][1] = 0.0;
            BBU[i][j][1][0] = 0.0;
            BBU[i][j][1][1] = 0.0;

            FFP[i][j][0]    = NOTHING;
            FFP[i][j][1]    = NOTHING;
            BBP[i][j][0]    = 0.0;
            BBP[i][j][1]    = 0.0;
            BPi[i][j][0]    = 0.0;
            BPi[i][j][1]    = 0.0;
        }
    }

//  inflow boundary : (u,v,w) = (1,0,0), dp/dx = 0

    for (int j = J1; j <= J2; j ++) {
        FFU[I1 - 1][j][0]    = BOUNDARY;
        BBU[I1 - 1][j][0][0] = 1;
        BBU[I1 - 1][j][0][1] = 0;

        FFP[I1 - 1][j][0]    = BOUNDARY;

    //  apply zero gradient condition to pressure

    }

//  outflow boundary : (u,v,w) = (outflow), dp/dx = 0

    for (int j = J1; j <= J2; j ++) {
        FFU[I2][j][0] = BOUNDARY;
        BBU[I2][j][0][0] = 1;
        BBU[I2][j][0][1] = 0;

    //  apply zero gradient to velocity

        FFP[I2][j][0] = BOUNDARY;
        BBP[I2][j][0] = 0;
        BPi[I2][j][0] = 0;
    }

//  north boundary : (du/dy,v,dw/dy) = (0,0,0), dp/dy = 0

    for (int i = I1; i <= I2; i ++) {
        FFU[i][J2][1]    = BOUNDARY;

    //  apply slip wall condition to velocity

        FFP[i][J2][1]    = BOUNDARY;

    //  apply zero gradient condition to pressure

    }

//  south boundary : (du/dy,v,dw/dy) = (0,0,0), dp/dy = 0

    for (int i = I1; i <= I2; i ++) {
        FFU[i][J1 - 1][1]    = BOUNDARY;

    //  apply slip wall condition to velocity

        FFP[i][J1 - 1][1]    = BOUNDARY;

    //  apply zero gradient condition to pressure

    }

//  square east face : (u,v,w) = (0,0,0), dp/dx = [dp/dx]wall

    for (int j = J3; j <= J4; j ++) {
            FFU[I4][j][0]    = BOUNDARY;
            BBU[I4][j][0][0] = 0.0;
            BBU[I4][j][0][1] = 0.0;

            FFP[I4][j][0]    = BOUNDARY;

    //  apply wall zero gradient to pressure

    }

//  square west face : (u,v,w) = (0,0,0), dp/dx = [dp/dx]wall

    for (int j = J3; j <= J4; j ++) {
            FFU[I3 - 1][j][0]    = BOUNDARY;
            BBU[I3 - 1][j][0][0] = 0.0;
            BBU[I3 - 1][j][0][1] = 0.0;

            FFP[I3 - 1][j][0]    = BOUNDARY;

    //  apply wall zero gradientto pressure

    }

//  square north face : (u,v,w) = (0,0,0), dp/dy = [dp/dy]wall

    for (int i = I3; i <= I4; i ++) {
            FFU[i][J4][1]    = BOUNDARY;
            BBU[i][J4][1][0] = 0.0;
            BBU[i][J4][1][1] = 0.0;

            FFP[i][J4][1]    = BOUNDARY;

    //  apply wall zero gradient to pressure

    }

//  square south face : (u,v,w) = (0,0,0), dp/dy = [dp/dy]wall

    for (int i = I3; i <= I4; i ++) {
            FFU[i][J3 - 1][1]    = BOUNDARY;
            BBU[i][J3 - 1][1][0] = 0.0;
            BBU[i][J3 - 1][1][1] = 0.0;

            FFP[i][J3 - 1][1]    = BOUNDARY;

    //  apply wall zero gradient to pressure

    }
}

void setupx(void) {
    double DX, DY;

    DX = (X2 - X1) / NX;
    DY = (Y2 - Y1) / NY;

    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            double XX, YY;

            XX = (i - I1) * DX + 0.5 * DX;
            YY = (j - J1) * DY + 0.5 * DY;

            X[i][j][0] = XX;
            X[i][j][1] = YY;
        }
    }

    for (int i = I1 - 1; i <= I2 + 1; i ++) {
        for (int j = J1 - 1; j <= J2 + 1; j ++) {
            double XC0, XE1, XW1;
            double YC0, YN1, YS1;
            double X1K1, X2K2;
            double X1K11, X2K22;
            double K1X1, K2X2;
            double C1, C2, C7, C8;
            double G11, G22;
            double DET;

            XC0   = X[i    ][j    ][0];
            YC0   = X[i    ][j    ][1];
            XE1   = X[i + 1][j    ][0];
            XW1   = X[i - 1][j    ][0];
            YN1   = X[i    ][j + 1][1];
            YS1   = X[i    ][j - 1][1];

            X1K1  = 0.5 * (XE1 - XW1);
            X2K2  = 0.5 * (YN1 - YS1);
            X1K11 = XE1 - 2 * XC0 + XW1;
            X2K22 = YN1 - 2 * YC0 + YS1;
            K1X1  = 1 / X1K1;
            K2X2  = 1 / X2K2;
            DET   = X1K1 * X2K2;

            C1    = 1 / (X1K1 * X1K1);
            C2    = 1 / (X2K2 * X2K2);
            C7    = X1K11 / (X1K1 * X1K1 * X1K1);
            C8    = X2K22 / (X2K2 * X2K2 * X2K2);
                
            G11   = DET * K1X1 * K1X1;
            G22   = DET * K2X2 * K2X2;

            KX[i][j][0] = K1X1;
            KX[i][j][1] = K2X2;
            C[i][j][0]  = C1;
            C[i][j][1]  = C2;
            C[i][j][2]  = C7;
            C[i][j][3]  = C8;
            G[i][j][0]  = G11;
            G[i][j][1]  = G22;
            J[i][j]     = DET;
        }
    }
}

void setupu(void) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            if (F[i][j] == FLUID) {
                U[i][j][0] = 1;
                U[i][j][1] = 0;
            }
        }
    }
}

void init(void) {
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
                U[i][j][0]      = 0;
                U[i][j][1]      = 0;
                UA[i][j][0]     = 0;
                UA[i][j][1]     = 0;
                UD[i][j][0]     = 0;
                UD[i][j][1]     = 0;
                UP[i][j][0]     = 0;
                UP[i][j][1]     = 0;
                UC[i][j][0]     = 0;
                UC[i][j][1]     = 0;
                UU[i][j][0]     = 0;
                UU[i][j][1]     = 0;
                UUA[i][j][0]    = 0;
                UUA[i][j][1]    = 0;
                UUD[i][j][0]    = 0;
                UUD[i][j][1]    = 0;
                UUP[i][j][0]    = 0;
                UUP[i][j][1]    = 0;
                FFU[i][j][0]    = 0;
                FFU[i][j][1]    = 0;
                BBU[i][j][0][0] = 0;
                BBU[i][j][0][1] = 0;
                BBU[i][j][1][0] = 0;
                BBU[i][j][1][1] = 0;

                P[i][j]         = 0;
                PD[i][j]        = 0;
                FFP[i][j][0]    = 0;
                FFP[i][j][1]    = 0;
                BBP[i][j][0]    = 0;
                BBP[i][j][1]    = 0;
                Pi[i][j]        = 0;
                PiD[i][j]       = 0;
                BPi[i][j][0]    = 0;
                BPi[i][j][1]    = 0;

                DVR[i][j]       = 0;
                DVA[i][j]       = 0;
                DVP[i][j]       = 0;

                X[i][j][0]      = 0;
                X[i][j][1]      = 0;
                KX[i][j][0]     = 0;
                KX[i][j][1]     = 0;
                J[i][j]         = 0;
                G[i][j][0]      = 0;
                G[i][j][1]      = 0;
                C[i][j][0]      = 0;
                C[i][j][1]      = 0;
                C[i][j][2]      = 0;
                C[i][j][3]      = 0;
                SGS[i][j]       = 0;
                F[i][j]         = 0;
        }
    }

    setupx();
    setupf();
    setupbb();
    setupu();
}

void initpi(void) {
    #pragma acc kernels loop present(Pi)
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
                Pi[i][j] = 0;
        }
    }
}

void initp(void) {
    #pragma acc kernels loop present(P)
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
                P[i][j] = 0;
        }
    }
}

void paraout(void) {
    FILE *fo;
    fo = fopen("para.csv", "w+t");
    if (fo == NULL) {
        return;
    }
    else {
        fprintf(fo, "x,y,z,f,fxu,fyu,fxp,fyp,bxu,bxv,byu,byv,bxp,byp,k1x1,k2x2,j,c1,c2,c7,c8,g11,g22\n");
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
            for (int i = I1 - 2; i <= I2 + 2; i ++) {
                double XX, YY;
                int    FF;
                int    FXU, FYU;
                int    FXP, FYP;
                double BXU, BXV;
                double BYU, BYV;
                double BXP, BYP;
                double KX1, KX2;
                double DET;
                double C1, C2, C7, C8;
                double G11, G22;

                XX  =   X[i][j][0];
                YY  =   X[i][j][1];
                FF  =   F[i][j];
                FXU = FFU[i][j][0];
                FYU = FFU[i][j][1];
                FXP = FFP[i][j][0];
                FYP = FFP[i][j][1];
                BXU = BBU[i][j][0][0];
                BXV = BBU[i][j][0][1];
                BYU = BBU[i][j][1][0];
                BYV = BBU[i][j][1][1];
                BXP = BBP[i][j][0];
                BYP = BBP[i][j][1];
                KX1 =  KX[i][j][0];
                KX2 =  KX[i][j][1];
                DET =   J[i][j];
                C1  =   C[i][j][0];
                C2  =   C[i][j][1];
                C7  =   C[i][j][2];
                C8  =   C[i][j][3];
                G11 =   G[i][j][0];
                G22 =   G[i][j][1];
                fprintf(fo, "%lf,%lf,%lf,%d,%d,%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", XX, YY, 0.0, FF, FXU, FYU, FXP, FYP, BXU, BXV, BYU, BYV, BXP, BYP, KX1, KX2, DET, C1, C2, C7, C8, G11, G22);
            }
        }
    }
}

void fio3(char* fname) {
    #pragma acc data present(X, U, UA, UU, UUA, P, Pi, DVR, DVA, DVP, J, C)
    // acc data starts
    {
    #pragma acc update host(X, U, UA, UU, UUA, P, Pi, DVR, DVA, DVP, J, C)

    double XX, YY, PP, PPD, DD, DA, DP, PPi;
    double U1 , U2 , UU1 , UU2 ;
    double UA1, UA2, UUA1, UUA2;
    double UP1, UP2, UUP1, UUP2;
    double UD1, UD2, UUD1, UUD2;
    double J0J, C1, C2, C7, C8;
    double K1X1, K2X2;
    FILE  *fo;

    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,ua,va,ju,jv,jua,jva,p,pi,div,dva,dvp,j,c1,c2,c7,c8,k1,k2\n");
        for (int j = J1; j <= J2; j ++) {
            for (int i = I1; i <= I2; i ++) {
                XX   =   X[i][j][0];
                YY   =   X[i][j][1];
                U1   =   U[i][j][0];
                U2   =   U[i][j][1];
                UU1  =  UU[i][j][0];
                UU2  =  UU[i][j][1];
                UA1  =  UA[i][j][0];
                UA2  =  UA[i][j][1];
                UUA1 = UUA[i][j][0];
                UUA2 = UUA[i][j][1];
                PP   =   P[i][j];
                PPi  =  Pi[i][j];
                DD   = DVR[i][j];
                DA   = DVA[i][j];
                DP   = DVP[i][j];
                J0J  =   J[i][j];
                C1   =   C[i][j][0];
                C2   =   C[i][j][1];
                C7   =   C[i][j][2];
                C8   =   C[i][j][3];
                K1X1 =  KX[i][j][0];
                K2X2 =  KX[i][j][1];
                fprintf(fo, "%.2lf,%.2lf,%.2lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf,%.18lf\n", XX, YY, 0.0, U1, U2, UA1, UA2, UU1, UU2, UUA1, UUA2, PP, PPi, DD, DA, DP, J0J, C1, C2, C7, C8, K1X1, K2X2);
            }
        }

        fclose(fo);
    }

    // acc data ends
    }
}

void fio2(char* fname) {
    #pragma acc data present(X, U, P, Pi, DVR, DVA, DVP)
    // acc data starts
    {
    #pragma acc update host(X, U, P, Pi, DVR, DVA, DVP)

    double XX, YY, PP, PPD, DD, DA, DP, PPi;
    double U1 , U2 , UU1 , UU2 ;
    double UA1, UA2, UUA1, UUA2;
    double UP1, UP2, UUP1, UUP2;
    double UD1, UD2, UUD1, UUD2;
    FILE  *fo;

    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,p,pi,div,dva,dvp\n");
        for (int j = J1; j <= J2; j ++) {
            for (int i = I1; i <= I2; i ++) {
                XX   =   X[i][j][0];
                YY   =   X[i][j][1];
                U1   =   U[i][j][0];
                U2   =   U[i][j][1];
                PP   =   P[i][j];
                PPi  =  Pi[i][j];
                DD   = DVR[i][j];
                DA   = DVA[i][j];
                DP   = DVP[i][j];
                fprintf(fo, "%.2lf,%.2lf,%.2lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf\n", XX, YY, 0.0, U1, U2, PP, PPi, DD, DA, DP);
            }
        }

        fclose(fo);
    }

    // acc data ends
    }
}

void fio(char* fname) {
    #pragma acc data present(X, U, UA, UP, UD, UU, UUA, UUP, UUD, P, PD, Pi, DVR, DVA, DVP)
    // acc data starts
    {
    #pragma acc update host(X, U, UA, UP, UD, UU, UUA, UUP, UUD, P, PD, Pi, DVR, DVA, DVP)

    double XX, YY, PP, PPD, DD, DA, DP, PPi;
    double U1 , U2 , UU1 , UU2 ;
    double UA1, UA2, UUA1, UUA2;
    double UP1, UP2, UUP1, UUP2;
    double UD1, UD2, UUD1, UUD2;
    FILE  *fo;

    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,ua,va,up,vp,ud,vd,ju,jv,jua,jva,jup,jvp,jud,jvd,p,pd,pi,div,dva,dvp\n");
        for (int j = J1; j <= J2; j ++) {
            for (int i = I1; i <= I2; i ++) {
                XX   =   X[i][j][0];
                YY   =   X[i][j][1];
                U1   =   U[i][j][0];
                U2   =   U[i][j][1];
                UU1  =  UU[i][j][0];
                UU2  =  UU[i][j][1];
                UA1  =  UA[i][j][0];
                UA2  =  UA[i][j][1];
                UP1  =  UP[i][j][0];
                UP2  =  UP[i][j][1];
                UUA1 = UUA[i][j][0];
                UUA2 = UUA[i][j][1];
                UUP1 = UUP[i][j][0];
                UUP2 = UUP[i][j][1];
                UD1  =  UD[i][j][0];
                UD2  =  UD[i][j][1];
                UUD1 = UUD[i][j][0];
                UUD2 = UUD[i][j][1];
                PP   =   P[i][j];
                PPD  =  PD[i][j];
                PPi  =  Pi[i][j];
                DD   = DVR[i][j];
                DA   = DVA[i][j];
                DP   = DVP[i][j];
                fprintf(fo, "%.2lf,%.2lf,%.2lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf,%.12lf\n", XX, YY, 0.0, U1, U2, UA1, UA2, UP1, UP2, UD1, UD2, UU1, UU2, UUA1, UUA2, UUP1, UUP2, UUD1, UUD2, PP, PPD, PPi, DD, DA, DP);
            }
        }

        fclose(fo);
    }

    // acc data ends
    }
}

void backup(void) {
    #pragma acc kernels loop independent collapse(2) present(U, UU, UD, UUD, P, PD)
    for (int i = I1 - 2; i <= I2 + 2; i ++) {
        for (int j = J1 - 2; j <= J2 + 2; j ++) {
                UD[i][j][0]  =  U[i][j][0];
                UD[i][j][1]  =  U[i][j][1];
                UUD[i][j][0] = UU[i][j][0];
                UUD[i][j][1] = UU[i][j][1];
                // PD[i][j]     =  P[i][j];
        }
    }
}

void p0avg(void) {
    double SUM = 0;
    int    CNT = 0;
    #pragma acc kernels loop independent collapse(2) reduction(+:SUM, CNT) present(P, F) copy(SUM, CNT)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (F[i][j] == FLUID) {
                SUM = SUM + P[i][j];
                CNT = CNT + 1;
            }
        }
    }
    double AVG = SUM / CNT;
    #pragma acc kernels loop independent collapse(2) present(P, F) copyin(AVG)
    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            if (F[i][j] == FLUID) {
                P[i][j] = P[i][j] - AVG;
            }
        }
    }
}

int main(void) {
    char   fname[128];
    int    NFile = 0;
    int    CNT   = 0;
    double ADP;

    init();

    #pragma acc enter data copyin(U, UA, UP, UD, UC, UU, UUA, UUP, UUD, FFU, BBU, P, PD, FFP, BBP, Pi, PiD, BPi, DVR, DVA, DVP, X, KX, J, G, C, SGS, F)
    // outflow(U, UU, BBU, KX, J, DT);
    slip(U, BBU);
    contra(U, UC, UU, FFU, BBU, X, KX, J);

    sprintf(fname, "./data/varsimple.csv.%d", NFile);
    fio2(fname);
    NFile++;

    //  time advancing
    for (STEP = 1; STEP <= NSTEP; STEP ++) {
        // if (STEP <= 2000) {
        //     DT = 0.005;
        // }
        // else {
        //     DT = 0.001;
        // }
        // inflow(BBU, STEP * DT);
        // backup();
        pseudo(U, UA, UU, FFU, BBU, SGS, KX, J, C, F, ALPHA, DT, RI);
        outflow(U, UU, BBU, X, J, DT);
        contra(UA, UC, UUA, FFU, BBU, X, KX, J);
        // diver(UUA, J, DVA, F, AD);

        correction1(UP, UA, P, FFP, BBP, KX, F, DT);
        // outflow(U, UU, BBU, KX, J, DT);
        slip(UP, BBU);
        correction2(UUP, UUA, FFU, BBU, P, FFP, BBP, X, KX, G, DT);
        // diver(UUP, J, DVR, F, AD);
        diver(UUP, J, DVP, F, AD);
        ADP = AD;

        initpi();
        p0gradient(Pi, BPi);
        // loop for divergence
        CNT = 0;
        do {
            ITER = 0;
            // initpi();
            // p0gradient(Pi, BPi);
            // loop for poisson
            do {
                // jacobi(Pi, PiD, FFP, BPi, DVR, C, F, OMEGA, DT, R);
                sor(Pi, FFP, BPi, DVP, C, F, OMEGA, DT, R);
                // sor(P, FFP, BBP, DVA, C, F, OMEGA, DT, R);
                // p0avg(P, F);
                // wallns(U, P, BBP, X, RI);
                p0gradient(Pi, BPi);
                // p0gradient(P, BBP);
                ITER += 1;
            } while (R > EPOI && ITER < MAXIT);
            // hsmac(Pi, DVR, C, F, DT);
            // p0gradient(Pi, BPi);

            correction1(U, UP, Pi, FFP, BPi, KX, F, DT);
            // outflow(U, UU, BBU, KX, J, DT);
            slip(U, BBU);
            correction2(UU, UUP, FFU, BBU, Pi, FFP, BPi, X, KX, G, DT);
            diver(UU, J, DVR, F, AD);
            // correction3(P, Pi);
            // p0avg();
            // p0gradient(P, BBP);
            
            CNT += 1;
            printf("\rs(%6d,%6d):p(%4d,%13.10lf),d(%13.10lf,%13.10lf)", STEP, CNT, ITER, R, ADP, AD);
            fflush(stdout);
            // if (CNT >= 100) {
            //     printf("\n");
            //     break;
            // }
            
        } while (AD > EDIV);

        correction3(P, Pi);
        p0avg();
        p0gradient(P, BBP);

        if (STEP % 100 == 0 || STEP <= 50) {
            // sprintf(fname, "var.csv.%d", NFile);
            // fio(fname);

            sprintf(fname, "./data/varsimple.csv.%d", NFile);
            fio2(fname);

            NFile++;
        }
        // if (STEP  == 2001) {
        //     // sprintf(fname, "var.csv.%d", NFile);
        //     // fio(fname);

        //     sprintf(fname, "var2001.csv");
        //     fio2(fname);
        // }
    }
END:
    printf("\n");

    sprintf(fname, "./data/varsimple.csv.%d", NFile);
    fio2(fname);

    printf("NFile=%d\n", NFile);

    #pragma acc exit data copyout(U, UA, UP, UD, UC, UU, UUA, UUP, UUD, FFU, BBU, P, PD, FFP, BBP, Pi, PiD, BPi, DVR, DVA, DVP, X, KX, J, G, C, SGS, F)

    paraout();
    // fio((char*)"var.csv");

    return 0;
}
