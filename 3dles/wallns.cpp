#include "topo.h"

void wallns(
    double   U[N1 + 4][N2 + 4][N3 + 4][3],
    double  UD[N1 + 4][N2 + 4][N3 + 4][3],
    double   P[N1 + 4][N2 + 4][N3 + 4],
    double BBP[N1 + 4][N2 + 4][N3 + 4][3],
    double   X[N1 + 4][N2 + 4][N3 + 4][3],
    double   RI,
    bool     Ph // whether the velocity is physical
) {

//  bottom boundary

    for (int i = I1; i <= I2; i ++) {
        for (int j = J1; j <= J2; j ++) {
            double LAP, DK3, DK33, C3, C9;
            double X3K3, X3K33;
            double ZF0, ZT1, ZB1;
            double WF0, WT1, WB1;
            double DPK3;

            ZF0   = Z1;
            ZT1   = X[i][j][K1    ][z];
            ZB1   = X[i][j][K1 - 1][z];
            X3K3  = ZT1 - ZB1;
            X3K33 = 4 * (ZT1 - 2 * ZF0 + ZB1);
            C3    = 1 / (X3K3 * X3K3);
            C9    = X3K33 / (X3K3 * X3K3 * X3K3);

            WF0   = 0;
            WT1   = U[i][j][K1][w];
            WB1   = (Ph)? (- WT1) : (- UD[i][j][K1][w]);
            DK3   = WT1 - WB1;
            DK33  = 4 * (WT1 - 2 * WF0 + WB1);
            LAP   = C3 * DK33 - C9 * DK3;
            DPK3  = X3K3 * RI * LAP;

            BBP[i][j][K1 - 1][z] = P[i][j][K1] - 0.5 * DPK3;
        }
    }

//  cube top face

    for (int i = I3; i <= I4; i ++) {
        for (int j = J3; j <= J4; j ++) {
            double LAP, DK3, DK33, C3, C9;
            double X3K3, X3K33;
            double ZF0, ZT1, ZB1;
            double WF0, WT1, WB1;
            double DPK3;

            ZF0   = Z4;
            ZT1   = X[i][j][K4 + 1][z];
            ZB1   = X[i][j][K4    ][z];
            X3K3  = ZT1 - ZB1;
            X3K33 = 4 * (ZT1 - 2 * ZF0 + ZB1);
            C3    = 1 / (X3K3 * X3K3);
            C9    = X3K33 / (X3K3 * X3K3 * X3K3);

            WF0   = 0;
            WT1   = U[i][j][K4 + 1][w];
            WB1   = (Ph)? (- WT1) : (- UD[i][j][K4 + 1][w]);
            DK3   = WT1 - WB1;
            DK33  = 4 * (WT1 - 2 * WF0 + WB1);
            LAP   = C3 * DK33 - C9 * DK3;
            DPK3  = X3K3 * RI * LAP;

            BBP[i][j][K4][z] = P[i][j][K4 + 1] - 0.5 * DPK3;
        }
    }

//  cube east face

    for (int j = J3; j <= J4; j ++) {
        for (int k = K3; k <= K4; k ++) {
            double LAP, DK1, DK11, C1, C7;
            double X1K1, X1K11;
            double XF0, XE1, XW1;
            double UF0, UE1, UW1;
            double DPK1;

            XF0   = X4;
            XE1   = X[I4 + 1][j][k][x];
            XW1   = X[I4    ][j][k][x];
            X1K1  = XE1 - XW1;
            X1K11 = 4 * (XE1 - 2 * XF0 + XW1);
            C1    = 1 / (X1K1 * X1K1);
            C7    = X1K11 / (X1K1 * X1K1 * X1K1);

            UF0   = 0;
            UE1   = U[I4 + 1][j][k][u];
            UW1   = (Ph)? (- UE1) : (- UD[I4 + 1][j][k][u]);
            DK1   = UE1 - UW1;
            DK11  = 4 * (UE1 - 2 * UF0 + UW1);
            LAP   = C1 * DK11 - C7 * DK1;
            DPK1  = X1K1 * RI * LAP;

            BBP[I4][j][k][x] = P[I4 + 1][j][k] - 0.5 * DPK1;
        }
    }

//  cube west face

    for (int j = J3; j <= J4; j ++) {
        for (int k = K3; k <= K4; k ++) {
            double LAP, DK1, DK11, C1, C7;
            double X1K1, X1K11;
            double XF0, XE1, XW1;
            double UF0, UE1, UW1;
            double DPK1;

            XF0   = X3;
            XE1   = X[I3    ][j][k][x];
            XW1   = X[I3 - 1][j][k][x];
            X1K1  = XE1 - XW1;
            X1K11 = 4 * (XE1 - 2 * XF0 + XW1);
            C1    = 1 / (X1K1 * X1K1);
            C7    = X1K11 / (X1K1 * X1K1 * X1K1);

            UF0   = 0;
            UW1   = U[I3 - 1][j][k][u];
            UE1   = (Ph)? (- UW1) : (- UD[I3 - 1][j][k][u]);
            DK1   = UE1 - UW1;
            DK11  = 4 * (UE1 - 2 * UF0 + UW1);
            LAP   = C1 * DK11 - C7 * DK1;
            DPK1  = X1K1 * RI * LAP;

            BBP[I3 - 1][j][k][x] = P[I3 - 1][j][k] + 0.5 * DPK1;
        }
    }

//  cube north face

    for (int i = I3; i <= I4; i ++) {
        for (int k = K3; k <= K4; k ++) {
            double LAP, DK2, DK22, C2, C8;
            double X2K2, X2K22;
            double YF0, YN1, YS1;
            double VF0, VN1, VS1;
            double DPK2;

            YF0   = Y4;
            YN1   = X[i][J4 + 1][k][y];
            YS1   = X[i][J4    ][k][y];
            X2K2  = YN1 - YS1;
            X2K22 = 4 * (YN1 - 2 * YF0 + YS1);
            C2    = 1 / (X2K2 * X2K2);
            C8    = X2K22 / (X2K2 * X2K2 * X2K2);

            VF0   = 0;
            VN1   = U[i][J4 + 1][k][v];
            VS1   = (Ph)? (- VN1) : (- UD[i][J4 + 1][k][v]);
            DK2   = VN1 - VS1;
            DK22  = 4 * (VN1 - 2 * VF0 + VS1);
            LAP   = C2 * DK22 - C8 * DK2;
            DPK2  = X2K2 * RI * LAP;

            BBP[i][J4][k][y] = P[i][J4 + 1][k] - 0.5 * DPK2;
        }
    }

//  cube south face

    for (int i = I3; i <= I4; i ++) {
        for (int k = K3; k <= K4; k ++) {
            double LAP, DK2, DK22, C2, C8;
            double X2K2, X2K22;
            double YF0, YN1, YS1;
            double VF0, VN1, VS1;
            double DPK2;

            YF0   = Y3;
            YN1   = X[i][J3    ][k][y];
            YS1   = X[i][J3 - 1][k][y];
            X2K2  = YN1 - YS1;
            X2K22 = 4 * (YN1 - 2 * YF0 + YS1);
            C2    = 1 / (X2K2 * X2K2);
            C8    = X2K22 / (X2K2 * X2K2 * X2K2);

            VF0   = 0;
            VS1   = U[i][J3 - 1][k][v];
            VN1   = (Ph)? (- VS1) : (- UD[i][J3 - 1][k][v]);
            DK2   = VN1 - VS1;
            DK22  = 4 * (VN1 - 2 * VF0 + VS1);
            LAP   = C2 * DK22 - C8 * DK2;
            DPK2  = X2K2 * RI * LAP;

            BBP[i][J3 - 1][k][y] = P[i][J3 - 1][k] + 0.5 * DPK2;
        }
    }

}
