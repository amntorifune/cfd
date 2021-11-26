#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <algorithm>

using namespace std;

#define u_ 0
#define v_ 1
#define x_ 0
#define y_ 1

const int NCX      = 256;
const int NCY      = 128;
const int GHST     = 2;
const int NNX      = NCX + 2 * GHST;
const int NNY      = NCY + 2 * GHST;
const int WBOUND   = GHST;
const int SBOUND   = GHST;
const int EBOUND   = NCX + GHST;
const int NBOUND   = NCY + GHST;
const int NGX      = NCX + 1;
const int NGY      = NCY + 1;

const int FLUID        = 0;
const int NONFLUID     = 1;

const int NONE         = 0;
const int NSLIP        = 1;
const int IFLOW        = 2;
const int OFLOW        = 3;

int cylinder[2][2] = {
    {50, 50},
    {82, 82}
};

int    FF[NNX][NNY];
int    BF[NNX][NNY][2];        

void setupGeom(int ff[NNX][NNY], int bf[NNX][NNY][2]) {
    for (int i = 0; i < NNX; i ++) {
        for (int k = 0; k < NNY; k ++) {
            ff[i][k] = NONFLUID;
        }
    }
    // the tunnel
    for (int i = WBOUND; i < EBOUND; i ++) {
        for (int k = SBOUND; k < NBOUND; k ++) {
            ff[i][k] = FLUID;
        }
    }
    // the square cylinder
    for (int i = cylinder[0][x_]; i < cylinder[1][x_]; i ++) {
        for (int k = cylinder[0][y_]; k < cylinder[1][y_]; k ++) {
            ff[i][k] = NONFLUID;
        }
    }
    
    for (int i = cylinder[0][x_]; i < cylinder[1][x_]; i ++) {
        for (int k = cylinder[0][y_]; k < cylinder[1][y_]; k ++) {
            bf[i][k][x_] = NONE;
            bf[i][k][y_] = NONE;
        }
    }
    for (int k = SBOUND; k < NBOUND; k ++) {
        // the in-flow boundary
        bf[WBOUND - 1][k][x_] = IFLOW;
        // the out-flow boundary
        bf[EBOUND - 1][k][x_] = OFLOW;
    }
    for (int i = WBOUND; i < EBOUND; i ++) {
        // the wall boundaries
        bf[i][SBOUND - 1][y_] = NSLIP;
        bf[i][NBOUND - 1][y_] = NSLIP;
    }
    // boundaries of the square cylinder
    for (int i = cylinder[0][x_]; i < cylinder[1][x_]; i ++) {
        bf[i][cylinder[0][y_] - 1][y_] = NSLIP;
        bf[i][cylinder[1][y_] - 1][y_] = NSLIP;
    }
    for (int k = cylinder[0][y_]; k < cylinder[1][y_]; k ++) {
        bf[cylinder[0][x_] - 1][k][x_] = NSLIP;
        bf[cylinder[1][x_] - 1][k][x_] = NSLIP;
    }
}

int main(void) {
    setupGeom(FF, BF);

    FILE* fo;
    fo = fopen("test.txt", "w+t");
    if (fo == NULL) {
        return 0;
    }
    else {
        for (int k = NNY - 1; k >= 0; k --) {
            for (int i = 0; i < NNX; i ++) {
                fprintf(fo, "%d", FF[i][k]);
            }
            fprintf(fo, "\n");
        }
        fprintf(fo, "\n");
        for (int k = NNY - 1; k >= 0; k --) {
            for (int i = 0; i < NNX; i ++) {
                fprintf(fo, "%d", BF[i][k][x_]);
            }
            fprintf(fo, "\n");
        }
        fprintf(fo, "\n");
        for (int k = NNY - 1; k >= 0; k --) {
            for (int i = 0; i < NNX; i ++) {
                fprintf(fo, "%d", BF[i][k][y_]);
            }
            fprintf(fo, "\n");
        }
    }
}