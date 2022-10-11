#ifndef _POISSON_H
#define _POISSON_H 1

#include <stdio.h>
#include "matrix.h"
#include "mesh.h"
#include "param.h"

static void poisson_prepare_eq(Mat<real_t> &a, Mat<real_t> &rhs, Mesh &mesh, Dom &dom, stencil_t s) {
    if (s == stencil_t::d1s3) {
        #pragma acc kernels loop independent present(a, rhs, mesh, dom)
        for (int i = 0; i < dom.num; i ++) {
            if (i == 0) {
                int    i0 = (i - 1 >= 0      )? (i - 1) : i;
                int    i1 = i;
                int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
                real_t h1 = mesh.h.get(i1, 0);
                real_t x1 = mesh.x.get(i1, 0);
                real_t x2 = mesh.x.get(i2, 0);
                real_t c2 = kappa / (h1 * (x2 - x1));
                real_t c1 = - (c2 + 2 * kappa / (h1 * h1));
                real_t c0 = 0;
                real_t r1 = - 2 * kappa * TL / (h1 * h1);
                a.get(i1, 0) = c0;
                a.get(i1, 1) = c1;
                a.get(i1, 2) = c2;
                rhs.get(i1)  = r1;
            } else if (i == dom.num - 1) {
                int    i0 = (i - 1 >= 0      )? (i - 1) : i;
                int    i1 = i;
                int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
                real_t h1 = mesh.h.get(i1, 0);
                real_t x0 = mesh.x.get(i0, 0);
                real_t x1 = mesh.x.get(i1, 0);
                real_t c0 = kappa / (h1 * (x1 - x0));
                real_t c1 = - (c0 + 2 * kappa / (h1 * h1));
                real_t c2 = 0;
                real_t r1 = - 2 * kappa * TR / (h1 * h1);
                a.get(i1, 0) = c0;
                a.get(i1, 1) = c1;
                a.get(i1, 2) = c2;
                rhs.get(i1)  = r1;
            } else {
                int    i0 = (i - 1 >= 0      )? (i - 1) : i;
                int    i1 = i;
                int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
                real_t h1 = mesh.h.get(i1, 0);
                real_t x0 = mesh.x.get(i0, 0);
                real_t x1 = mesh.x.get(i1, 0);
                real_t x2 = mesh.x.get(i2, 0);
                real_t c2 = kappa / (h1 * (x2 - x1));
                real_t c0 = kappa / (h1 * (x1 - x0));
                real_t c1 = - (c2 + c0);
                real_t r1 = 0;
                a.get(i1, 0) = c0;
                a.get(i1, 1) = c1;
                a.get(i1, 2) = c2;
                rhs.get(i1)  = r1;
            }
        }
    }
}

static void poisson_sor(Mat<real_t> &a, Mat<real_t> &t, Mat<real_t> &rhs, Mat<real_t> &res, real_t &norm, Mesh &mesh, Dom &dom) {
    poisson_prepare_eq(a, rhs, mesh, dom, stencil_t::d1s3);
    int it = 0;
    do {
        #pragma acc kernels loop independent present(a, t, rhs, res, mesh, dom)
        for (int i = 0; i <  dom.num; i += 2) {
            int    i0 = (i - 1 >= 0      )? (i - 1) : i;
            int    i1 = i;
            int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
            real_t _0 = t.get(i0);
            real_t _1 = t.get(i1);
            real_t _2 = t.get(i2);
            real_t c0 = a.get(i1, 0);
            real_t c1 = a.get(i1, 1);
            real_t c2 = a.get(i1, 2);
            real_t _d = (rhs.get(i1) - c0 * _0 - c1 * _1 - c2 * _2) / c1;
            t.get(i1) = _1 + omega * _d;
        }
        #pragma acc kernels loop independent present(a, t, rhs, res, mesh, dom)
        for (int i = 1; i <  dom.num; i += 2) {
            int    i0 = (i - 1 >= 0      )? (i - 1) : i;
            int    i1 = i;
            int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
            real_t _0 = t.get(i0);
            real_t _1 = t.get(i1);
            real_t _2 = t.get(i2);
            real_t c0 = a.get(i1, 0);
            real_t c1 = a.get(i1, 1);
            real_t c2 = a.get(i1, 2);
            real_t _d = (rhs.get(i1) - c0 * _0 - c1 * _1 - c2 * _2) / c1;
            t.get(i1) = _1 + omega * _d;
        }
        calc_res(a, t, rhs, res, norm, dom, stencil_t::d1s3);
        printf("\r%6d %10.3e", ++it, norm);
    } while (norm > 1e-6);
    printf("\n");
}

#endif