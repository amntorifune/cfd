#ifndef _POISSON_H
#define _POISSON_H 1

#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "mesh.h"
#include "param.h"

static real_t dot2(Matrix<real_t> &a, Matrix<real_t> &b, Dom &dom) {
    real_t r = 0;
    #pragma acc kernels loop independent collapse(3) reduction(+:r) present(a, b, dom) copy(r)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                r += a.get(Util::id(i,j,k,dom.size)) * b.get(Util::id(i,j,k,dom.size));
            }
        }
    }
    return r;
}

static void poisson_prepare_eq(Matrix<real_t> &a, Matrix<real_t> &b, Mesh &mesh, Dom &dom, stencil_t s) {
    if (s == stencil_t::d1s3) {
        #pragma acc kernels loop independent collapse(3) present(a, b, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if (i == 0) {
                        int    cur = Util::id(i,j,k,dom.size);
                        real_t dxr = mesh.x.get(Util::id(i+1,j,k,dom.size), 0) - mesh.x.get(cur, 0);
                        real_t dx  = mesh.h.get(cur, 0);
                        real_t cr  = kappa / (dx * dxr);
                        real_t cl  = 0;
                        real_t cc  = - (cr + 2 * kappa / (dx * dx));
                        real_t rc  = - 2 * kappa * TL / (dx * dx);
                        a.get(cur, 0) = cl;
                        a.get(cur, 1) = cc;
                        a.get(cur, 2) = cr;
                        b.get(cur)    = rc;
                    } else if (i == dom.size[0] - 1) {
                        int    cur = Util::id(i,j,k,dom.size);
                        real_t dxl = mesh.x.get(cur, 0) - mesh.x.get(Util::id(i-1,j,k,dom.size), 0);
                        real_t dx  = mesh.h.get(cur, 0);
                        real_t cr  = 0;
                        real_t cl  = kappa / (dx * dxl);
                        real_t cc  = - (cl + 2 * kappa / (dx * dx));
                        real_t rc  = - 2 * kappa * TR / (dx * dx);
                        a.get(cur, 0) = cl;
                        a.get(cur, 1) = cc;
                        a.get(cur, 2) = cr;
                        b.get(cur)    = rc;
                    } else {
                        int    cur = Util::id(i,j,k,dom.size);
                        real_t dxr = mesh.x.get(Util::id(i+1,j,k,dom.size), 0) - mesh.x.get(cur, 0);
                        real_t dxl = mesh.x.get(cur, 0) - mesh.x.get(Util::id(i-1,j,k,dom.size), 0);
                        real_t dx  = mesh.h.get(cur, 0);
                        real_t cr  = kappa / (dx * dxr);
                        real_t cl  = kappa / (dx * dxl);
                        real_t cc  = - (cr + cl);
                        real_t rc  = 0;
                        a.get(cur, 0) = cl;
                        a.get(cur, 1) = cc;
                        a.get(cur, 2) = cr;
                        b.get(cur)    = rc;
                    }
                }
            }
        }
    }
}

static void scale_eq(Matrix<real_t> &a, Matrix<real_t> &b, Dom &dom) {
    real_t max_diag = 0;
    #pragma acc kernels loop independent collapse(3) reduction(max:max_diag) present(a, dom) copy(max_diag)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int    cur = Util::id(i,j,k,dom.size);
                real_t cc  = a.get(cur, 1);
                if (fabs(cc) > max_diag) {
                    max_diag = fabs(cc);
                }
            }
        }
    }
    #pragma acc kernels loop independent collapse(3) present(a, b, dom) copyin(max_diag)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int    cur = Util::id(i,j,k,dom.size);
                for (int m = 0; m < a.col; m ++) {
                    a.get(cur, m) = a.get(cur, m) / max_diag;
                }
                b.get(cur) = b.get(cur) / max_diag;
            }
        }
    }
}

static void calc_res(Matrix<real_t> &a, Matrix<real_t> &x, Matrix<real_t> &b, Matrix<real_t> &r, real_t &norm, Mesh &mesh, Dom &dom, stencil_t s) {
    real_t sum = 0;
    if (s == stencil_t::d1s3) {
        /* #pragma acc kernels loop independent collapse(2) reduction(+:sum) present(a, x, b, r, dom) copy(sum)
        for (int i = 0; i < dom.num; i ++) {
            for (int m = 0; m < x.col; m ++) {
                int    i0 = (i - 1 >= 0      )? (i - 1) : i;
                int    i1 = i;
                int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
                real_t _0 = x.get(i0, m);
                real_t _1 = x.get(i1, m);
                real_t _2 = x.get(i2, m);
                real_t c0 = a.get(i1, 0);
                real_t c1 = a.get(i1, 1);
                real_t c2 = a.get(i1, 2);
                real_t r1 = b.get(i1) - (c0 * _0 + c1 * _1 + c2 * _2);
                r.get(i1, m) = r1;
                sum += r1 * r1;
            }
        } */
        #pragma acc kernels loop independent collapse(3) reduction(+:sum) present(a, x, b, r, dom, mesh) copy(sum)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    int    cur = Util::id(i,j,k,dom.size);
                    real_t cr  = a.get(cur, 2);
                    real_t cc  = a.get(cur, 1);
                    real_t cl  = a.get(cur, 0);
                    int    ir  = mesh.map.get(cur, 3);
                    int    il  = mesh.map.get(cur, 0);
                    int    ic  = cur;
                    for (int m = 0; m < x.col; m ++) {
                        real_t _r = x.get(ir, m);
                        real_t _l = x.get(il, m);
                        real_t _c = x.get(ic, m);
                        real_t rc = b.get(ic, m) - (_r * cr + _c * cc + _l * cl);
                        r.get(ic, m) = rc;
                        sum += rc * rc;
                    }
                }
            }
        }
    }
    norm = sqrt(sum / x.num);
}

static void poisson_sor(Matrix<real_t> &a, Matrix<real_t> &t, Matrix<real_t> &rhs, Matrix<real_t> &res, real_t &norm, Mesh &mesh, Dom &dom) {
    int it = 0;
    do {
        /* #pragma acc kernels loop independent present(a, t, rhs, res, mesh, dom)
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
        } */
        #pragma acc kernels loop independent collapse(3) present(a, t, rhs, res, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if ((i + j + k) % 2 == 0) {
                        int    cur = Util::id(i,j,k,dom.size);
                        real_t cr  = a.get(cur, 2);
                        real_t cc  = a.get(cur, 1);
                        real_t cl  = a.get(cur, 0);
                        int    ir  = mesh.map.get(cur, 3);
                        int    il  = mesh.map.get(cur, 0);
                        int    ic  = cur;
                        real_t _r  = t.get(ir);
                        real_t _l  = t.get(il);
                        real_t _c  = t.get(ic);
                        real_t _d  = (rhs.get(ic) - (_r * cr + _c * cc + _l * cl)) / cc;
                        t.get(ic)  = _c + sor_omega * _d;
                    }

                }
            }
        }
        #pragma acc kernels loop independent collapse(3) present(a, t, rhs, res, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if ((i + j + k) % 2 == 1) {
                        int    cur = Util::id(i,j,k,dom.size);
                        real_t cr  = a.get(cur, 2);
                        real_t cc  = a.get(cur, 1);
                        real_t cl  = a.get(cur, 0);
                        int    ir  = mesh.map.get(cur, 3);
                        int    il  = mesh.map.get(cur, 0);
                        int    ic  = cur;
                        real_t _r  = t.get(ir);
                        real_t _l  = t.get(il);
                        real_t _c  = t.get(ic);
                        real_t _d  = (rhs.get(ic) - (_r * cr + _c * cc + _l * cl)) / cc;
                        t.get(ic)  = _c + sor_omega * _d;
                    }

                }
            }
        }
        calc_res(a, t, rhs, res, norm, mesh, dom, stencil_t::d1s3);
        
        if (it % 1000 == 0) {
            t.update_self();
            FILE *fo;
            char fname[128];
            sprintf(fname, "temperature.csv.%d", it / 1000);
            fo = fopen(fname, "w+t");
            if (fo) {
                fprintf(fo, "x,y,z,t\n");
                real_t x, y, z, tp;
                for (int k = 0; k < dom.size[2]; k ++) {
                    for (int j = 0; j < dom.size[1]; j ++) {
                        for (int i = 0; i < dom.size[0]; i ++) {
                            int cur = Util::id(i,j,k,dom.size);
                            x  = mesh.x.get(cur, 0);
                            y  = mesh.x.get(cur, 1);
                            z  = mesh.x.get(cur, 2);
                            tp = t.get(cur);
                            fprintf(fo, "%10.3e,%10.3e,%10.3e,%10.3e\n", x, y, z, tp);
                        }
                    }
                }
            }
        }

        printf("\r%6d %10.3e", ++it, norm);
    } while (norm > 1e-9);
    printf("\n");
}

static void calc_ax(Matrix<real_t> &a, Matrix<real_t> &x, Matrix<real_t> &ax, Mesh &mesh, Dom &dom, stencil_t s) {
    if (s == stencil_t::d1s3) {
        #pragma acc kernels loop independent collapse(3) present(a, x, ax, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    int    cur = Util::id(i,j,k,dom.size);
                    real_t cr  = a.get(cur, 2);
                    real_t cc  = a.get(cur, 1);
                    real_t cl  = a.get(cur, 0);
                    int    ir  = mesh.map.get(cur, 3);
                    int    il  = mesh.map.get(cur, 0);
                    int    ic  = cur;
                    for (int m = 0; m < x.col; m ++) {
                        real_t _r = x.get(ir, m);
                        real_t _l = x.get(il, m);
                        real_t _c = x.get(ic, m);
                        ax.get(ic, m) = cr * _r + cc * _c + cl * _l;
                    }
                }
            }
        }
    }
}

static void pbicgstab_1(Matrix<real_t> &p, Matrix<real_t> &q, Matrix<real_t> &r, real_t beta, real_t omega, Dom &dom) {
    #pragma acc kernels loop independent collapse(3) present(p, q, r, dom) copyin(beta, omega)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int cur = Util::id(i,j,k,dom.size);
                p.get(cur) = r.get(cur) + beta * (p.get(cur) - omega * q.get(cur));
            }
        }
    }
}

static void pbicgstab_2(Matrix<real_t> &s, Matrix<real_t> &q, Matrix<real_t> &r, real_t alpha, Dom &dom) {
    #pragma acc kernels loop independent collapse(3) present(s, q, r, dom) copyin(alpha)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int cur = Util::id(i,j,k,dom.size);
                s.get(cur) = r.get(cur) - alpha * q.get(cur);
            }
        }
    }
}

static void pbicgstab_3(Matrix<real_t> &x, Matrix<real_t> &_p, Matrix<real_t> &_s, real_t alpha, real_t omega, Dom &dom) {
    #pragma acc kernels loop independent collapse(3) present(x, _p, _s, dom) copyin(alpha, omega)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int cur = Util::id(i,j,k,dom.size);
                x.get(cur) = x.get(cur) + alpha * _p.get(cur) + omega * _s.get(cur);
            }
        }
    }
}

static void pbicgstab_4(Matrix<real_t> &r, Matrix<real_t> &s, Matrix<real_t> &t, real_t omega, Dom &dom) {
    #pragma acc kernels loop independent collapse(3) present(r, s, t, dom) copyin(omega)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int cur = Util::id(i,j,k,dom.size);
                r.get(cur) = s.get(cur) - omega * t.get(cur);
            }
        }
    }
}

static void precondition_sor(Matrix<real_t> &a, Matrix<real_t> &x, Matrix<real_t> &b, Mesh &mesh, Dom &dom, int maxit) {
    for (int it = 0; it < maxit; it ++) {
        #pragma acc kernels loop independent collapse(3) present(a, x, b, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if ((i + j + k) % 2 == 0) {
                        int    cur = Util::id(i,j,k,dom.size);
                        real_t cr  = a.get(cur, 2);
                        real_t cc  = a.get(cur, 1);
                        real_t cl  = a.get(cur, 0);
                        int    ir  = mesh.map.get(cur, 3);
                        int    il  = mesh.map.get(cur, 0);
                        int    ic  = cur;
                        real_t _r  = x.get(ir);
                        real_t _l  = x.get(il);
                        real_t _c  = x.get(ic);
                        real_t _d  = (b.get(ic) - (_r * cr + _c * cc + _l * cl)) / cc;
                        x.get(ic)  = _c + 0.9 * _d;
                    }
                }
            }
        }
        #pragma acc kernels loop independent collapse(3) present(a, x, b, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if ((i + j + k) % 2 == 1) {
                        int    cur = Util::id(i,j,k,dom.size);
                        real_t cr  = a.get(cur, 2);
                        real_t cc  = a.get(cur, 1);
                        real_t cl  = a.get(cur, 0);
                        int    ir  = mesh.map.get(cur, 3);
                        int    il  = mesh.map.get(cur, 0);
                        int    ic  = cur;
                        real_t _r  = x.get(ir);
                        real_t _l  = x.get(il);
                        real_t _c  = x.get(ic);
                        real_t _d  = (b.get(ic) - (_r * cr + _c * cc + _l * cl)) / cc;
                        x.get(ic)  = _c + 0.9 * _d;
                    }
                }
            }
        }
    }
}

static real_t calc_norm(Matrix<real_t> &x) {
    real_t sum = 0;
    #pragma acc kernels loop independent reduction(+:sum) present(x) copy(sum)
    for (int i = 0; i < x.num; i ++) {
        real_t value = x.get(i);
        sum += ((value * value) / x.num);
    }
    return sqrt(sum);
}

static void poisson_pbicgstab(Matrix<real_t> &a, Matrix<real_t> &x, Matrix<real_t> &b, Matrix<real_t> &r, real_t &norm, Mesh &mesh, Dom &dom, stencil_t stencil) {
    Matrix<real_t> _r(dom, 1);
    Matrix<real_t>  p(dom, 1);
    Matrix<real_t>  q(dom, 1);
    Matrix<real_t>  s(dom, 1);
    Matrix<real_t> _p(dom, 1);
    Matrix<real_t> _s(dom, 1);
    Matrix<real_t>  t(dom, 1);
    _r.to_device();
    p.to_device();
    q.to_device();
    s.to_device();
    _p.to_device();
    _s.to_device();
    t.to_device();
    real_t rho, _rho, alpha, beta, omega, dummy;
    int it = 0;
    calc_res(a, x, b, r, dummy, mesh, dom, stencil);

    // t.update_self();
    // for(int i = 0; i < dom.num; i ++) {
    //     printf("%10.3e ", x.get(i));
    // }
    // r.update_self();
    // printf("\n");
    // for(int i = 0; i < dom.num; i ++) {
    //     printf("%10.3e ", r.get(i));
    // }
    // printf("\n");

    _r = r;
    // _r.update_self();
    // for(int i = 0; i < dom.num; i ++) {
    //     printf("%10.3e ", _r.get(i));
    // }
    // printf("\n");

    q     = 0;
    _rho  = 1;
    alpha = 0;
    omega = 1;
    do {
        rho = dot2(_r, r, dom);
        if (fabs(rho) < __FLT_MIN__) {
            printf("\nsmall rho\n");
            // r.update_self();
            // for(int i = 0; i < dom.num; i ++) {
            //     printf("%10.3e ", r.get(i));
            // }
            // printf("\n");

            // _r.update_self();
            // for(int i = 0; i < dom.num; i ++) {
            //     printf("%10.3e ", _r.get(i));
            // }
            printf("\n");
            printf("%10.3e\n", rho);
            break;
        }

        if (it == 0) {
            p = r;
        } else {
            beta = (rho / _rho) * (alpha / omega);
            pbicgstab_1(p, q, r, beta, omega, dom);
        }

        _p = 0;
        precondition_sor(a, _p, p, mesh, dom, 5);
        calc_ax(a, _p, q, mesh, dom, stencil);
        alpha = rho / dot2(_r, q, dom);
        pbicgstab_2(s, q, r, alpha, dom);
        _s = 0;
        precondition_sor(a, _s, s, mesh, dom, 5);
        calc_ax(a, _s, t, mesh, dom, stencil);
        omega = dot2(t, s, dom) / dot2(t, t, dom);
        pbicgstab_3(x, _p, _s, alpha, omega, dom);
        pbicgstab_4(r, s, t, omega, dom);
        _rho = rho;
        norm = calc_norm(r);
        it ++;
        printf("\r%6d %10.3e %10.3e", it, norm, rho);
        fflush(stdout);

        /* if (it % 1 == 0) {
            x.update_self();
            FILE *fo;
            char fname[128];
            sprintf(fname, "temperature.csv.%d", it);
            fo = fopen(fname, "w+t");
            if (fo) {
                fprintf(fo, "x,y,z,t\n");
                real_t x1, x2, x3, tp;
                for (int k = 0; k < dom.size[2]; k ++) {
                    for (int j = 0; j < dom.size[1]; j ++) {
                        for (int i = 0; i < dom.size[0]; i ++) {
                            int cur = Util::id(i,j,k,dom.size);
                            x1 = mesh.x.get(cur, 0);
                            x2 = mesh.x.get(cur, 1);
                            x3 = mesh.x.get(cur, 2);
                            tp = x.get(cur);
                            fprintf(fo, "%10.3e,%10.3e,%10.3e,%10.3e\n", x1, x2, x3, tp);
                        }
                    }
                }
            }
        } */
    } while (norm > 1e-9);
    printf("\n");
    _r.off_device();
    p.off_device();
    q.off_device();
    s.off_device();
    _p.off_device();
    _s.off_device();
    t.off_device();
}

#endif