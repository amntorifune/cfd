#include "poisson.h"

Dom dom(N, 1, 1, 0, 0, 0);

int main() {
    Mesh mesh(dom);
    real_t h = L / N;
    real_t x = 0;
    for (int i = 0; i < dom.num; i ++) {
        mesh.h.get(i, 0) = h;
        mesh.h.get(i, 1) = h;
        mesh.h.get(i, 2) = h;
        mesh.x.get(i, 0) = x + 0.5 * h;
        mesh.x.get(i, 1) = 0;
        mesh.x.get(i, 2) = 0;
        mesh.v.get(i)    = h * h * h;
        x += h;
    }
    Mat<real_t> t(dom, 1);
    Mat<real_t> a(dom, 3);
    Mat<real_t> rhs(dom, 1);
    Mat<real_t> res(dom, 1);
    real_t norm;

    mesh.to_device();
    dom.to_device();
    t.to_device();
    a.to_device();
    rhs.to_device();
    res.to_device();

    poisson_prepare_eq(a, rhs, mesh, dom, stencil_t::d1s3);
    poisson_sor(a, t, rhs, res, norm, mesh, dom);

    t.update_self();
    mesh.off_device();
    dom.off_device();
    t.off_device();
    a.off_device();
    rhs.off_device();
    res.off_device();

    for(int i = 0; i < dom.num; i ++) {
        printf("%10.3e ", t.get(i));
    }
    printf("\n");
}
