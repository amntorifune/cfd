#ifndef _MESH_H
#define _MESH_H

#include "matrix.h"

class Mesh {
public:
    Mat<real_t> x;
    Mat<real_t> h;
    Mat<real_t> v;
public:
    Mesh(Dom &dom);
    void to_device();
    void off_device();
};

Mesh::Mesh(Dom &dom) : x(dom.num, 3), h(dom.num, 3), v(dom.num, 1) {}

void Mesh::to_device() {
    #pragma acc enter data copyin(this[0:1])
    x.to_device();
    h.to_device();
    v.to_device();
}

void Mesh::off_device() {
    x.off_device();
    h.off_device();
    v.off_device();
    #pragma acc exit data delete(this[0:1])
}

#endif