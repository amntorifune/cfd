#ifndef _MESH_H
#define _MESH_H

#include "matrix.h"

class Mesh {
public:
    Matrix<real_t>   x;
    Matrix<real_t>   h;
    Matrix<real_t>   v;
    Matrix<real_t>   map;
    Matrix<unsigned> flag;
public:
    Mesh(Dom &dom);
    void to_device();
    void off_device();
};

Mesh::Mesh(Dom &dom) : x(dom, 3), h(dom, 3), v(dom, 1), map(dom, 6), flag(dom, 1) {}

void Mesh::to_device() {
    #pragma acc enter data copyin(this[0:1])
    x.to_device();
    h.to_device();
    v.to_device();
    map.to_device();
}

void Mesh::off_device() {
    x.off_device();
    h.off_device();
    v.off_device();
    map.off_device();
    #pragma acc exit data delete(this[0:1])
}

#endif