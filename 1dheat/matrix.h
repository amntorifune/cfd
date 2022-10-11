#ifndef _MATRIX_H
#define _MATRIX_H 1

#include <string.h>
#include <math.h>
#include "type.h"
#include "dom.h"
#include "util.h"

template<class T>
class Mat {
public:
    T  *mat;
    int row;
    int col;
    int num;
public:
    Mat(int row, int col);
    Mat(Dom &dom, int col);
    T &get(int i, int j);
    T &get(int i);
    void to_device();
    void off_device();
    void update_device();
    void update_self();
    ~Mat();
};

template<class T>
Mat<T>::Mat(int row, int col) : row(row), col(col) {
    num = row * col;
    mat = new T[num];
    memset(mat, 0, num * sizeof(T));
}

template<class T>
Mat<T>::Mat(Dom &dom, int col) : row(dom.num), col(col) {
    num = row * col;
    mat = new T[num];
    memset(mat, 0, num * sizeof(T));
}

template<class T>
inline T &Mat<T>::get(int i, int j) {
    return mat[i * col + j];
}

template<class T>
inline T &Mat<T>::get(int i) {
    return mat[i];
}

template<class T>
void Mat<T>::to_device() {
    #pragma acc enter data copyin(this[0:1], mat[0:num])
}

template<class T>
void Mat<T>::off_device() {
    #pragma acc exit data delete(mat[0:num], this[0:1])
}

template<class T>
void Mat<T>::update_device() {
    #pragma acc update device(mat[0:num])
}

template<class T>
void Mat<T>::update_self() {
    #pragma acc update self(mat[0:num])
}

template<class T>
Mat<T>::~Mat() {
    delete[] mat;
}

static void calc_res(Mat<real_t> &a, Mat<real_t> &x, Mat<real_t> &b, Mat<real_t> &r, real_t &norm, Dom &dom, stencil_t s) {
    real_t sum = 0;
    if (s == stencil_t::d1s3) {
        #pragma acc kernels loop independent collapse(2) reduction(+:sum) present(a, x, b, r, dom) copy(sum)
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
        }
    }
    norm = sqrt(sum / x.num);
}

#endif