#ifndef _MATRIX_H
#define _MATRIX_H 1

#include <cstring>
#include <cmath>
#include <iostream>
#include "type.h"
#include "dom.h"
#include "util.h"

template<class T>
class Matrix {
public:
    T  *mat;
    int row;
    int col;
    int num;
public:
    Matrix(int row, int col);
    Matrix(Dom &dom, int col);
    T &get(int i, int j);
    T &get(int i);
    void to_device();
    void off_device();
    void update_device();
    void update_self();
    Matrix<T> &operator=(const Matrix<T> &b);
    Matrix<T> &operator=(const T &b);
    ~Matrix();
};

template<class T>
Matrix<T>::Matrix(int row, int col) : row(row), col(col) {
    num = row * col;
    mat = new T[num];
    memset(mat, 0, num * sizeof(T));
}

template<class T>
Matrix<T>::Matrix(Dom &dom, int col) : row(dom.num), col(col) {
    num = row * col;
    mat = new T[num];
    memset(mat, 0, num * sizeof(T));
}

template<class T>
inline T &Matrix<T>::get(int i, int j) {
    return mat[i * col + j];
}

template<class T>
inline T &Matrix<T>::get(int i) {
    return mat[i];
}

template<class T>
void Matrix<T>::to_device() {
    #pragma acc enter data copyin(this[0:1], mat[0:num])
}

template<class T>
void Matrix<T>::off_device() {
    #pragma acc exit data delete(mat[0:num], this[0:1])
}

template<class T>
void Matrix<T>::update_device() {
    #pragma acc update device(mat[0:num])
}

template<class T>
void Matrix<T>::update_self() {
    #pragma acc update self(mat[0:num])
}

template<class T>
Matrix<T>::~Matrix() {
    delete[] mat;
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &b) {
    #pragma acc kernels loop independent present(this[0:1], b)
    for (int i = 0; i < num; i ++) {
        mat[i] = b.mat[i];
    }
    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const T &b) {
    #pragma acc kernels loop independent present(this[0:1]) copyin(b)
    for (int i = 0; i < num; i ++) {
        mat[i] = b;
    }
    return *this;
}

#endif