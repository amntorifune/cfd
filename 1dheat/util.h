#ifndef _UTIL_H
#define _UTIL_H 1

namespace Util {

static inline int id(const int i, const int j, const int k, const int* size) {
    return i * (size[1] * size[2]) + j * size[2] + k;
}

template<typename T>
static T max(T a, T b) {
    return (a > b)? (a) : (b);
}

template<typename T>
static T min(T a, T b) {
    return (a < b)? (a) : (b);
}

static inline int eucmod(int a, int b) {
    return (a < 0)? (((a % b) + b) % b) : (a % b);
}

static unsigned ibsee(unsigned bits, unsigned i, unsigned mask) {
    return (bits >> i) & mask;
}

static unsigned ibset(unsigned bits, unsigned i, unsigned mask, unsigned value) {
    return (bits & ~(mask << i)) | (value << i);
}

}
#endif