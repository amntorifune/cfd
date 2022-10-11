#ifndef _TYPE_H
#define _TYPE_H 1

typedef double real_t;

enum class stencil_t {d1s3, d1s5, d2s5, d2s9, d3s7, d3s13};
enum class sgs_t {smagorinsky, csm};

#endif