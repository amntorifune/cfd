#ifndef _DOM_H
#define _DOM_H 1

struct Dom {
    int   size[3];
    int origin[3];
    int    num;
    Dom(const int cx, const int cy, const int cz, const int ox, const int oy, const int oz);
    // void set(const int cx, const int cy, const int cz, const int ox, const int oy, const int oz);
    void to_device();
    void off_device();
};

Dom::Dom(const int x, const int y, const int z, const int ox, const int oy, const int oz) {
    size[0]   =  x;
    size[1]   =  y;
    size[2]   =  z;
    origin[0] = ox;
    origin[1] = oy;
    origin[2] = oz;
    num = x * y * z;
}

// void Dom::set(const int x, const int y, const int z, const int ox, const int oy, const int oz) {
//     size[0]   =  x;
//     size[1]   =  y;
//     size[2]   =  z;
//     origin[0] = ox;
//     origin[1] = oy;
//     origin[2] = oz;
//     num = x * y * z;
// }

void Dom::to_device() {
    #pragma acc enter data copyin(this[0:1])
}

void Dom::off_device() {
    #pragma acc exit data delete(this[0:1])
}

#endif