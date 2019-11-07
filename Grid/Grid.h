#ifndef CPP_RZ_PIC_GRID_H
#define CPP_RZ_PIC_GRID_H

#include <iostream>

typedef double type_double; // c++0x

struct Grid {
    size_t Nz, Nr;
    type_double dz, dr;
    Grid() : Nz(0), Nr(0), dz(0), dr(0) {}
    Grid(size_t _Nz, size_t _Nr, type_double _dz, type_double  _dr) : Nz(_Nz), Nr(_Nr), dz(_dz), dr(_dr) {}
};


#endif //CPP_RZ_PIC_GRID_H
