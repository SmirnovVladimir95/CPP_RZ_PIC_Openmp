//
// Created by Vladimir Smirnov on 22.10.2019.
//

#ifndef CPP_RZ_PIC_FIELD_H
#define CPP_RZ_PIC_FIELD_H
#include "../Tools/Matrix.h"
typedef double my_double; // c++0x

class Field {
private:
    my_double _dz, _dr;
    int _Nz, _Nr;
    Matrix _cell_type, _r;
    Matrix set_radia(my_double dr);
    bool convergence();

public:
    void PoissonSolverJacobi(Matrix phi, Matrix rho);
    void PoissonSolverSOR(Matrix phi, Matrix rho);
};


#endif //CPP_RZ_PIC_FIELD_H
