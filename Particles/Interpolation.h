//
// Created by Vladimir Smirnov on 06.11.2019.
//

#ifndef CPP_RZ_PIC_INTERPOLATION_H
#define CPP_RZ_PIC_INTERPOLATION_H

#include "../Tools/Matrix.h"
#include "../Grid/Grid.h"

void LinearFieldInterpolation(type_double efz[], type_double efr[], const type_double z[], const type_double r[],
                              const type_double Ez[], const type_double Er[], const Grid& grid, const size_t Ntot);

void LinearChargeInterpolation(type_double rho[], const type_double z[], const type_double r[], const Grid& grid,
                               type_double charge, const size_t Ntot, const type_double node_volume[]);


#endif //CPP_RZ_PIC_INTERPOLATION_H
