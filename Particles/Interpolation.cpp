//
// Created by Vladimir Smirnov on 06.11.2019.
//

#include "Interpolation.h"
#include "../Grid/Grid.h"
#include <cmath>
#include <omp.h>
#define NUM_THREADS 100

void LinearFieldInterpolation(type_double efz[], type_double efr[], const type_double z[], const type_double r[],
                              const type_double Ez[], const type_double Er[], const Grid& grid, const size_t Ntot) {
    int cell_z, cell_r, Nr=grid.Nr;
    type_double hz, hr;
    #pragma omp parallel for private(hz, hr, cell_z, cell_r) num_threads(NUM_THREADS)
    for (int i = 0; i < Ntot; i++) {
        cell_z = floor(z[i]/grid.dz);
        cell_r = floor(r[i]/grid.dr);
        hz = (z[i] - cell_z*grid.dz) / grid.dz;
        hr = (r[i] - cell_r*grid.dr) / grid.dr;

        efz[i] = Ez[cell_z*Nr + cell_r] * (1 - hz) * (1 - hr);
        efz[i] += Ez[(cell_z+1)*Nr + cell_r] * hz * (1 - hr);
        efz[i] += Ez[(cell_z+1)*Nr + cell_r+1] * hz * hr;
        efz[i] += Ez[cell_z*Nr + cell_r+1] * (1 - hz) * hr;

        efr[i] = Er[cell_z*Nr + cell_r] * (1 - hz) * (1 - hr);
        efr[i] += Er[(cell_z+1)*Nr + cell_r] * hz * (1 - hr);
        efr[i] += Er[(cell_z+1)*Nr + cell_r+1] * hz * hr;
        efr[i] += Er[cell_z*Nr + cell_r+1] * (1 - hz) * hr;
    }
}

void LinearChargeInterpolation(type_double rho[], const type_double z[], const type_double r[], const Grid& grid,
                               type_double charge, const size_t Ntot, const type_double node_volume[]) {
    int cell_z, cell_r, Nr = grid.Nr;
    type_double hz, hr;
    #pragma omp parallel private(hz, hr, cell_z, cell_r) num_threads(NUM_THREADS)
    {
    #pragma omp for
    for (int i = 0; i < grid.Nz*grid.Nr; i++) {
        rho[i] = 0;
    }
    #pragma omp for
    for (int i = 0; i < Ntot; i++) {
        cell_z = floor(z[i] / grid.dz);
        cell_r = floor(r[i] / grid.dr);
        hz = (z[i] - cell_z * grid.dz) / grid.dz;
        hr = (r[i] - cell_r * grid.dr) / grid.dr;

        rho[cell_z * Nr + cell_r] += charge * (1 - hz) * (1 - hr) / node_volume[cell_z * Nr + cell_r];
        rho[(cell_z + 1) * Nr + cell_r] += charge * hz * (1 - hr) / node_volume[(cell_z + 1) * Nr + cell_r];
        rho[(cell_z + 1) * Nr + cell_r + 1] += charge * hz * hr / node_volume[(cell_z + 1) * Nr + cell_r + 1];
        rho[cell_z * Nr + cell_r + 1] += charge * (1 - hz) * hr / node_volume[cell_z * Nr + cell_r + 1];
    }
    }
}
