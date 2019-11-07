#ifndef CPP_RZ_PIC_POISSONSOLVER_H
#define CPP_RZ_PIC_POISSONSOLVER_H

#include "../Tools/Matrix.h"

void set_radii(Matrix& radii, int Nz, int Nr, double dr);
void compute_b(Matrix& b, Matrix& rho);
void phi_init(Matrix& phi, const int CathodeR, const double CathodeV, const double AnodeV, double value=0);
bool convergence(Matrix& phi, Matrix& b, Matrix& r, const double dz, const double dr, const double tolerance);
void PoissonSolverJacobi(Matrix& phi, Matrix& rho, Matrix& radii, const double dz, const double dr, const int CathodeR,
                         const double tolerance, const int max_iter=1e6, int convergence_check=10);
void PoissonSolverSOR(Matrix& phi, Matrix& rho, Matrix& radii, const double dz, const double dr, const int CathodeR,
                      const double tolerance, const int max_iter=1e6, const double betta=1.5, int convergence_check=10);

#endif //CPP_RZ_PIC_POISSONSOLVER_H
