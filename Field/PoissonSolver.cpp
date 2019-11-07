#include "../Tools/Matrix.h"
#include <cmath>
#define epsilon_0 8.854187817620389e-12

void compute_b(Matrix& b, Matrix& rho) {
    for (int i = 0; i < rho.rows(); i++) {
        for (int j = 0; j < rho.columns(); j++) {
            b(i, j) = rho(i, j) / epsilon_0;
        }
    }
}

void set_radii(Matrix& radii, int Nz, int Nr, double dr) {
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nr; j++) {
            radii(i, j) = j*dr;
        }
    }
}

bool convergence(Matrix& phi, Matrix& b, Matrix& r, const double dz, const double dr,
        const double tolerance) {
    double res;
    double dz2 = dz*dz;
    double dr2 = dr*dr;
    for (int i = 1; i < phi.rows()-1; i++) {
        for (int j = 1; j < phi.columns()-1; j++) {
            res = phi(i, j) - (b(i, j) + (phi(i, j-1) + phi(i, j+1))/dr2 +
                               (phi(i, j+1) - phi(i, j-1))/(2*dr*r(i, j)) +
                               (phi(i-1, j) + phi(i+1, j))/dz2) / (2/dr2 + 2/dz2);
            if (abs(res) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

void phi_init(Matrix& phi, const int CathodeR, const double CathodeV, const double AnodeV, double value=0) {
    int Nz = phi.rows(), Nr = phi.columns();
    for(int i = 0; i < Nz; i++) {
        for(int j = 0; j < Nr; j++) {
            if ((i == Nz-1 || i == 0) && j < CathodeR) {
                phi(i, j) = CathodeV;
            }
            else if (j == Nr-1) {
                phi(i, j) = AnodeV;
            }
            else {
                phi(i, j) = value;
            }
        }
    }
}

void PoissonSolverJacobi(Matrix& phi, Matrix& rho, Matrix& radii, const double dz, const double dr, const int CathodeR,
        const double tolerance, const int max_iter=1e6, int convergence_check=10) {
    int Nz = phi.rows();
    int Nr = phi.columns();
    Matrix g(Nz, Nr), b(Nz, Nr);
    double dz2 = dz*dz;
    double dr2 = dr*dr;
    compute_b(b, rho);
    g.copy(phi);
    for (int it = 0; it < max_iter; it++) {
        for (int i = 1; i < Nz-1; i++) {
            for (int j = 1; j < Nr-1; j++) {
                g(i, j) = (b(i, j) + (phi(i, j+1) + phi(i, j-1))/dr2 +
                                     (phi(i+1, j) + phi(i-1, j))/dz2 +
                                     (phi(i, j+1) - phi(i, j-1))/(radii(i, j)*2*dr))/(2/dr2 + 2/dz2);
            }
        }
        // Neumann boundaries on the axis r = 0
        for(int i = 0; i < Nz; i++) {
            g(i, 0) = g(i, 1);
        }
        // Neumann boundaries on the left/right wall
        for (int i = CathodeR; i < Nr; i++) {
            g(0, i) = g(1, i);
            g(Nz-1, i) = g(Nz-2, i);
        }
        // Dirichlet boundaries
        phi.copy(g);
        // Convergence check
        if (it % convergence_check == 0 and convergence(phi, b, radii, dz, dr, tolerance)) {
            cout << "Convergence achieved at iteration: " << it << endl;
            break;
        }
    }
}

void PoissonSolverSOR(Matrix& phi, Matrix& rho, Matrix& radii, const double dz, const double dr, const int CathodeR,
        const double tolerance, const int max_iter=1e6, const double betta=1.5, int convergence_check=10) {
    int Nz = phi.rows();
    int Nr = phi.columns();
    Matrix g(Nz, Nr), b(Nz, Nr);
    double dz2 = dz*dz;
    double dr2 = dr*dr;
    compute_b(b, rho);
    g.copy(phi);
    for (int it = 0; it < max_iter; it++) {
        for (int i = 1; i < Nz-1; i++) {
            for (int j = 1; j < Nr-1; j++) {
                g(i, j) = betta*((b(i, j) + (phi(i, j+1) + g(i, j-1))/dr2 + (phi(i+1, j) +
                        g(i-1, j))/dz2 + (phi(i, j+1) - g(i, j-1))/(radii(i, j)*2*dr)) /
                          (2/dr2 + 2/dz2)) + (1-betta)*phi(i, j);
            }
        }
        // Neumann boundaries on the axis r = 0
        //cout << g(1, 1) << endl;
        for(int i = 0; i < Nz; i++) {
            g(i, 0) = g(i, 1);
        }
        // Neumann boundaries on the left/right wall
        for (int i = CathodeR; i < Nr; i++) {
            g(0, i) = g(1, i);
            g(Nz-1, i) = g(Nz-2, i);
        }
        // Dirichlet boundaries
        phi.copy(g);
        // Convergence check
        if (it % convergence_check == 0 and convergence(phi, b, radii, dz, dr, tolerance)) {
            cout << "Convergence achieved at iteration: " << it << endl;
            break;
        }
    }
}