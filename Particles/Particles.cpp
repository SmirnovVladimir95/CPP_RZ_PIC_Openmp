#include "Particles.h"
#include <random>
#include <cmath>
#include "../Tools/Matrix.h"
#include "Pusher.h"
#include "Interpolation.h"

Particles::Particles(const type_double m, const type_double q, const size_t N, const Grid& init_grid) {
    mass = m;
    charge = q;
    z.resize(N, 0);
    r.resize(N, 0);
    vz.resize(N, 0);
    vr.resize(N, 0);
    vy.resize(N, 0);
    efz.resize(N, 0);
    efr.resize(N, 0);
    Ntot = N;
    grid = init_grid;
    node_volume.resize(init_grid.Nz, init_grid.Nr);
    init_node_volume(node_volume);
    rho.resize(init_grid.Nz, init_grid.Nr);
}

void Particles::generate_velocities(const type_double energy, const int seed) {
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0, sqrt(2*energy/(3*mass)));
    for (int i = 0; i < Ntot; i++) {
        vz[i] = distribution(generator);
        vr[i] = distribution(generator);
        vy[i] = distribution(generator);
    }
}

void Particles::generate_positions(const array<type_double, 2> &z_bounds, const array<type_double, 2> &r_bounds,
                                   const int seed) {
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution_z(z_bounds[0],z_bounds[1]);
    std::uniform_real_distribution<double> distribution_r(r_bounds[0],r_bounds[1]);
    for (int i = 0; i < Ntot; i++) {
        z[i] = distribution_z(generator);
        r[i] = distribution_r(generator);
    }
}

void Particles::vel_pusher(const type_double  dt) {
    auto vel_z = vz.data();
    auto vel_r = vr.data();
    auto vel_y = vy.data();
    auto Ez = efz.data();
    auto Er = efr.data();
    auto Bz = mfz.data();
    auto Br = mfr.data();
    UpdateVelocity(vel_z, vel_r, vel_y, Ez, Er, Bz, Br, dt, charge, mass, Ntot);
}

void Particles::pusher(const type_double dt) {
    auto pos_z = z.data();
    auto pos_r = r.data();
    auto vel_z = vz.data();
    auto vel_r = vr.data();
    auto vel_y = vy.data();
    auto Ez = efz.data();
    auto Er = efr.data();
    auto Bz = mfz.data();
    auto Br = mfr.data();
    ParticlePush(pos_z, pos_r, vel_z, vel_r, vel_y, Ez, Er, Bz, Br, dt, charge, mass,
                 Ntot, grid.dr);
    //UpdateVelocity(vel_z, vel_r, vel_y, Ez, Er, Bz, Br, dt, charge, mass, Ntot);
    //UpdatePosition(pos_z, pos_r, vel_z, vel_r, vel_y, dt, Ntot, grid.dr);
}

void Particles::init_node_volume(Matrix node_volume) {
    int j_min, j_max;
    float a;
    for (int i = 0; i < grid.Nz; i++) {
        for (int j = 0; j < grid.Nr; j++) {
            j_min = j - 0.5;
            j_max = j + 0.5;
            if (j_min < 0)
                j_min = 0;
            if (j_max > grid.Nr - 1)
                j_max = grid.Nr - 1;
            if (i == 0 or i == grid.Nz - 1)
                a = 0.5;
            else
                a = 1;
            node_volume(i, j) = a*grid.dz*((j_max*grid.dr)*(j_max*grid.dr) - (j_min*grid.dr)*(j_min*grid.dr))*M_PI;
        }
    }
}

void Particles::electric_field_interpolation(Matrix& Ez, Matrix& Er) {
    LinearFieldInterpolation(efz.data(), efr.data(), z.data(), r.data(), Ez.data_ptr(), Er.data_ptr(),
                             grid, Ntot);
}

void Particles::charge_interpolation() {
    LinearChargeInterpolation(rho.data_ptr(), z.data(), r.data(), grid, charge, Ntot, node_volume.data_ptr());
}

void Particles::set_const_magnetic_field(const vector<type_double>& Bz, const vector<type_double>& Br) {
    mfz = Bz;
    mfr = Br;
}