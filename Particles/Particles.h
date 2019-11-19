#ifndef CPP_RZ_PIC_PARTICLES_H
#define CPP_RZ_PIC_PARTICLES_H

#include <cstddef>
#include <vector>
#include <array>
#include <random>
#include "../Tools/Matrix.h"
#include "../Grid/Grid.h"
using namespace std;

typedef double type_double; // c++0x

class Particles {
private:
    Grid grid;
    void init_node_volume(Matrix node_volume);
public:
    size_t Ntot;
    type_double mass;
    type_double charge;
    vector<type_double> z;
    vector<type_double> r;
    vector<type_double> vz;
    vector<type_double> vr;
    vector<type_double> vy;
    vector<type_double> efz;
    vector<type_double> efr;
    vector<type_double> mfz;
    vector<type_double> mfr;
    Matrix rho;
    Matrix node_volume;
    Particles(const type_double m, const type_double q, const size_t N, const Grid& grid);
    void generate_velocities(const type_double energy, const int seed=time(nullptr));
    void generate_positions(const array<type_double, 2>& z_bounds, const array<type_double, 2>& r_bounds, const int seed=time(nullptr));
    vector<vector<type_double>> get_positions() const;
    vector<vector<type_double>> get_velocities() const;
    void pusher(const type_double  dt);
    void vel_pusher(const type_double  dt);
    void electric_field_interpolation(Matrix& Ez, Matrix& Er);
    void magnetic_field_interpolation(Matrix& Bz, Matrix& Br);
    void set_const_magnetic_field(const vector<type_double>& Bz, const vector<type_double>& Br);
    void charge_interpolation();
};


#endif //CPP_RZ_PIC_PARTICLES_H