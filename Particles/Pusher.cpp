#include <cmath>
#include <omp.h>
#define NUM_THREADS 100
#include "Pusher.h"

void CrossProduct(const type_double v1[], const type_double v2[], type_double result[]) {
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - 1*v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void RotateToRZ(type_double& pos_r, type_double& pos_y, type_double& vel_r, type_double& vel_y, type_double dr) {
    type_double r, sin_theta, cos_theta;
    r = sqrt(pos_r*pos_r + pos_y*pos_y);
    pos_r = r;
    if (r > dr) {
        sin_theta = pos_y / r;
        cos_theta = sqrt(1-sin_theta*sin_theta);
    } else {
        sin_theta = 0;
        cos_theta = 1;
    }
    vel_r = cos_theta*vel_r - sin_theta*vel_y;
    vel_y = sin_theta*vel_r + cos_theta*vel_y;
}

void UpdateSingleVelocityBoris(type_double& vel_z, type_double& vel_r, type_double& vel_y, type_double Ez,
                               type_double Er, type_double Bz, type_double Br, type_double dt, type_double q,
                               type_double m) {
    int dim = 3;
    type_double t[dim], s[dim], v_minus[dim], v_minus_cross[dim], v_prime[dim], v_prime_cross[dim], v_plus[dim];
    type_double q_div_m = q/m;
    t[0] = q_div_m*Bz*0.5*dt;
    t[1] = q_div_m*Br*0.5*dt;
    t[2] = 0;
    double t_mag2 = t[0]*t[0] + t[1]*t[1];
    s[0] = 2*t[0]/(1+t_mag2);
    s[1] = 2*t[1]/(1+t_mag2);
    s[2] = 0;
    v_minus[0] = vel_z + q_div_m*Ez*0.5*dt;
    v_minus[1] = vel_r + q_div_m*Er*0.5*dt;
    v_minus[2] = vel_y;
    CrossProduct(v_minus, t, v_minus_cross);
    for (int i = 0; i < dim; i++) {
        v_prime[i] = v_minus[i] + v_minus_cross[i];
    }
    CrossProduct(v_prime, s, v_prime_cross);
    for (int i = 0; i < dim; i++) {
        v_plus[i] = v_minus[i] + v_prime_cross[i];
    }
    vel_z = v_plus[0] + q_div_m*Ez*0.5*dt;
    vel_r = v_plus[1] + q_div_m*Er*0.5*dt;
    vel_y = v_plus[2];
}

void UpdateVelocity(type_double vel_z[], type_double vel_r[], type_double vel_y[], const type_double Ez[],
                    const type_double Er[], const type_double Bz[], const type_double Br[],
                    const type_double dt, const type_double q, const type_double m, const size_t Ntot) {
    // Loop over ptcls
    #pragma omp for
    for (int ip = 0; ip < Ntot; ip++) {
        UpdateSingleVelocityBoris(vel_z[ip], vel_r[ip], vel_y[ip], Ez[ip], Er[ip], Bz[ip], Br[ip], dt, q, m);
    }
}

void UpdatePosition(type_double pos_z[], type_double pos_r[], type_double vel_z[], type_double vel_r[],
                    type_double vel_y[], const type_double dt, const size_t Ntot, const type_double dr) {
    type_double pos_y_ip=0;
    // Loop over ptcls
    #pragma omp for
    for (int ip = 0; ip < Ntot; ip++) {
        pos_z[ip] += vel_z[ip]*dt;
        pos_r[ip] += vel_r[ip]*dt;
        pos_y_ip += vel_y[ip]*dt;
        RotateToRZ(pos_r[ip], pos_y_ip, vel_r[ip], vel_y[ip], dr);
    }
}

void ParticlePush(type_double pos_z[], type_double pos_r[], type_double vel_z[], type_double vel_r[],
                  type_double vel_y[], const type_double Ez[], const type_double Er[], const type_double Bz[],
                  const type_double Br[], const type_double dt, const type_double q, const type_double m,
                  const size_t Ntot, const type_double dr) {
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        UpdateVelocity(vel_z, vel_r, vel_y, Ez, Er, Bz, Br, dt, q, m, Ntot);
        UpdatePosition(pos_z, pos_r, vel_z, vel_r, vel_y, dt, Ntot, dr);
    }
}
