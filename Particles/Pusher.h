//
// Created by Vladimir Smirnov on 05.11.2019.
//

#ifndef CPP_RZ_PIC_PUSHER_H
#define CPP_RZ_PIC_PUSHER_H

using namespace std;

typedef double type_double; // c++0x

void CrossProduct(const type_double v1[], const type_double v2[], type_double result[]);

void RotateToRZ(type_double& pos_r, type_double& pos_y, type_double& vel_r, type_double& vel_y);

void UpdateSingleVelocityBoris(type_double& vel_z, type_double& vel_r, type_double& vel_y, type_double Ez,
                               type_double Er, type_double Bz, type_double Br, type_double dt, type_double q,
                               type_double m);

void UpdateVelocity(type_double vel_z[], type_double vel_r[], type_double vel_y[], const type_double Ez[],
                    const type_double Er[], const type_double Bz[], const type_double Br[],
                    const type_double dt, const type_double q, const type_double m, const size_t Ntot);

void UpdatePosition(type_double pos_z[], type_double pos_r[], type_double vel_z[], type_double vel_r[],
                    type_double vel_y[], const type_double dt, const size_t Ntot, const type_double dr);

void ParticlePush(type_double pos_z[], type_double pos_r[], type_double vel_z[], type_double vel_r[],
                  type_double vel_y[], const type_double Ez[], const type_double Er[], const type_double Bz[],
                  const type_double Br[], const type_double dt, const type_double q, const type_double m,
                  const size_t Ntot, const type_double dr);


#endif //CPP_RZ_PIC_PUSHER_H
