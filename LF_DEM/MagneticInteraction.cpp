//
//  MagneticInteraction.cpp
//  LF_DEM
//
//  Created by Zhiyuan Zhao on 2020/9/1.
//  Copyright © 2020 Ryohei Seto. All rights reserved.
//

#include <vector>
#include "MagneticInteraction.h"

namespace Interactions
{

MagneticInteraction::MagneticInteraction(PairwiseInteraction* interaction_, struct MagneticInteractionParams params) :
PotentialForce(interaction_),
p(params)
{
    calcForce();
}

void MagneticInteraction::calcForce()
{
    double gap = interaction->getGap();
    vec3d e_force_dd = interaction->nvec*dot(interaction->dipole_moment_p0,interaction->dipole_moment_p1) + interaction->dipole_moment_p0*dot(interaction->nvec,interaction->dipole_moment_p1) + interaction->dipole_moment_p1*dot(interaction->nvec,interaction->dipole_moment_p0) - 5*interaction->nvec*dot(interaction->nvec,interaction->dipole_moment_p1)*dot(interaction->nvec,interaction->dipole_moment_p0);
    vec3d e_torque_dd = dot(interaction->dipole_moment_p1,interaction->nvec)*cross(interaction->dipole_moment_p0,interaction->nvec) + cross(interaction->dipole_moment_p1,interaction->dipole_moment_p0)/3;
    if (gap > 0) {
        force_vector = p.beta_dipoledipole*(-e_force_dd)/pow(interaction->r, 4);
        torque_vector = p.beta_dipoledipole*e_torque_dd/pow(interaction->r, 3);
    } else {
        force_vector = p.beta_dipoledipole*(-e_force_dd)/pow(interaction->a0+interaction->a1, 4);
        torque_vector = p.beta_dipoledipole*e_torque_dd/pow(interaction->a0+interaction->a1, 3);
    }
}

double MagneticInteraction::calcEnergy() const
{
    double energy;
    double gap = interaction->getGap();
    if (gap > 0) {
        energy = p.beta_dipoledipole*(3*dot(interaction->dipole_moment_p0,interaction->nvec)*dot(interaction->dipole_moment_p1,interaction->nvec) - dot(interaction->dipole_moment_p0,interaction->dipole_moment_p1))/(4*pow(interaction->r,3));
    } else {
        energy = p.beta_dipoledipole*(3*dot(interaction->dipole_moment_p0,interaction->nvec)*dot(interaction->dipole_moment_p1,interaction->nvec) - dot(interaction->dipole_moment_p0,interaction->dipole_moment_p1))/(4*pow(interaction->a0+interaction->a1,3));
    }
    return energy;
}

}
