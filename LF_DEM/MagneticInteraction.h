//
//  MagneticInteraction.h
//  LF_DEM
//
//  Created by Zhiyuan Zhao on 2020/9/1.
//  Copyright © 2020 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__MagneticInteraction__
#define __LF_DEM__MagneticInteraction__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "PotentialForce.h"
#include "MagneticInteractionParams.h"
#include "System.h"

class system;

namespace Interactions
{

class MagneticInteraction : public PotentialForce {
private:
    std::shared_ptr<System> sys;
    struct MagneticInteractionParams p;
    double geometric_factor;
    double reduced_force_norm;
    void MagneticDipoleMoment();
    vec3d orient_p0;                                    // Orientation vector of the magnetic dipole moment of particle 0
    vec3d orient_p1;                                    // Orientation vector of the magnetic dipole moment of particle 1
public:
    MagneticInteraction(PairwiseInteraction* interaction_, struct MagneticInteractionParams params);
    void calcForce();
    double calcEnergy() const;
};

}

#endif /* __LF_DEM__MagneticInteraction__ */
