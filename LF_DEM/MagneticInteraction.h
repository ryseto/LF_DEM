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
#include <memory>
#include <vector>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "PotentialForce.h"
#include "ParticleConfig.h"
#include "MagneticInteractionParams.h"

namespace Interactions
{

class MagneticInteraction : public PotentialForce {
private:
    struct MagneticInteractionParams p;
    std::shared_ptr<ParticleConfig> conf;
    void calcDipoleOrientation();
    vec3d dipole_orient_p0;
    vec3d dipole_orient_p1;
public:
    MagneticInteraction(PairwiseInteraction* interaction_, struct MagneticInteractionParams params);
//    void calcDipoleOrientation(vec3d orient_p0, vec3d orient_p1);
    void calcForce();
    double calcEnergy() const;
};

}

#endif /* __LF_DEM__MagneticInteraction__ */
