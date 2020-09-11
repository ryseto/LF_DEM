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
#include "ParticleConfig.h"
#include "MagneticInteractionParams.h"

namespace Interactions
{

class MagneticInteraction : public PotentialForce {
private:
    struct MagneticInteractionParams p;
public:
    MagneticInteraction(PairwiseInteraction* interaction_, struct MagneticInteractionParams params);
    void calcForce();
    double calcEnergy() const;
};

}

#endif /* __LF_DEM__MagneticInteraction__ */
