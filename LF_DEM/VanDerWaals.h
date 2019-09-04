//
//  RepulsiveForce.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/11/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class RepulsiveForce
 \brief Exponentially decaying repulsion, to be called from an Interaction object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__vanDerWaalsForce__
#define __LF_DEM__vanDerWaalsForce__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "PotentialForce.h"
#include "VanDerWaalsParams.h"

namespace Interactions
{

class vanDerWaalsForce : public PotentialForce {
private:
	struct vanDerWaalsForceParams p;
	double geometric_factor;
	double reduced_force_norm;
	void calcReducedForceNorm();
	void calcScaledForce();

public:
	vanDerWaalsForce(PairwiseInteraction* interaction_, struct vanDerWaalsForceParams params);
	void calcForce();
	double calcEnergy() const;
};

} // namespace Interactions

#endif /* defined(__LF_DEM__vanDerWaalsForce__) */
