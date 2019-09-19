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

#ifndef __LF_DEM__RepulsiveForce__
#define __LF_DEM__RepulsiveForce__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "PotentialForce.h"
#include "RepulsiveForceParams.h"

namespace Interactions
{

class RepulsiveForce : public PotentialForce {
private:
	struct RepulsiveForceParams p;
	//===== forces and stresses ==================== //
	double geometric_factor;
	double reduced_force_norm;
	void calcReducedForceNorm();
	void calcScaledForce();
public:
	RepulsiveForce(PairwiseInteraction* interaction_, struct RepulsiveForceParams params);
	void calcForce();
	double calcEnergy() const;
};

inline double calcRepulsiveForceMaxGap(RepulsiveForceParams input_p)
{
	double max_gap;
	if (input_p.max_gap == -1) {
		max_gap = 7*input_p.screening_length;
	} else {
		max_gap = input_p.max_gap;
	}
	return max_gap;
}

inline double calcRepulsiveForceRange(RepulsiveForceParams input_p, double a0, double a1)
{
	auto max_gap = calcRepulsiveForceMaxGap(input_p);
	return (a0 + a1)*(1 + 0.5*max_gap);
}

} // namespace Interactions

#endif /* defined(__LF_DEM__RepulsiveForce__) */
