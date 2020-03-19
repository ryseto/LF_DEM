//
//  ActivatedAdhesion.h
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class TimeActivatedAdhesion
 \brief Adhesive force activated only after a given contact time.
 \author Romain Mari
 */

#ifndef __LF_DEM__ActivatedAdhesion__
#define __LF_DEM__ActivatedAdhesion__

#include "vec3d.h"
#include "Sym2Tensor.h"
#include "PairwiseInteraction.h"
#include "ActivatedAdhesion_Params.h"
#include "ActivatedAdhesionState.h"
#include "MersenneTwister.h"

namespace Interactions {

namespace ActAdhesion {

class ActivatedAdhesion {

public:
	ActivatedAdhesion(PairwiseInteraction* interaction_, struct Params &p);
	
	void update(double dt);
	void deactivate();
	void setState(struct State st);
	struct State getState() const {return state;};
	void addUpForce(vec3d &force_p0, vec3d &force_p1) const;
	void addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1, const vec3d &rvec) const;
	double ratioUptimeToActivation() const;
	double getForceNorm() const;

private:

	PairwiseInteraction* interaction;
	struct Params params;
	struct State state;

	double force_amplitude;
	vec3d force_on_p0;
	double stress_split_p0;
	double initial_time;

	static MTRand *r_gen;
};

} // namespace ActAdhesion

}

#endif // ifndef __LF_DEM__ActivatedAdhesion__
