//
//  TimeActivatedAdhesion.h
//  LF_DEM
//
//  Created by Romain Mari on 01/11/17.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class TimeActivatedAdhesion
 \brief Adhesive force activated only after a given contact time.
 \author Ryohei Seto
 \author Romain Mari
 */


#ifndef __LF_DEM__TimeActivatedAdhesion__
#define __LF_DEM__TimeActivatedAdhesion__

#include "vec3d.h"
#include "Sym2Tensor.h"
#include "TAAParams.h"


namespace TActAdhesion {

// explicit numbering as it is used in output file
enum class TAAActivity : unsigned {
	inactive = 0,
	dormant = 1,
	active = 2
};

struct TAAState {
	TAAActivity activity;
	double contacting_time;
};


class TimeActivatedAdhesion{

public:
	TimeActivatedAdhesion(struct TAAParams &p, double r0, double r1) 
	: params(p), 
	state({TAAActivity::inactive, 0}),
	force_amplitude(0),
	force_on_p0(0),
	stress_split_p0(r0/(r0+r1)) {};
	void update(double time_now, double gap, vec3d &nvec);
	void deactivate();

	void addUpForce(vec3d &force_p0, vec3d &force_p1) const;
	void addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1, const vec3d &rvec) const;

private:
	struct TAAParams params;
	TAAState state;
	double force_amplitude;
	vec3d force_on_p0;
	double stress_split_p0;

};

} // namespace TActAdhesion

#endif // ifndef __LF_DEM__TimeActivatedAdhesion__