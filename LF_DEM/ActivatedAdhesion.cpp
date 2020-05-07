//
//  ActivatedAdhesion.cpp
//  LF_DEM
//
//  Created by Romain Mari on 14/02/20.
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//

#include <sstream>
#include "ActivatedAdhesion.h"

namespace Interactions {

namespace ActAdhesion {

ActivatedAdhesion::ActivatedAdhesion(PairwiseInteraction* interaction_, struct Params &p):
interaction(interaction_),
params(p),
state({Activity::dormant, interaction_->p0, interaction_->p1}),
force_amplitude(0),
force_on_p0(0),
stress_split_p0(interaction_->a0/(interaction_->a0+interaction_->a1)) {
}

MTRand *ActivatedAdhesion::r_gen = new MTRand;

void ActivatedAdhesion::update(double dt)
{
	if (state.activity == Activity::dormant) {
		double prob = params.activation_rate*dt;
		if (prob > 0.2) {
			std::ostringstream error_str;
			error_str << " ActivatedAdhesion:: Time step too large! Activation rate = " << params.activation_rate << ", dt = " << dt;
			std::runtime_error(error_str.str());
		}
		if (r_gen->rand() < prob) {
			state.activity = Activity::active;
		}
	} else {
		double prob = params.deactivation_rate*dt;
		if (prob > 0.2) {
			std::ostringstream error_str;
			error_str << " ActivatedAdhesion:: Time step too large! Deactivation rate = " << params.deactivation_rate << ", dt = " << dt << std::endl;
			std::runtime_error(error_str.str());
		}
		if (r_gen->rand() < prob) {
			state.activity = Activity::dormant;
		}
	}
	if (state.activity == Activity::active) {
		if (interaction->reduced_gap > 0) {// we also know gap < range 
			force_amplitude = params.max_force*interaction->reduced_gap/params.range; 
		} else {
			force_amplitude = 0;
		}
		force_on_p0 = force_amplitude*interaction->nvec;
	}
}

void ActivatedAdhesion::addUpForce(vec3d &force_p0, vec3d &force_p1) const 
{
	if (state.activity == Activity::active) {
		force_p0 += force_on_p0;
		force_p1 -= force_on_p0;
	}
}

void ActivatedAdhesion::addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1, const vec3d &rvec) const
{
	if (state.activity == Activity::active) {
		auto sc = outer_sym(rvec, force_on_p0);
		stress_p0 += stress_split_p0*sc;
		stress_p1 += (1-stress_split_p0)*sc;
	}
}

void ActivatedAdhesion::setState(struct State st)
{
	state = st;
}

double ActivatedAdhesion::getForceNorm() const
{
	return force_amplitude;
}

} // namespace ActAdhesion

}