#include "TimeActivatedAdhesion.h"

namespace TActAdhesion {

void TimeActivatedAdhesion::update(double time_now, double gap, vec3d &nvec)
{
	if (gap > params.adhesion_range) {
		if (state.activity != Activity::inactive) {
			state.activity = Activity::inactive;
			force_amplitude = 0;
		}
	} else {
		switch (state.activity) {
			case Activity::inactive:
				state.activity = Activity::dormant;
				initial_time = time_now;
				// no break, as we want to allow for immediate activation
			case Activity::dormant:
				state.uptime = time_now - initial_time;
				if (state.uptime >= params.activation_time) {
					state.activity = Activity::active;
				}
			default:
				break;

		}
		if (state.activity == Activity::active) {
			if (gap > 0) {// we also know gap < range 
				force_amplitude = params.adhesion_max_force*gap/params.adhesion_range; 
			} else {
				force_amplitude = 0;
			}
			force_on_p0 = -force_amplitude*nvec;
		}
	}
}

void TimeActivatedAdhesion::deactivate() 
{
	state.activity = Activity::inactive;
	force_on_p0 = 0;
	force_amplitude = 0;
}

void TimeActivatedAdhesion::addUpForce(vec3d &force_p0, vec3d &force_p1) const 
{
	if (state.activity == Activity::active) {
		force_p0 += force_on_p0;
		force_p1 -= force_on_p0;
	}
}

void TimeActivatedAdhesion::addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1, const vec3d &rvec) const
{
	auto sc = outer_sym(rvec, force_on_p0);
	stress_p0 += stress_split_p0*sc;
	stress_p1 += (1-stress_split_p0)*sc;
}

void TimeActivatedAdhesion::setState(struct State st, double time_now)
{
	state = st;
	initial_time = time_now - state.uptime;
}



} // namespace TActAdhesion