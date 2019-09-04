#include "VelocityAssignor.h"

namespace Dynamics
{

FixedVelocityVector::FixedVelocityVector(std::vector<vec3d> fixed_velocities) :
fixed_vel(fixed_velocities)
{
	np_fixed = fixed_velocities.size();
}

void FixedVelocityVector::set(std::vector<vec3d>::iterator vel, std::vector<vec3d>::iterator ang_vel)
{
	for (auto v: fixed_vel) {
		*vel = v;
		(*ang_vel).reset();
		vel++;
		ang_vel++;
	}
}


}