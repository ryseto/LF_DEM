#include "VelocityAssignor.h"

namespace Dynamics
{

VelocityAssignor::VelocityAssignor(ParticleVelocity fixed_velocities) :
fixed_vel(fixed_velocities)
{
	np_fixed = fixed_velocities.vel.size();
}

FixedVelocityVector::FixedVelocityVector(ParticleVelocity fixed_velocities) :
VelocityAssignor(fixed_velocities)
{
	
}


void FixedVelocityVector::set(std::vector<vec3d>::iterator vel, std::vector<vec3d>::iterator ang_vel)
{
	for (auto v: fixed_vel.vel) {
		*vel = v;
		vel++;
	}
	for (auto av: fixed_vel.ang_vel) {
		*ang_vel = av;
		ang_vel++;
	}
}


}