#ifndef __LF_DEM__VelocityAssignor__
#define __LF_DEM__VelocityAssignor__

#include <vector>
#include "vec3d.h"
#include "ParticleConfig.h"

namespace Dynamics {

class VelocityAssignor {
public:
	VelocityAssignor(ParticleVelocity fixed_velocities);
	virtual	void set(std::vector<vec3d>::iterator vel, std::vector<vec3d>::iterator ang_vel) = 0;
	ParticleVelocity getFixedVel() {return fixed_vel;};
	
	unsigned np_fixed;
protected:
	ParticleVelocity fixed_vel;
};


class FixedVelocityVector : public VelocityAssignor {
public:
	FixedVelocityVector(ParticleVelocity fixed_velocities);
	void set(std::vector<vec3d>::iterator vel, std::vector<vec3d>::iterator ang_vel);

private:
};

}

#endif