#ifndef __LF_DEM__VelocityAssignor__
#define __LF_DEM__VelocityAssignor__

#include <vector>
#include "vec3d.h"

namespace Dynamics {

class VelocityAssignor {
public:
	virtual	void set(std::vector<vec3d>::iterator vel, std::vector<vec3d>::iterator ang_vel) = 0;
	unsigned np_fixed;
};


class FixedVelocityVector : public VelocityAssignor {
public:
	FixedVelocityVector(std::vector<vec3d> fixed_velocities);
	void set(std::vector<vec3d>::iterator vel, std::vector<vec3d>::iterator ang_vel);
private:
	std::vector<vec3d> fixed_vel;
};

}

#endif