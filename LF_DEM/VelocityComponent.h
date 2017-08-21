#ifndef __LF_DEM__VelocityComponent__
#define __LF_DEM__VelocityComponent__
#include <vector>
#include "vec3d.h"

struct VelocityComponent
{
	unsigned int rate_dependence;
	std::vector<vec3d> vel;
	std::vector<vec3d> ang_vel;
	
	VelocityComponent(){};
	VelocityComponent(std::size_t size,
					  unsigned int _rate_dependence):
	rate_dependence(_rate_dependence)
	{
		vel.resize(size);
		ang_vel.resize(size);
		reset();
	}
	
	void reset()
	{
		for (auto &v: vel) {
			v.reset();
		}
		for (auto &av: ang_vel) {
			av.reset();
		}
	}
	
	struct VelocityComponent&	operator*=(double d)
	{
		for (auto &v: vel) {
			v *= d;
		}
		for (auto &av: ang_vel) {
			av *= d;
		}
		return *this;
	}
};

#endif/* defined(__LF_DEM__VelocityComponent__) */
