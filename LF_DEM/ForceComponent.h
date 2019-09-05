#ifndef __LF_DEM__ForceComponent__
#define __LF_DEM__ForceComponent__
#include <vector>
#include "Sym2Tensor.h"
#include "RateDependence.h"

struct ForceComponent
{
	RateDependence rate_dependence;
	bool has_torque;
	std::vector<vec3d> force;
	std::vector<vec3d> torque;
	// sysGetForceTorque getForceTorque;

	ForceComponent(){};
	ForceComponent(std::size_t size,
				   RateDependence _rate_dependence,
				   bool _has_torque) :
				   // sysGetForceTorque _getForceTorque):
	rate_dependence(_rate_dependence),
	has_torque(_has_torque)
	// getForceTorque(_getForceTorque)
	{
		force.resize(size);
		torque.resize(size);
		for (auto &f: force) {
			f.reset();
		}
		for (auto &t: torque) {
			t.reset();
		}
	}

	void reset()
	{
		for (auto &f: force) {
			f.reset();
		}
		if (has_torque) {
			for (auto &t: torque) {
				t.reset();
			}
		}
	}

	template <typename T>
	struct ForceComponent& operator*=(const T& a)
	{
		for (auto &f: force) {
			f *= a;
		}
		if (has_torque) {
			for (auto &t: torque) {
				t *= a;
			}
		}
		return *this;
	}
};

#endif/* defined(__LF_DEM__ForceComponent__) */
