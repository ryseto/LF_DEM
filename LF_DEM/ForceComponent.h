#ifndef __LF_DEM__ForceComponent__
#define __LF_DEM__ForceComponent__
#include <vector>
#include "vec3d.h"

class System;

typedef void (System::*sysGetForceTorque)(std::vector<vec3d> &force,
                                          std::vector<vec3d> &torque);

struct ForceComponent
{
  unsigned int rate_dependence;
  bool has_torque;
  std::vector<vec3d> force;
  std::vector<vec3d> torque;
  sysGetForceTorque getForceTorque;

  ForceComponent(){};
  ForceComponent(std::size_t size,
                 unsigned int _rate_dependence,
                 bool _has_torque,
                 sysGetForceTorque _getForceTorque):
                 rate_dependence(_rate_dependence),
                 has_torque(_has_torque),
                 getForceTorque(_getForceTorque)
  {
    force.resize(size);
    torque.resize(size);
    reset();
  }

  void reset()
  {
    for(auto &f: force) {
      f.reset();
    }
    for(auto &t: torque) {
      t.reset();
    }
  }

  struct ForceComponent&	operator*=(double d)
	{
    for(auto &f: force) {
      f *= d;
    }
    for(auto &t: torque) {
      t *= d;
    }
		return *this;
	}
};

#endif/* defined(__LF_DEM__ForceComponent__) */
