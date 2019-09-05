#ifndef __LF_DEM__ParticleConfig__
#define __LF_DEM__ParticleConfig__

#include <vector>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "RateDependence.h"

class ParticleConfig {
public:
	std::vector <vec3d> position;
	std::vector <double> radius;
	std::vector <double> angle;  // 2d only
	ParticleConfig(){};
	ParticleConfig(unsigned np) :
	position(np),
	radius(np),
	angle(np) {};
};

enum class VelocityType : unsigned {
	total,
	nonaffine
};

class ParticleVelocity
{
public:
	VelocityType type;
	RateDependence rate_dependence;
	std::vector<vec3d> vel;
	std::vector<vec3d> ang_vel;
	ParticleVelocity(){};	
	ParticleVelocity(std::size_t size,
					 VelocityType _type,
					 RateDependence _rate_dependence):
	type(_type),
	rate_dependence(_rate_dependence),
	vel(size, 0),
	ang_vel(size, 0)
	{}
	
	void reset()
	{
		for (auto &v: vel) {
			v.reset();
		}
		for (auto &av: ang_vel) {
			av.reset();
		}
	}
	
	class ParticleVelocity&	operator*=(double d)
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

class ParticleVelocityGrad
{
public:
	std::vector<Sym2Tensor> E;
	ParticleVelocityGrad() {};
	ParticleVelocityGrad(std::size_t size) :
	E(size, 0) {};
};

// class ParticleKinematics {
// public:
// 	std::vector<vec3d> U;  // translational vel
// 	std::vector<vec3d> O;  // rotational vel
// 	std::vector<vec3d> U_na;  // non-affine translational vel
// 	std::vector<vec3d> O_na;  // non-affine rotational vel
// 	std::vector<vec3d> Uinf;  // fluid velocity at particle position "as if particle was not there"
// 	std::vector<vec3d> Oinf;      // fluid vorticity
// 	std::vector<Sym2Tensor> Einf;     // fluid symmetrized velocity grad
// 	ParticleKinematics(unsigned np) :
// 	U(np),
// 	O(np),
// 	U_na(np),
// 	O_na(np),
// 	Uinf(np),
// 	Oinf(np),
// 	Einf(np) {};
// };

#endif