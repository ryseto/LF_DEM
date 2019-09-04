#ifndef __LF_DEM__PairVelocity__
#define __LF_DEM__PairVelocity__

#include <array>
#include "vec3d.h"
#include "Sym2Tensor.h"

struct PairVelocity
{
	// with PBC applied such that velocities are the ones seen from part 0
	std::array<vec3d, 2> U;
	std::array<vec3d, 2> O;
};

struct PairVelocityGrad
{
	std::array<vec3d, 2> Uinf;
	std::array<Sym2Tensor, 2> Einf;
	std::array<vec3d, 2> Oinf;
};

#endif // __LF_DEM__PairVelocity__