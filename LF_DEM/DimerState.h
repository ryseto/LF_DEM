#ifndef __LF_DEM__DimerState__
#define __LF_DEM__DimerState__

#include "vec3d.h"

namespace Interactions {

namespace Dimer {

struct SpringState {
	vec3d stretch;
	double relaxed_length;
};

struct DimerState {
	unsigned int p0;
	unsigned int p1;
	SpringState sliding_st;
	SpringState rotation_st;
};

struct UnloadedDimerState {
	unsigned int p0;
	unsigned int p1;
	double relaxed_length;
};

}

}

#endif