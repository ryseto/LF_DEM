#ifndef __LF_DEM__vanDerWaalsParams__
#define __LF_DEM__vanDerWaalsParams__

namespace Interactions
{

struct vanDerWaalsForceParams
{
	double singularity_cutoff;
	double coefficient;
	double amplitude;
};

}

#endif