#ifndef __LF_DEM__RepulsiveForceParams__
#define __LF_DEM__RepulsiveForceParams__

namespace Interactions {

struct RepulsiveForceParams
{
	double screening_length;    ///<  Repulsive force proportional to exp(-reduced_gap/screening_length) [0.05]
	double max_length;			///<  Assume force is 0 for gap>max_length, if max_length == -1, max_length set to 7*screening_length [-1]
	double repulsion;           ///<  Amplitude of the repulsive force at contact  [0 guarranted unit]
	double smoothing;			///<  Smooth the exponential down to 0 at max_length with a tanh of width smoothing [0]
};

bool has_repulsion(struct RepulsiveForceParams p)
{
	return p.repulsion > 0;
}

}

#endif /* defined(__LF_DEM__RepulsiveForceParams__) */
