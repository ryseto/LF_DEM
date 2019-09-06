#ifndef __LF_DEM__TimeActivatedAdhesion_Params__
#define __LF_DEM__TimeActivatedAdhesion_Params__

namespace Interactions {

namespace TActAdhesion {
	struct Parameters {
		double activation_time;
		double adhesion_max_force;
		double adhesion_range;
	};
}

inline bool has_delayed_adhesion(const TActAdhesion::Parameters &p)
{
	return p.adhesion_max_force > 0;
}

}
#endif //__LF_DEM__TimeActivatedAdhesion_Params__
