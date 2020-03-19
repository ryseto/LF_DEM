#ifndef __LF_DEM__ActivatedAdhesion_Params__
#define __LF_DEM__ActivatedAdhesion_Params__

namespace Interactions {

namespace ActAdhesion {
	struct Params {
		double activation_rate;
		double deactivation_rate;
		double max_force;
		double range;
	};
}

inline bool has_activated_adhesion(const ActAdhesion::Params &p)
{
	return p.max_force > 0;
}

}
#endif //__LF_DEM__ActivatedAdhesion_Params__
