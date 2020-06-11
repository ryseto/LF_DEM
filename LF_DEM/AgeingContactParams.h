#ifndef __LF_DEM__AgeingContactParams__
#define __LF_DEM__AgeingContactParams__

namespace Interactions {

struct AgeingContactParams {
	double amplitude;
	double timescale;
};

inline bool has_ageing_contacts(const AgeingContactParams &p)
{
	return p.timescale != 0;
}

} // namespace Interaction

#endif