#ifndef __LF_DEM__ActAdhesionState__
#define __LF_DEM__ActAdhesionState__

namespace Interactions {

namespace ActAdhesion {

// explicit numbering as it is used in output file
enum class Activity : unsigned {
	dormant = 0,
	active = 1
};

struct State {
	Activity activity;
	unsigned p0;
	unsigned p1;
};

} // namespace ActAdhesion

} // namespace Interactions

#endif