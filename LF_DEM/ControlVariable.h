#ifndef __LF_DEM__ControlVariable__
#define __LF_DEM__ControlVariable__

namespace Parameters {
	enum class ControlVariable : unsigned {
		rate,
		stress,
		pressure_drop,
		force
		// viscnb
	};
}

#endif