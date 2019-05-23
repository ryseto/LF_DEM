#ifndef __LF_DEM__Confinement_Params__
#define __LF_DEM__Confinement_Params__

#include "DimensionalQty.h"

namespace Confinement {
	struct Parameters {
		bool on;           								///< True: in use, False: no confinement [false]
		double y_min;	   								///< y_min 
		double y_max;									///< y_max
		Dimensional::DimensionalQty<double> k;          ///< potential curvature
	};
}

#endif //__LF_DEM__Confinement_Params__
