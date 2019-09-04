#ifndef __LF_DEM__BoundaryCondition__
#define __LF_DEM__BoundaryCondition__

#include "vec3d.h"
#include "Box3d.h"

namespace BC {

class BoundaryCondition {
public:
	virtual void periodize(vec3d& pos) const = 0;
	virtual vec3d periodized(const vec3d &pos) const = 0;
	virtual void periodizeDiff(vec3d& pos_diff) const = 0;

	BoundaryCondition(struct Geometry::box3d sysbox);
	struct Geometry::box3d getContainer() const {return container;};

protected:
	struct Geometry::box3d container;
};

} // namespace BC
#endif /* defined(__LF_DEM__Geometry__) */
