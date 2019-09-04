#include <iostream>
#include "BoundaryCondition.h"

namespace BC {
BoundaryCondition::BoundaryCondition(struct Geometry::box3d sysbox) :
	container(sysbox) 
{}
} // namespace BC