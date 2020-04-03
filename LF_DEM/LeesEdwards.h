#ifndef __LF_DEM__LeesEdwards__
#define __LF_DEM__LeesEdwards__

#include <memory>
#include "BoundaryCondition.h"
#include "Matrix.h"
#include "Box3d.h"
#include "ImposedDeformation.h"

namespace BC {

class LeesEdwardsBC : public BoundaryCondition {
public:
	void periodize(vec3d& pos) const;
	vec3d periodized(const vec3d &pos) const;
	void periodizeDiff(vec3d& pos_diff) const;

	vec3d getShearDisp() const {return shear_disp;};
	vec3d getShearStrain() const {return shear_strain;};
	vec3d getVelDifference(const vec3d &separation) const;
	void incrementStrain(double dt);

	
	LeesEdwardsBC(Geometry::box3d sysbox, 
				  vec3d lees_shear_disp, 
				  bool reset_strain, 
				  std::shared_ptr<Geometry::ImposedDeformation> imposed_deformation);
private:
	struct Geometry::box3d container_half;
	bool twodimension;
	vec3d shear_disp; // lees-edwards shift between top and bottom. only shear_disp.x, shear_disp.y is used
	vec3d shear_strain;

	std::shared_ptr<Geometry::ImposedDeformation> deformation;
};

matrix flowShapeSimpleShear(double theta_shear, bool radians=true);

} // namespace BC
#endif /* defined(__LF_DEM__LeesEdwards__) */
