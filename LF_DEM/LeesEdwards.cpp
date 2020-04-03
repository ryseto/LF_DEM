#include <iostream>
#include "LeesEdwards.h"

namespace BC {

void LeesEdwardsBC::periodize(vec3d &pos) const
{
	/* Lees-Edwards boundary condition
	 *
	 */
	while (pos.z >= container.lz) {
		pos.z -= container.lz;
		pos -= shear_disp;
	} 
	while (pos.z < 0) {
		pos.z += container.lz;
		pos += shear_disp;
	}
	while (pos.x >= container.lx) {
		pos.x -= container.lx;
	}
	while (pos.x < 0) {
		pos.x += container.lx;
	}
	if (!twodimension) {
		while (pos.y >= container.ly) {
			pos.y -= container.ly;
		} 
		while (pos.y < 0) {
			pos.y += container.ly;
		}
	}
}

vec3d LeesEdwardsBC::periodized(const vec3d &pos) const
{
	vec3d periodized_pos = pos;
	periodize(periodized_pos);
	return periodized_pos;
}

void LeesEdwardsBC::periodizeDiff(vec3d &pos_diff) const
{
	/** Periodize a separation vector with Lees-Edwards boundary condition
	 */
	while (pos_diff.z > container_half.lz) {
		pos_diff.z -= container.lz;
		pos_diff -= shear_disp;
	} 
	while (pos_diff.z < -container_half.lz) {
		pos_diff.z += container.lz;
		pos_diff += shear_disp;
	}
	while (pos_diff.x > container_half.lx) {
		pos_diff.x -= container.lx;
	}
	while (pos_diff.x < -container_half.lx) {
		pos_diff.x += container.lx;
	}
	if (!twodimension) {
		while (pos_diff.y > container_half.ly) {
			pos_diff.y -= container.ly;
		} 
		while (pos_diff.y < -container_half.ly) {
			pos_diff.y += container.ly;
		}
	}
}

LeesEdwardsBC::LeesEdwardsBC(Geometry::box3d sysbox, 
							 vec3d lees_shear_disp, 
							 bool keep_strain, 
							 std::shared_ptr<Geometry::ImposedDeformation> imposed_deformation) :
BoundaryCondition(sysbox),
twodimension(sysbox.ly == 0),
shear_disp(lees_shear_disp),
deformation(imposed_deformation)
{
	container_half = {0.5*container.lx, 0.5*container.ly, 0.5*container.lz};
	if (keep_strain) {
		shear_strain = shear_disp/container.lz;
	} else {
		shear_strain = {0, 0, 0};
	}
}

void LeesEdwardsBC::incrementStrain(double dt)
{
	vec3d shear_strain_increment = 2*dot(deformation->getSymGradU(), {0, 0, 1})*dt;
	shear_strain += shear_strain_increment;
	shear_disp += shear_strain_increment*container.lz;
	int m = (int)(shear_disp.x/container.lx);
	if (shear_disp.x < 0) {
		m--;
	}
	shear_disp.x = shear_disp.x-m*container.lx;
	if (!twodimension) {
		m = (int)(shear_disp.y/container.ly);
		if (shear_disp.y < 0) {
			m--;
		}
		shear_disp.y = shear_disp.y-m*container.ly;
	}
}

vec3d LeesEdwardsBC::getVelDifference(const vec3d &separation) const
{
	return deformation->getGradU()*separation;
}

matrix flowShapeSimpleShear(double theta_shear, bool radians) {
	if (!radians) {
		theta_shear *= M_PI/180.;
	}
	double costheta_shear = cos(theta_shear);
	double sintheta_shear = sin(theta_shear);
	if (std::abs(sintheta_shear) < 1e-15) {
		sintheta_shear = 0;
		if (costheta_shear > 0) {
			costheta_shear = 1;
		} else {
			costheta_shear = -1;
		}
	} else if (std::abs(costheta_shear) < 1e-15) {
		costheta_shear = 0;
		if (sintheta_shear > 0) {
			sintheta_shear = 1;
		} else {
			sintheta_shear = -1;
		}
	}
	matrix flow_shape (0, 0, costheta_shear, 
					   0, 0, sintheta_shear,
					   0, 0, 0              );
	return flow_shape;
}

} // namespace BC