#ifndef __LF_DEM__KraynikReinelt__
#define __LF_DEM__KraynikReinelt__

#include <memory>
#include "vec3d.h"
#include "Box3d.h"
#include "Matrix.h"
#include "ImposedDeformation.h"

namespace Boxing {
	class ExtensionalShearBoxSet;
}

namespace BC {

struct ExtFlowAxis {
	vec3d box_axis1;  // x (6--7)
	vec3d box_axis2;  // z (6--10)
	vec3d box_diagonal_7_10; //
	vec3d box_diagonal_6_11;
};

class KraynikReineltBC { // ideally this should inherit from BoundaryCondition
public:
	KraynikReineltBC(Geometry::box3d sysbox, double magic_angle,
					 std::shared_ptr<Geometry::ImposedDeformation> imposed_deformation);

	/****************************************************************************************************
	 * Extensional flow using Kraynik-Reinelt Method was originally implemented                         *
	 * by Antonio Martiniello and Giulio Giuseppe Giusteri from Auguest to November 2016 at OIST.       *
	 ****************************************************************************************************/
	matrix deform_forward; // Extension flow
	matrix deform_backward; // Extension flow
	struct ExtFlowAxis ext_ax;

	struct Geometry::box3d getContainer() const {return container;};
	void periodize(unsigned i, vec3d &pos, bool &pd_transport) const;
	void periodize(unsigned i, vec3d &pos) const;

	vec3d getVelDifference(const vec3d &separation) const;
	double getStrainRetrim() const {return strain_retrim;};
	double getStrainRetrimInterval() const {return strain_retrim_interval;};

	void retrimProcess(std::vector<vec3d> &position, double cumulated_strain);
	void setBoxSet(Boxing::ExtensionalShearBoxSet *bxset, struct Geometry::box3d *_container_ext_flow);

private:
	/*****************************
	 * Domains in the simulation box are numbered as follows.
	 *    10 - 8 --11
	 *    |         |
	 *    2    1    3
	 *    |         |
	 *    6 -- 4 ---7
	 *************************/
	double strain_retrim; // APR
	double strain_retrim_interval; // APR

	double sq_cos_ma; // magic angle @@@@
	double sq_sin_ma; // magic angle @@@@
	double cos_ma_sin_ma; // magic angle @@@@
	std::shared_ptr<Geometry::ImposedDeformation> deformation;

	void updateH(double cumulated_strain);
	struct Geometry::box3d container;
	void retrim(vec3d& pos);
	Boxing::ExtensionalShearBoxSet *boxset;
	struct Geometry::box3d *container_ext_flow;
};


} // namespace BC
#endif /* defined(__LF_DEM__KraynikReinelt__) */
