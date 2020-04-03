#include <iostream>
#include "KraynikReinelt.h"
#include "ExtensionalShearBoxSet.h"

namespace BC {

KraynikReineltBC::KraynikReineltBC(Geometry::box3d sysbox, double magic_angle, 
								   std::shared_ptr<Geometry::ImposedDeformation> imposed_deformation) :
deformation(imposed_deformation),
container(sysbox)
{
// extensional flow
	// @@@ Some variables need to be set before initializeBoxing() @@@
	// @@@ This part is very messy. Needs to be rewritten careffuly.
	strain_retrim_interval = 2*log(0.5*(3+sqrt(5))); // every this strain, the main simulation box is retrimmed.
	strain_retrim = strain_retrim_interval; // Setting the first value of strain to retrim.
	double cos_ma = cos(magic_angle);
	double sin_ma = sin(magic_angle);
	sq_cos_ma = cos_ma*cos_ma;
	sq_sin_ma = sin_ma*sin_ma;
	cos_ma_sin_ma = cos_ma*sin_ma;
	updateH(0); // cumulated_strain = 0
}

void KraynikReineltBC::setBoxSet(Boxing::ExtensionalShearBoxSet *bxset, struct Geometry::box3d *_container_ext_flow)
{
	boxset = bxset;
	container_ext_flow = _container_ext_flow;
}

void KraynikReineltBC::updateH(double cumulated_strain)
{
	/*
	 * strainH is twice of extensional strain.
	 *
	 */
	/**** extensional flow ****/
	/* H_orig = (exp(epsilon_dot t) 0 0 / 0 1 0 / 0 0 exp(-epsilon_dot t))
	 * H = R(-magicangle) H_orig R(magicangle)
	 */
	double diff_ext_strain = 0.5*(cumulated_strain-(strain_retrim-strain_retrim_interval));
	double exp_strain_x = exp(diff_ext_strain); // extensional strain = 0.5*(shear strain)
	double exp_strain_z = 1.0/exp_strain_x;
	/******** Set H     *****************************************/
	// 00, 01, 02, 12, 11, 22
	deform_forward.setSymmetric(exp_strain_x*sq_cos_ma+exp_strain_z*sq_sin_ma, // 00
								0, // 01
								-(exp_strain_x-exp_strain_z)*cos_ma_sin_ma, // 02
								0, // 12
								1, // 11
								exp_strain_x*sq_sin_ma+exp_strain_z*sq_cos_ma); // 22
	deform_backward = deform_forward.inverse();
	/************** for boxing **************************************/
	ext_ax.box_axis1 = container.lx*deform_forward.getLine(0); // 6--7 = 10--11
	ext_ax.box_diagonal_6_11 = ext_ax.box_axis1+ext_ax.box_axis2; // 6--11
	ext_ax.box_diagonal_7_10 = ext_ax.box_axis2-ext_ax.box_axis1; // 7--10
}

void KraynikReineltBC::retrim(vec3d& pos)
{
	while (pos.z >= container.lz) {
		pos.z -= container.lz;
	}
	while (pos.z < 0) {
		pos.z += container.lz;
	}
	while (pos.x >= container.lx) {
		pos.x -= container.lx;
	}
	while (pos.x < 0) {
		pos.x += container.lx;
	}
}

void KraynikReineltBC::retrimProcess(std::vector<vec3d> &position,
									 double cumulated_strain)
{
	std::cerr << "retrim" << std::endl;
	strain_retrim += strain_retrim_interval;
	updateH(cumulated_strain);
	for (unsigned i=0; i<position.size(); i++) {
		retrim(position[i]);
		boxset->box(i, position[i]);
	}
	boxset->updateExtFlow(*container_ext_flow, ext_ax);
}

void KraynikReineltBC::periodize(unsigned i, vec3d &pos, bool &pd_transport) const
{
	if (boxset->boxType(i) != 1) {
		vec3d s = deform_backward*pos;
		pd_transport = false;
		if (s.z >= container.lz) {
			s.z -= container.lz;
			pd_transport = true;
		} else if (s.z < 0) {
			s.z += container.lz;
			pd_transport = true;
		}
		if (s.x >= container.lx) {
			s.x -= container.lx;
			pd_transport = true;
		} else if (s.x < 0) {
			s.x += container.lx;
			pd_transport = true;
		}
		if (pd_transport) {
			pos = deform_forward*s;
		}
	}
	if (pos.y >= container.ly) {
		pos.y -= container.ly;
	} else if (pos.y < 0) {
		pos.y += container.ly;
	}
}

void KraynikReineltBC::periodize(unsigned i, vec3d &pos) const
{
	if (boxset->boxType(i) != 1) {
		vec3d s = deform_backward*pos;
		bool pd_transport = false;
		if (s.z >= container.lz) {
			s.z -= container.lz;
			pd_transport = true;
		} else if (s.z < 0) {
			s.z += container.lz;
			pd_transport = true;
		}
		if (s.x >= container.lx) {
			s.x -= container.lx;
			pd_transport = true;
		} else if (s.x < 0) {
			s.x += container.lx;
			pd_transport = true;
		}
		if (pd_transport) {
			pos = deform_forward*s;
		}
	}
	if (pos.y >= container.ly) {
		pos.y -= container.ly;
	} else if (pos.y < 0) {
		pos.y += container.ly;
	}
}

vec3d KraynikReineltBC::getVelDifference(const vec3d &separation) const
{
	return deformation->getGradU()*separation;
}

} // namespace BC