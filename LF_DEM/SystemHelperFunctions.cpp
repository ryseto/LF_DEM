//
//  SystemHelperFunctions.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <stdexcept>

using namespace std;


void System::evaluateMaxContactVelocity()
{
	max_contact_velo_tan = 0;
	max_contact_velo_normal = 0;
	max_relative_velocity = 0;
	double sum_contact_velo_tan = 0;
	double sum_contact_velo_normal = 0;
	double sum_sliding_velocity = 0;
	int cnt_contact = 0;
	int cnt_sliding = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			cnt_contact++;
			sum_contact_velo_tan += interaction[k].getContactVelocity();
			sum_contact_velo_normal += abs(interaction[k].getNormalVelocity());
			if (interaction[k].getContactVelocity() > max_contact_velo_tan) {
				// relative_surface_velocity for both static friction and sliding state.
				max_contact_velo_tan = interaction[k].getContactVelocity();
			}
			if (abs(interaction[k].getNormalVelocity()) > max_contact_velo_normal) {
				max_contact_velo_normal = abs(interaction[k].getNormalVelocity());
			}
			if (interaction[k].getRelativeVelocity() > max_relative_velocity) {
				max_relative_velocity = interaction[k].getRelativeVelocity();
			}
			if (interaction[k].contact.state == 3) {
				/*
				 * relative_surface_velocity for only sliding state.
				 */
				cnt_sliding++;
				sum_sliding_velocity += interaction[k].getContactVelocity();
			}
		}
	}
	if (cnt_contact > 0) {
		ave_contact_velo_tan = sum_contact_velo_tan/cnt_contact;
		ave_contact_velo_normal = sum_contact_velo_normal/cnt_contact;
	} else {
		ave_contact_velo_tan = 0;
		ave_contact_velo_normal = 0;
	}
	if (cnt_sliding > 0) {
		ave_sliding_velocity = sum_sliding_velocity/cnt_sliding;
	} else {
		ave_sliding_velocity = 0;
	}
}

double System::evaluateMaxVelocity()
{
	double sq_max_velocity = 0;
	for (int i = 0; i < np; i++) {
		vec3d na_velocity_tmp = velocity[i];
		if (!zero_shear) {
			if (!p.cross_shear) {
				na_velocity_tmp.x -= shear_rate*position[i].z;
			} else {
				na_velocity_tmp.x -= costheta_shear*shear_rate*position[i].z;
				na_velocity_tmp.y -= sintheta_shear*shear_rate*position[i].z;
			}
		}
		if (na_velocity_tmp.sq_norm() > sq_max_velocity) {
			sq_max_velocity = na_velocity_tmp.sq_norm();
		}
	}
	return sqrt(sq_max_velocity);
}

double System::evaluateMaxAngVelocity()
{
	double _max_ang_velocity = 0;
	for (int i = 0; i < np; i++) {
		vec3d na_ang_velocity_tmp = ang_velocity[i];
		if (!zero_shear) {
			if (!p.cross_shear) {
				na_ang_velocity_tmp.y -= 0.5*shear_rate;
			}
			else {
				na_ang_velocity_tmp.y -= 0.5*costheta_shear*shear_rate;
				na_ang_velocity_tmp.x += 0.5*sintheta_shear*shear_rate;
			}
		}
		if (na_ang_velocity_tmp.norm() > _max_ang_velocity) {
			_max_ang_velocity = na_ang_velocity_tmp.norm();
		}
	}
	return _max_ang_velocity;
}

double System::evaluateMinGap()
{
	double _min_reduced_gap = p.lub_max_gap;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].get_reduced_gap() < _min_reduced_gap) {
			_min_reduced_gap = interaction[k].get_reduced_gap();
		}
	}
	return _min_reduced_gap;
}

double System::evaluateAvgContactGap()
{
	double _avg_reduced_gap = 0;
	int nb=0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.state > 0) {
			_avg_reduced_gap += interaction[k].get_reduced_gap();
			nb++;
		}
	}
	_avg_reduced_gap /= nb;
	return _avg_reduced_gap;
}

double System::evaluateMaxContactGap()
{
	double _max_contact_gap = 0;
	for (int k= 0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.state > 0 &&
			interaction[k].get_reduced_gap() > _max_contact_gap) {
			_max_contact_gap = interaction[k].get_reduced_gap();
		}
	}
	return _max_contact_gap;
}

double System::evaluateMaxDispTan()
{
	double _max_disp_tan = 0;
	for (int k= 0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.disp_tan.norm() > _max_disp_tan) {
			_max_disp_tan = interaction[k].contact.disp_tan.norm();
		}
	}
	return _max_disp_tan;
}

double System::evaluateMaxDispRolling()
{
	double _max_disp_rolling = 0;
	for (int k= 0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.disp_rolling.norm() > _max_disp_rolling) {
			_max_disp_rolling = interaction[k].contact.disp_rolling.norm();
		}
	}
	return _max_disp_rolling;
}

double System::evaluateMaxFcNormal()
{
	double max_fc_normal_ = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].is_contact() &&
			interaction[k].contact.get_f_contact_normal_norm() > max_fc_normal_) {
			max_fc_normal_ = interaction[k].contact.get_f_contact_normal_norm();
		}
	}
	return max_fc_normal_;
}

double System::evaluateMaxFcTangential()
{
	double max_fc_tan_ = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].is_contact() &&
			interaction[k].contact.get_f_contact_tan_norm() > max_fc_tan_) {
			max_fc_tan_ = interaction[k].contact.get_f_contact_tan_norm();
		}
	}
	return max_fc_tan_;
}

void System::countNumberOfContact()
{
	contact_nb = 0;
	fric_contact_nb = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].is_contact()) {
			contact_nb ++;
			if (interaction[k].is_friccontact()) {
				fric_contact_nb ++;
			}
		}
	}
}

void System::analyzeState()
{
	//	max_velocity = evaluateMaxVelocity();
	computeMaxNAVelocity();
	max_ang_velocity = evaluateMaxAngVelocity();
	if (friction) {
		evaluateMaxContactVelocity();
	}
	min_reduced_gap = evaluateMinGap();
	if (cohesion) {
		max_contact_gap = evaluateMaxContactGap();
	}
	max_disp_tan = evaluateMaxDispTan();
	if (rolling_friction) {
		max_disp_rolling = evaluateMaxDispRolling();
	}
	max_fc_normal = evaluateMaxFcNormal();
	max_fc_tan = evaluateMaxFcTangential();
	countNumberOfContact();
	calcPotentialEnergy();
}

void System::calcPotentialEnergy()
{
	total_energy = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (interaction[k].is_contact()){
				total_energy += interaction[k].contact.calcEnergy();
			}
			if (repulsiveforce) {
				total_energy += interaction[k].repulsion.calcEnergy();
			}
		}
	}
	if (magnetic) {
		calcMagneticEnergy();
	}
}

void System::calcMagneticEnergy()
{
	/*
	 Magnetic energy is given in the thrmal unit.

	 E_M = -(mu_f)/(4 pi r**3) (3 (m1.n)(m2.n)-m1.m2)
	 \tilde{E}_M = E_M/kT = (2/3)*Pe_M*\hat{E}_M
	 \hat{E}_M = E_M / E_M^{*} = 8*(3*m0.n*m1.n-m0.m1)/r**3;
	 */

	magnetic_dd_energy = 0;
	vec3d pos_diff;
	vec3d nvec;
	for (int p0=0; p0<np-1; p0++) {
		for (const int p1: magnetic_pair[p0]) {
			pos_diff = position[p1]-position[p0];
			periodize_diff(pos_diff);
			double r = pos_diff.norm();
			nvec = pos_diff/r;
			/*
			 * amplitudes.magnetic = Pe_magnetic
			 */
			double m0_n = dot(magnetic_moment[p0], nvec);
			double m1_n = dot(magnetic_moment[p1], nvec);
			double m0_m1 = dot(magnetic_moment[p0], magnetic_moment[p1]);
			double energy_mag = -8*(3*m0_n*m1_n-m0_m1)/(r*r*r); // = E_M/E_M0
			energy_mag *= (2.0/3)*amplitudes.magnetic; // = E_M/E_B0
			magnetic_dd_energy += energy_mag; //
			total_energy += energy_mag;
		}
	}
	if (p.magnetic_type == 1) {
		if (external_magnetic_field.is_not_zero()) {
			throw runtime_error( "not implemented yet @ calcMagneticEnergy");
			for (int i=0; i<np; i++) {
				double tmp_magnetic_energy_ex = -dot(magnetic_moment[i], external_magnetic_field);
				total_energy += tmp_magnetic_energy_ex;
				magnetic_dd_energy += tmp_magnetic_energy_ex;
			}
		}
	}
	magnetic_dd_energy = magnetic_dd_energy/np;
}
