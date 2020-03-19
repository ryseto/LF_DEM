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


double evaluateMaxNAVelocityComponent(const System & sys, std::string component)
{
	double sq_max_na_velocity = 0;
	for (unsigned i=0; i<sys.get_np(); i++) {
		auto sq_na_velocity = sys.na_velo_components.at(component).vel[i].sq_norm();
		if (sq_na_velocity > sq_max_na_velocity) {
				sq_max_na_velocity = sq_na_velocity;
		}
	}
	return std::sqrt(sq_max_na_velocity);
}

double evaluateMaxAngVelocity(const System & sys)
{
	double _max_ang_velocity = 0;
	for (unsigned i = 0; i < sys.get_np(); i++) {
		vec3d na_ang_velocity_tmp = sys.na_velocity.ang_vel[i];
		if (na_ang_velocity_tmp.norm() > _max_ang_velocity) {
			_max_ang_velocity = na_ang_velocity_tmp.norm();
		}
	}
	return _max_ang_velocity;
}


double evaluateMinGap(const System & sys)
{
	double _min_reduced_gap = 1000;
	unsigned int p0, p1;
	for (const auto &inter: *(sys.interaction)) {
		std::tie(p0, p1) = inter->get_par_num();
		if ((int)p0 < sys.np_mobile && // exclude fixed-fixed
			inter->getReducedGap() < _min_reduced_gap) {
			_min_reduced_gap = inter->getReducedGap();
		}
	}
	return _min_reduced_gap;
}

double evaluateAvgContactGap(const System & sys)
{
	double _avg_reduced_gap = 0;
	int nb=0;
	for (const auto &inter: *(sys.interaction)) {
		if (inter->contact) {
			_avg_reduced_gap += inter->getReducedGap();
			nb++;
		}
	}
	_avg_reduced_gap /= nb;
	return _avg_reduced_gap;
}

double evaluateMaxContactGap(const System & sys)
{
	double _max_contact_gap = 0;
	for (const auto &inter: *(sys.interaction)) {
		if (inter->contact &&
			inter->getReducedGap() > _max_contact_gap) {
			_max_contact_gap = inter->getReducedGap();
		}
	}
	return _max_contact_gap;
}

double evaluateMaxVelocity(const System &sys)
{
	double sq_max_velocity = 0;
	for (unsigned i = 0; i < sys.get_np(); i++) {
		vec3d na_velocity_tmp = sys.na_velocity.vel[i];
		if (na_velocity_tmp.sq_norm() > sq_max_velocity) {
			sq_max_velocity = na_velocity_tmp.sq_norm();
		}
	}
	return sqrt(sq_max_velocity);
}

double evaluateMaxDispTan(const System &sys)
{
	double _max_disp_tan = 0;
	for (const auto &inter: *(sys.interaction)) {
		if (inter->contact) {
			if (inter->contact->disp_tan.norm() > _max_disp_tan) {
				_max_disp_tan = inter->contact->disp_tan.norm();
			}
		}
	}
	return _max_disp_tan;
}

double evaluateMaxDispRolling(const System &sys)
{
	double _max_disp_rolling = 0;
	for (const auto &inter: *(sys.interaction)) {
		if (inter->contact) {
			if (inter->contact->disp_rolling.norm() > _max_disp_rolling) {
				_max_disp_rolling = inter->contact->disp_rolling.norm();
			}
		}
	}
	return _max_disp_rolling;
}


std::pair<unsigned int, unsigned int> countNumberOfContact(const System &sys)
{
	unsigned int contact_nb = 0;
	unsigned int fric_contact_nb = 0;
	for (const auto &inter: *(sys.interaction)) {
		if (inter->contact) {
			contact_nb ++;
			if (inter->contact->is_frictional()) {
				fric_contact_nb ++;
			}
		}
	}
	return std::make_pair(contact_nb, fric_contact_nb);
}

std::pair<unsigned, double> getActAdhesionActivityStatistics(const System &sys)
{
	assert(Interactions::has_activated_adhesion(sys.p.activated_adhesion));
	unsigned active_nb = 0;
	unsigned total_nb = 0;
	for (const auto &inter: *(sys.interaction)) {
		if (inter->act_adhesion) {
			total_nb++;
			if (inter->act_adhesion->getState().activity == Interactions::ActAdhesion::Activity::active) {
				active_nb++;
			}
		}
	}
	double active_ratio = 0;
	if (total_nb>0) {
		active_ratio = (double)(active_nb)/total_nb;
	}
	return std::make_pair(active_nb, active_ratio);
}

double getPotentialEnergy(const System &sys)
{
	double total_energy = 0;
	for (const auto &inter: *(sys.interaction)) {
		if (inter->contact){
			total_energy += inter->contact->calcEnergy();
		}
		if (inter->repulsion) {
			total_energy += inter->repulsion->calcEnergy();
		}
	}
	return total_energy;
}

void isInContact(const System &sys, std::vector<int> &isincontact)
{
	isincontact.resize(sys.get_np(), 0);
	for (const auto &inter: *(sys.interaction)) {
		if (inter->contact) {
			unsigned p0, p1;
			std::tie(p0, p1) = inter->get_par_num();	
			isincontact[p0] = 1;
			isincontact[p1] = 1;
		}
	}
}