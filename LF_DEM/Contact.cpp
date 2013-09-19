//
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Contact.h"
#include "Interaction.h"

void
Contact::init(System *sys_, Interaction *interaction_){
	sys = sys_;
	interaction = interaction_;
	active = false;
}

void
Contact::getInteractionData(){
	interaction->get_par_num(i, j);
	double ro_12 = interaction->ro_12;
	kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
	kt_scaled = ro_12*sys->get_kt(); // F = kt_scaled * disp_tan <-- disp is not scaled
	mu = sys->get_mu_static();
}

/*
 * This is just used in the parameter determination code.
 */
void
Contact::updateContactModel(){
	exit(1);
	if (interaction->active) {
		double ro_12 = interaction->ro_12;
		kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
		kt_scaled = ro_12*sys->get_kt(); // F = kt_s
		if (active) {
			double lub_coeff = sys->get_lub_coeff_contact();
			double log_lub_coeff = sys->get_log_lub_coeff_dynamicfriction();
			interaction->lubrication.setResistanceCoeff(lub_coeff, log_lub_coeff);
		}
	}
}

void
Contact::resetObservables(){
	// observables
	strain_contact_start = sys->get_shear_strain();
	duration_contact = 0; // for output
	cnt_sliding = 0;
}

void
Contact::activate(){
	// r < a0 + a1
	active = true;
	staticfriction = true;
	disp_tan.reset();
	double lub_coeff = sys->get_lub_coeff_contact();
	double log_lub_coeff = sys->get_log_lub_coeff_dynamicfriction();
	interaction->lubrication.setResistanceCoeff(lub_coeff, log_lub_coeff);
}

void
Contact::deactivate(){
	// r > a0 + a1
#ifdef RECORD_HISTORY
	outputHistory();
#endif
	active = false;
	disp_tan.reset();
	f_contact_normal_norm = 0;
	f_contact_normal.reset();
	f_contact_tan.reset();
	duration_contact += sys->get_shear_strain()-strain_contact_start; // for output
}

/*********************************
 *                                *
 *	   Contact Forces Methods    *
 *                                *
 *********************************/

void
Contact::incrementTangentialDisplacement(){
	disp_tan += interaction->relative_surface_velocity*sys->get_dt();
}

/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * rvec = p[1] - p[0]
 * Fc_normal_norm is positive (for overlapping particles r < ro)
 */
void
Contact::calcContactInteraction(){
	frictionlaw();
	f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim(); // gap_nondim is negative, therefore it is allways positive.
	f_contact_normal = -f_contact_normal_norm*interaction->nvec;
	disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
	f_contact_tan = kt_scaled*disp_tan;
}

void
Contact::frictionlaw(){
	interaction->lubrication.calcLubricationForce();
	supportable_tanforce = mu*(f_contact_normal_norm+interaction->lubrication.get_lubforce_normal());
	if (supportable_tanforce < 0){
		supportable_tanforce = 0;
	}
	double sq_f_tan = f_contact_tan.sq_norm();
	if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
		staticfriction = false;
		disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		cnt_sliding++; // for output
	} else {
		staticfriction = true;
	}
}



void
Contact::calcContactInteractionRelax(){
	f_contact_normal_norm = -kn_scaled*(interaction->get_gap_nondim()-0.02);
	f_contact_normal = -f_contact_normal_norm*interaction->nvec;
}

void
Contact::addUpContactForceTorque(){
	if (active) {
		sys->contact_force[i] += f_contact_normal;
		sys->contact_force[j] -= f_contact_normal;
		if (sys->friction) {
			sys->contact_force[i] += f_contact_tan;
			sys->contact_force[j] -= f_contact_tan;
			vec3d t_ij = cross(interaction->nvec, f_contact_tan);
			sys->contact_torque[i] += interaction->a0*t_ij;
			sys->contact_torque[j] += interaction->a1*t_ij;
		}
	}
}

void
Contact::addContactStress(){
	/*
	 * Fc_normal_norm = -kn_scaled*gap_nondim; --> positive
	 * Fc_normal = -Fc_normal_norm*nvec;
	 * This force acts on particle 1.
	 * stress1 is a0*nvec[*]force.
	 * stress2 is (-a1*nvec)[*](-force) = a1*nvec[*]force
	 */
	if (active) {
		/*
		 * Fc_normal_norm = -kn_scaled*gap_nondim; --> positive
		 * Fc_normal = -Fc_normal_norm*nvec;
		 * This force acts on particle 1.
		 * stress1 is a0*nvec[*]force.
		 * stress2 is (-a1*nvec)[*](-force) = a1*nvec[*]force
		 */
		contact_stresslet_XF_normal.set(interaction->rvec, f_contact_normal);
		contact_stresslet_XF_tan.set(interaction->rvec, f_contact_tan);
	} else {
		contact_stresslet_XF_normal.reset();
		contact_stresslet_XF_tan.reset();
	}

}

double
Contact::getContactVelocity(){
	if (!active) {
		return 0;
	}
	return interaction->relative_surface_velocity.norm();
}

