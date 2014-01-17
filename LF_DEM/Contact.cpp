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
	if (sys->friction_model == 1) {
		frictionlaw = &Contact::frictionlaw_coulomb;
	} else if (sys->friction_model == 2) {
		frictionlaw = &Contact::frictionlaw_criticalload;
	} else if (sys->friction_model == 3) {
		frictionlaw = &Contact::frictionlaw_criticalload_mu_inf;
	} else {
		frictionlaw = &Contact::frictionlaw_null;
	}
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
	strain_contact_start = sys->get_shear_strain();
}

void
Contact::deactivate(){
	// r > a0 + a1
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
	f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim(); // gap_nondim is negative, therefore it is allways positive.
	f_contact_normal = -f_contact_normal_norm*interaction->nvec;
	disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
	f_contact_tan = kt_scaled*disp_tan;
	(this->*frictionlaw)();
}

void
Contact::frictionlaw_coulomb(){
	interaction->lubrication.calcLubricationForce_normal(); // dashpot for the squeezed mode.
	double supportable_tanforce = f_contact_normal_norm;
	supportable_tanforce += interaction->lubrication.get_lubforce_normal_fast();
	supportable_tanforce *=	mu;
	if (supportable_tanforce < 0){
		if (staticfriction) {
			sys->incrementCounter_static_to_dynamic();
		}
		staticfriction = false;
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		double sq_f_tan = f_contact_tan.sq_norm();
		if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
			if (staticfriction) {
				sys->incrementCounter_static_to_dynamic();
			}
			staticfriction = false;
			disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		} else {
			staticfriction = true;
		}
	}
	return;
}

void
Contact::frictionlaw_criticalload(){
	interaction->lubrication.calcLubricationForce_normal();
	/* Since gap_nondim < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 * 
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm;
	supportable_tanforce += interaction->lubrication.get_lubforce_normal_fast();
	supportable_tanforce -= sys->critical_normal_force; // critical load model.
	if (supportable_tanforce < 0){
		if (staticfriction) {
			sys->incrementCounter_static_to_dynamic();
		}
 		staticfriction = false;
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		supportable_tanforce *= mu;
		double sq_f_tan = f_contact_tan.sq_norm();
		if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
			if (staticfriction) {
				sys->incrementCounter_static_to_dynamic();
			}
			staticfriction = false;
			disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		} else {
			staticfriction = true;
		}
	}
	return;
}

void
Contact::frictionlaw_criticalload_mu_inf(){
	interaction->lubrication.calcLubricationForce_normal();
	/* Since gap_nondim < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm;
	supportable_tanforce += interaction->lubrication.get_lubforce_normal_fast();
	supportable_tanforce -= sys->critical_normal_force; // critical load model.
	if (supportable_tanforce < 0){
		if (staticfriction) {
			sys->incrementCounter_static_to_dynamic();
		}
 		staticfriction = false;
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		staticfriction = true;
	}
	return;
}

void
Contact::frictionlaw_null(){;}


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
	
	/* When we compose stress tensor,
	 * even individual leverl, this part calculates symmetric tensor.
	 * This symmetry is expected in the average ensumble.
	 * I'm not sure this is allowed or not.
	 */
	if (active) {
		/*
		 * Fc_normal_norm = -kn_scaled*gap_nondim; --> positive
		 * Fc_normal = -Fc_normal_norm*nvec;
		 * This force acts on particle 1.
		 * stress1 is a0*nvec[*]force.
		 * stress2 is (-a1*nvec)[*](-force) = a1*nvec[*]force
		 * stress1 + stress2 = (a1+a2)*nvec[*]force
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

