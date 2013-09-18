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
	if (interaction->active) {
		double ro_12 = interaction->ro_12;
		kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
		kt_scaled = ro_12*sys->get_kt(); // F = kt_s
		if (active) {
			double lub_coeff = sys->get_lub_coeff_contact();
			double log_lub_coeff;
			if (staticfriction) {
				log_lub_coeff = sys->get_log_lub_coeff_staticfriction();
			} else {
				log_lub_coeff = sys->get_log_lub_coeff_dynamicfriction();
			}
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
	staticfriction = false;
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
	staticfriction = false;
	disp_tan.reset();
	f_contact_normal_norm = 0;
	f_contact_normal.reset();
	f_contact_tan.reset();
	duration_contact += sys->get_shear_strain()-strain_contact_start; // for output


	previous_f_test = 0;
	previous_supportable_tanforce = 0;
	old_state = 0;
	old_relative_velocity.reset();
	old_f_test_vec.reset();
	old_lubforce_tan.reset();
	old_dashpot.reset();
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

	interaction->lubrication.calcLubricationForce();

	f_contact_normal_norm = -kn_scaled*interaction->get_gap_nondim(); // gap_nondim is negative, therefore it is allways positive.
	f_contact_normal = -f_contact_normal_norm*interaction->nvec;
	if (sys->friction) {
		/* disp_tan is orthogonal to the normal vector.
		 */
		if (staticfriction) {
			/* Static friction:
			 *
			 */
			disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
			f_contact_tan = kt_scaled*disp_tan;
		} else {
			/* Dynamic friction:
			 * The spring force is set with the new value of mu*F_n.
			 * The total resistance force consists of lubrication force and dashpot.
			 * In the contact model, the dashpot part is considered to set the tangential spring.
			 * As approximation, the previously calculated resistance force is used here.
			 */

			disp_tan = (1/kt_scaled)*(supportable_tanforce*tvec-sys->get_ratio_dashpot_total()*resforce_tan);
			disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
			f_contact_tan = kt_scaled*disp_tan;
		}
	}
}

void
Contact::frictionlaw(){
	interaction->lubrication.calcLubricationForce();
	supportable_tanforce = mu*(f_contact_normal_norm+interaction->lubrication.get_lubforce_normal());
	double f_test;

	if (staticfriction){
		/* Friction force is the sum of spring force and dashpot force.
		 * The lubrication forces obtained by solving the force balance equations
		 * are considered as dashpot force.
		 */

		/* resforce_tan is dashpot
		 */
		resforce_tan = interaction->lubrication.get_lubforce_tan();
		f_test_vec = f_contact_tan+resforce_tan;
		f_test = f_test_vec.norm();
		tvec = f_test_vec/f_test;
		if (f_test > supportable_tanforce){

			/* If f_est becomes larger than the supportable tangential force,
			 * the state is switched into the dynamics friction.
			 */
			staticfriction = false;
			sys->fout_dfric << previous_f_test << ' ' << previous_supportable_tanforce << endl;
			sys->fout_dfric << f_test << ' ' << supportable_tanforce << endl << endl;
			double log_lub_coeff = sys->get_log_lub_coeff_dynamicfriction();
			interaction->lubrication.setResistanceCoeffTang(log_lub_coeff);
		}
	} else {

		/* resforce_tan is lubforce + dashpot
		 */ 
		resforce_tan = interaction->lubrication.get_lubforce_tan();
		f_test_vec = f_contact_tan+resforce_tan;
		f_test = f_test_vec.norm();
		tvec = f_test_vec/f_test;
		if (f_test < supportable_tanforce) {
			staticfriction = true;
			disp_tan = (1/kt_scaled)*(supportable_tanforce*tvec-resforce_tan);
			double log_lub_coeff = sys->get_log_lub_coeff_staticfriction();
			interaction->lubrication.setResistanceCoeffTang(log_lub_coeff);
			if (previous_f_test != 0){
				sys->fout_sfric << previous_f_test << ' ' << previous_supportable_tanforce << endl;
				sys->fout_sfric << f_test << ' ' << supportable_tanforce << endl << endl;
			}
		}
	}
	previous_f_test = f_test;
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

