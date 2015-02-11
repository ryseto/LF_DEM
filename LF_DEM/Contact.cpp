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
	state = 0;
	if (sys->friction_model == 1) {
		frictionlaw = &Contact::frictionlaw_standard;
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
	interaction->get_par_num(p0, p1);
	double &ro_12 = interaction->ro_12;
	kn_scaled = ro_12*ro_12*sys->kn; // F = kn_scaled * _reduced_gap;  <-- gap is scaled @@@@ Why use reduced_gap? Why not gap?
	kt_scaled = ro_12*sys->kt; // F = kt_scaled * disp_tan <-- disp is not scaled
	if (sys->rolling_friction) {
		kr_scaled = ro_12; // F = kt_scaled * disp_tan <-- disp is not scaled
	}
	mu = sys->get_mu_static();
}

void
Contact::activate(){
	// r < a0 + a1
	/* state = 1 means frictionless contact.
	 * In frictional particle simulations,
	 * this value will be updated to 2 or 3 in friction law.
	 * In critical load model, the value can take 1 as well.
	 */
	state = 1;
	disp_tan.reset();
	disp_rolling.reset();
}

void
Contact::deactivate(){
	// r > a0 + a1
	state = 0;
	disp_tan.reset();
	disp_rolling.reset();
	f_contact_normal_norm = 0;
	f_contact_normal.reset();
	f_contact_tan.reset();
}

/*********************************
 *                                *
 *	   Contact Forces Methods    *
 *                                *
 *********************************/
void
Contact::incrementTangentialDisplacement(){
	if (sys->in_predictor) {
		/*
		 * relative_surface_velocity is true velocity in predictor and corrector.
		 * Thus, previous disp_tan is saved to use in corrector.
		 */
		prev_disp_tan = disp_tan;
	}
	disp_tan = prev_disp_tan+interaction->relative_surface_velocity*sys->dt; // always disp(t+1) = disp(t) + v*dt, no predictor-corrector headache :)
}

void
Contact::incrementRollingDisplacement(){
	if (sys->in_predictor) {
		/*
		 * relative_surface_velocity is true velocity in predictor and corrector.
		 * Thus, previous disp_tan is saved to use in corrector.
		 */
		prev_disp_rolling = disp_rolling;
	}
	disp_rolling = prev_disp_rolling+interaction->rolling_velocity*sys->dt; // always disp(t+1) = disp(t) + v*dt, no predictor-corrector headache :)
}

void
Contact::incrementDisplacements(){
	/**
	   \brief Increment the tangential and rolling spring stretches from relative velocities, @b without checking the friction laws.
	   
	   This should be called @b BEFORE updating the relative positions (ie normal and tangential vectors in the interaction). 
	   This is because it needs the relative velocities at time t, which depend on a variable zshift at time t which deals with Lees-Edwards PBC.
	   This zshift is updated to their value at time t+1 whenever the relative positions are computed, so updating relative positions should be done after incrementing stretches.
	 */
	if (sys->friction) {
		interaction->calcRelativeVelocities();
		incrementTangentialDisplacement();
		if (sys->rolling_friction) {
			interaction->calcRollingVelocities();
			incrementRollingDisplacement();
		}
	}
}



void
Contact::calcContactInteraction(){
	/** 
		\brief Compute the contact forces and apply friction law, by rescaling forces and stretches if necessary.
	*/
	
	/* h < 0
	 * f_contact_normal_norm > 0 ..... repulsive force
	 * h > 0
	 * f_contact_normal_norm < 0 ..... attractive force
	 *
	 *
	 reduced_gap is negative,
	 positive. */
	f_contact_normal_norm = -kn_scaled*interaction->get_reduced_gap();
	f_contact_normal = -f_contact_normal_norm*interaction->nvec;
	disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
	f_contact_tan = kt_scaled*disp_tan;
	if (sys->rolling_friction) {
		f_rolling = kr_scaled*disp_rolling;
	}
	(this->*frictionlaw)();
}

void
Contact::frictionlaw_standard(){
	double supportable_tanforce;
	if (!sys->cohesion) {
		supportable_tanforce = mu*f_contact_normal_norm;
	} else {
		supportable_tanforce = mu*(f_contact_normal_norm+sys->dimensionless_cohesive_force);
	}
	double sq_f_tan = f_contact_tan.sq_norm();
	if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
		state = 3; // sliding
		disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		f_contact_tan = kt_scaled*disp_tan;
	} else {
		state = 2; // static friction
	}
	if (sys->rolling_friction) {
		double sq_f_rolling = f_rolling.sq_norm();
		if (sq_f_rolling > supportable_tanforce*supportable_tanforce){
			disp_rolling *= supportable_tanforce/sqrt(sq_f_rolling);
			f_rolling = kr_scaled*disp_rolling;
		}
	}
	return;
}

void
Contact::frictionlaw_criticalload(){
	/* Since reduced_gap < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_reduced_gap(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 * 
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm-sys->critical_normal_force; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		supportable_tanforce *= mu;
		double sq_f_tan = f_contact_tan.sq_norm();
		if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
			state = 3; // sliding
			disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
			f_contact_tan = kt_scaled*disp_tan;
		} else {
			state = 2; // static friction
		}
	}
	return;
}

void
Contact::frictionlaw_criticalload_mu_inf(){
	/* Since reduced_gap < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_reduced_gap(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm-sys->critical_normal_force; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		// Never rescaled.
		// [note]
		// The tangential spring constant may be rescaled to control maximum strain
		state = 2; // static friction
	}
	return;
}

void
Contact::frictionlaw_null(){}

void
Contact::addUpContactForceTorque(){
	if (state > 0) {
		sys->contact_force[p0] += f_contact_normal;
		sys->contact_force[p1] -= f_contact_normal;
		if (state >= 2) {
			sys->contact_force[p0] += f_contact_tan;
			sys->contact_force[p1] -= f_contact_tan;
			vec3d t_ij = cross(interaction->nvec, f_contact_tan);
			sys->contact_torque[p0] += interaction->a0*t_ij;
			sys->contact_torque[p1] += interaction->a1*t_ij;
			if (sys->rolling_friction) {
				vec3d t_rolling = cross(interaction->nvec, f_rolling);
				sys->contact_torque[p0] += interaction->a0*t_rolling;
				sys->contact_torque[p1] -= interaction->a1*t_rolling;
			}
		}
	}
}

void
Contact::calcContactStress(){
	/*
	 * Fc_normal_norm = -kn_scaled*reduced_gap; --> positive
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
	if (state > 0) {
		/*
		 * Fc_normal_norm = -kn_scaled*reduced_gap; --> positive
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
