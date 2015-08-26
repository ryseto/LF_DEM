//
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Contact.h"
#include "Interaction.h"

void Contact::init(System *sys_, Interaction *interaction_)
{
	sys = sys_;
	interaction = interaction_;
	state = 0;
	f_contact_normal_norm = 0;
	if (sys->p.friction_model != 0) {
		if (sys->p.friction_model == 1) {
			frictionlaw = &Contact::frictionlaw_standard;
		} else if (sys->p.friction_model == 2) {
			frictionlaw = &Contact::frictionlaw_criticalload;
		} else if (sys->p.friction_model == 3) {
			frictionlaw = &Contact::frictionlaw_criticalload_mu_inf;
		} else if (sys->p.friction_model == 5) {
	 	frictionlaw = &Contact::frictionlaw_ft_max;
			ft_max = sys->p.ft_max;
		} else if (sys->p.friction_model == 6) {
			frictionlaw = &Contact::frictionlaw_coulomb_max;
			ft_max = sys->p.ft_max;
		}
	}
}

void Contact::setSpringConstants()
{
	const double &ro_12 = interaction->ro_12;
	kn_scaled = ro_12*ro_12*sys->p.kn; // F = kn_scaled * _reduced_gap;  <-- gap is scaled @@@@ Why use reduced_gap? Why not gap?
	if (sys->friction) {
		kt_scaled = ro_12*sys->p.kt; // F = kt_scaled * disp_tan <-- disp is not scaled
	}
	if (sys->rolling_friction) {
		kr_scaled = ro_12*sys->p.kr; // F = kt_scaled * disp_tan <-- disp is not scaled
	}
}

void Contact::setInteractionData()
{
	interaction->get_par_num(p0, p1);
	setSpringConstants();
	if (sys->friction) {
		mu_static = sys->p.mu_static;
		mu_dynamic = sys->p.mu_dynamic;
		if (sys->rolling_friction) {
			mu_rolling = sys->p.mu_rolling;
		}
	}
}

void Contact::activate()
{
	// r < a0 + a1
	/* state = 1 means frictionless contact.
	 * In frictional particle simulations,
	 * this value will be updated to 2 or 3 in friction law.
	 * In critical load model, the value can take 1 as well.
	 */
	if (sys->friction) {
		if (sys->p.friction_model == 2 || sys->p.friction_model == 3) {
			state = 1; // critical load model
		} else {
			state = 2; // static friction
		}
		disp_tan.reset();
		disp_rolling.reset();
	} else {
		state = 1;
	}
}

void Contact::deactivate()
{
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
void Contact::incrementTangentialDisplacement()
{
	if (sys->in_predictor) {
		/*
		 * relative_surface_velocity is true velocity in predictor and corrector.
		 * Thus, previous disp_tan is saved to use in corrector.
		 */
		prev_disp_tan = disp_tan;
	}
	disp_tan = prev_disp_tan+interaction->relative_surface_velocity*sys->dt; // always disp(t+1) = disp(t) + v*dt, no predictor-corrector headache :)
}

void Contact::incrementRollingDisplacement()
{
	if (sys->in_predictor) {
		/*
		 * relative_surface_velocity is true velocity in predictor and corrector.
		 * Thus, previous disp_tan is saved to use in corrector.
		 */
		prev_disp_rolling = disp_rolling;
	}
	disp_rolling = prev_disp_rolling+interaction->rolling_velocity*sys->dt; // always disp(t+1) = disp(t) + v*dt, no predictor-corrector headache :)
}

void Contact::incrementDisplacements()
{
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

void Contact::calcContactInteraction()
{
	/** 
		\brief Compute the contact forces and apply friction law, by rescaling forces and stretches if necessary.
	
	 f_something is the force acting on the particle 0.
	 disp_something is the relative displacement of the particle 1 from the particle 0.
	 Therefore, the sign of force is same as the one of the displacement.
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
	if (sys->friction) {
		disp_tan.vertical_projection(interaction->nvec);
		f_contact_tan = kt_scaled*disp_tan;
		if (sys->rolling_friction) {
			f_rolling = kr_scaled*disp_rolling;
		}
		(this->*frictionlaw)();
	}
	//	calcScaledForce();
}

void Contact::frictionlaw_standard()
{
	/**
	 \brief Friction law
	 */
	double supportable_tanforce = 0;
	double sq_f_tan = f_contact_tan.sq_norm();
	normal_load = f_contact_normal_norm;
	if (sys->cohesion) {
		normal_load += sys->amplitudes.cohesion;
	}
	if (state == 2) {
		// static friction in previous step
		supportable_tanforce = mu_static*normal_load;
	} else {
		// dynamic friction in previous step
		supportable_tanforce = mu_dynamic*normal_load;
	}
	if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
		state = 3; // dynamic friction
		supportable_tanforce = mu_dynamic*normal_load;
	} else {
		state = 2; // static friction
	}
	if (state == 3) {
		// adjust the sliding spring for dynamic friction law
		disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		f_contact_tan = kt_scaled*disp_tan;
	}
	if (sys->rolling_friction) {
		double supportable_rollingforce = mu_rolling*normal_load;
		double sq_f_rolling = f_rolling.sq_norm();
		if (sq_f_rolling > supportable_rollingforce*supportable_rollingforce) {
			disp_rolling *= supportable_rollingforce/sqrt(sq_f_rolling);
			f_rolling = kr_scaled*disp_rolling;
		}
	}
	return;
}

void Contact::frictionlaw_criticalload()
{
	/* Since reduced_gap < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_reduced_gap(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 * 
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm-sys->amplitudes.critical_normal_force; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_contact_tan.reset();
	} else {
		supportable_tanforce *= mu_static;
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

void Contact::frictionlaw_criticalload_mu_inf()
{
	/* Since reduced_gap < 0, f_contact_normal_norm is always positive.
	 * f_contact_normal_norm = -kn_scaled*interaction->get_reduced_gap(); > 0
	 * F_normal = f_contact_normal_norm(positive) + lubforce_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_contact_normal_norm-sys->amplitudes.critical_normal_force; // critical load model.
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

void Contact::frictionlaw_ft_max()
{
 	/**
	   \brief Friction law
	*/
 	double sq_f_tan = f_contact_tan.sq_norm();
 	if (sq_f_tan > ft_max*ft_max) {
 		state = 3; // dynamic friction
 		disp_tan *= ft_max/sqrt(sq_f_tan);
 		f_contact_tan = kt_scaled*disp_tan;
 	} else {
 		state = 2; // static friction
 	}
 	return;
}

void Contact::frictionlaw_coulomb_max()
{
	/**
	 \brief Friction law
	 *
	 * [NOTE]
	 * In the rate controlled simulation, the results look weird.
	 * There should be a bug.
	 *
	 * The calculated forces in calcContactInteraction() seem to be different
	 * between the rate-controlled and stress-controlled simulation.
	 * The current implementation is quite confusing, and should be fixed.
	 */
	double supportable_tanforce = 0;
	normal_load = f_contact_normal_norm;
	if (state == 2) {
		// static friction in previous step
		supportable_tanforce = mu_static*normal_load;
	} else if (state == 3) {
		// dynamic friction in previous step
		supportable_tanforce = mu_dynamic*normal_load;
	} else {
		exit(1);
	}
	if (supportable_tanforce > ft_max) {
	 	supportable_tanforce = ft_max;
	}
	double sq_f_tan = f_contact_tan.sq_norm();
	if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
		state = 3; // dynamic friction
		disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		f_contact_tan = kt_scaled*disp_tan;
	} else {
		state = 2; // static friction
	}
	return;
}

void Contact::addUpContactForceTorque()
{
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

void Contact::calcContactStress()
{
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
		if (state >= 2) {
			contact_stresslet_XF_tan.set(interaction->rvec, f_contact_tan);
		}
	} else {
		contact_stresslet_XF_normal.reset();
		contact_stresslet_XF_tan.reset();
	}
}

double Contact::calcEnergy()
{
	double overlap = -interaction->get_reduced_gap();
	double sq_tan_norm = disp_tan.sq_norm();
	double sq_disp_rolling = disp_rolling.sq_norm();
	double energy = 0.5*kn_scaled*overlap*overlap;
	if (state >= 2) {
		energy += 0.5*kt_scaled*sq_tan_norm;
	}
	if (sys->rolling_friction) {
		energy += 0.5*kr_scaled*sq_disp_rolling;
	}
	return energy;
}
