//
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Contact.h"
#include "Interaction.h"
#include "System.h"

void Contact::init(System* sys_, Interaction* interaction_)
{
	sys = sys_;
	interaction = interaction_;
	dashpot.init(sys_, interaction_);
	state = 0;
	f_spring_normal_norm = 0;
	f_spring_normal.reset();
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
    double ro_12 = interaction->ro_12;
    kn_scaled = ro_12*ro_12*sys->p.kn; // F = kn_scaled * _reduced_gap;  <-- gap is scaled @@@@ Why use reduced_gap? Why not gap?
    if (sys->friction) {
        kt_scaled = ro_12*sys->p.kt; // F = kt_scaled * disp_tan <-- disp is not scaled
        if (sys->rolling_friction) {
            kr_scaled = ro_12*sys->p.kr; // F = kt_scaled * disp_tan <-- disp is not scaled
        }
    }
}

void Contact::setInteractionData()
{
	std::tie(p0, p1) = interaction->get_par_num();
	setSpringConstants();
	if (sys->friction) {
		mu_static = sys->p.mu_static;
		mu_dynamic = sys->p.mu_dynamic;
		if (sys->rolling_friction) {
			mu_rolling = sys->p.mu_rolling;
		}
	}
	dashpot.setParticleData();
	dashpot.setDashpotResistanceCoeffs(kn_scaled, kt_scaled,
		                                 sys->p.contact_relaxation_time, sys->p.contact_relaxation_time_tan);
}

void Contact::activate()
{
	// r < a0 + a1
	/* state = 1 means frictionless contact.
	 * In frictional particle simulations,
	 * this value will be updated to 2 or 3 in friction law.
	 * In critical load model, the value can take 1 as well.
	 */
	f_spring_normal_norm = 0;
	f_spring_normal.reset();
	normal_load = 0;
	f_spring_total.reset();
	f_contact_total.reset();
	if (sys->friction) {
		if (sys->p.friction_model == 2 || sys->p.friction_model == 3) {
			state = 1; // critical load model
		} else {
			state = 2; // static friction
		}
		disp_tan.reset();
		disp_rolling.reset();
		prev_disp_tan.reset();
		prev_disp_rolling.reset();
		f_rolling.reset();
	} else {
		state = 1;
	}
	dashpot.activate();
	sys->updateNumberOfContacts(p0, p1, 1);
}

void Contact::deactivate()
{
	// r > a0 + a1
	state = 0;
	f_spring_normal_norm = 0;
	f_spring_normal.reset();
	f_spring_total.reset();
	f_contact_total.reset();
	if (sys->friction) {
		disp_tan.reset();
		disp_rolling.reset();
		f_spring_tan.reset();
	}
	dashpot.deactivate();
	sys->updateNumberOfContacts(p0, p1, -1);
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

void Contact::calcContactSpringForce()
{
	/**
		\brief Compute the contact forces and apply friction law, by rescaling forces and stretches if necessary.

	 f_something is the force acting on the particle 0.
	 disp_something is the relative displacement of the particle 1 from the particle 0.
	 Therefore, the sign of force is same as the one of the displacement.
	 */

	/* h < 0
	 * f_spring_normal_norm > 0 ..... repulsive force
	 * h > 0
	 * f_spring_normal_norm < 0 ..... attractive force
	 */
	f_spring_normal_norm = -kn_scaled*interaction->get_reduced_gap();
	f_spring_normal = -f_spring_normal_norm*interaction->nvec;
	if (sys->friction) {
		disp_tan.vertical_projection(interaction->nvec);
		f_spring_tan = kt_scaled*disp_tan;
		if (sys->rolling_friction) {
			f_rolling = kr_scaled*disp_rolling;
		}
		(this->*frictionlaw)();
	}
	//	calcScaledForce();
}

void Contact::calcTotalForce()
{
	/**
		\brief Compute the total contact forces (spring+dashpot).
		!!! Needs to have the correct velocities in the System class!
		(This contains the dashpot, it is NOT only a static force.)
		*/
	calcContactSpringForce();
	f_contact_total = f_spring_normal + f_spring_tan + dashpot.getPairwiseForce();
}

void Contact::frictionlaw_standard()
{
	/**
	 \brief Friction law
	 */
	double supportable_tanforce = 0;
	double sq_f_tan = f_spring_tan.sq_norm();
	normal_load = f_spring_normal_norm;
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
		f_spring_tan = kt_scaled*disp_tan;
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
	/* Since reduced_gap < 0, f_spring_normal_norm is always positive.
	 * f_spring_normal_norm = -kn_scaled*interaction->get_reduced_gap(); > 0
	 * F_normal = f_spring_normal_norm(positive) + lubforce_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_spring_normal_norm-sys->amplitudes.critical_normal_force; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_spring_tan.reset();
	} else {
		supportable_tanforce *= mu_static;
		double sq_f_tan = f_spring_tan.sq_norm();
		if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
			state = 3; // sliding
			disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
			f_spring_tan = kt_scaled*disp_tan;
		} else {
			state = 2; // static friction
		}
	}
	return;
}

void Contact::frictionlaw_criticalload_mu_inf()
{
	/* Since reduced_gap < 0, f_spring_normal_norm is always positive.
	 * f_spring_normal_norm = -kn_scaled*interaction->get_reduced_gap(); > 0
	 * F_normal = f_spring_normal_norm(positive) + lubforce_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_spring_normal_norm-sys->amplitudes.critical_normal_force; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_spring_tan.reset();
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
 	double sq_f_tan = f_spring_tan.sq_norm();
 	if (sq_f_tan > ft_max*ft_max) {
 		state = 3; // dynamic friction
 		disp_tan *= ft_max/sqrt(sq_f_tan);
 		f_spring_tan = kt_scaled*disp_tan;
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
	normal_load = f_spring_normal_norm;
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
	double sq_f_tan = f_spring_tan.sq_norm();
	if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
		state = 3; // dynamic friction
		disp_tan *= supportable_tanforce/sqrt(sq_f_tan);
		f_spring_tan = kt_scaled*disp_tan;
	} else {
		state = 2; // static friction
	}
	return;
}

void Contact::addUpContactForceTorque()
{
    /* Force
	 */
	if (state <= 1) {
    f_spring_total = f_spring_normal;
	} else {
		f_spring_total = f_spring_normal+f_spring_tan;
	}
	sys->contact_force[p0] += f_spring_total;
	sys->contact_force[p1] -= f_spring_total;
	/* Torque
	 */
	if (state >= 2) {
		vec3d t_ij = cross(interaction->nvec, f_spring_tan);
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
	/**
	 * The "xF" contact stress.
	 * The contact force F includes both the spring force and the dashpot force.
	 */
	if (is_active() > 0) {
		contact_stresslet_XF.set(interaction->rvec, f_spring_total + dashpot.getPairwiseForce());
	} else {
		contact_stresslet_XF.reset();
	}
}

double Contact::calcEnergy()
{
	double overlap = -interaction->get_reduced_gap();
	double sq_tan_norm = disp_tan.sq_norm();
	double sq_disp_rolling = disp_rolling.sq_norm();
	/* normal */
	double energy = 0.5*kn_scaled*overlap*overlap;
	if (state >= 2) {
			/* sliding */
		energy += 0.5*kt_scaled*sq_tan_norm;
		if (sys->rolling_friction) {
			/* roling */
			energy += 0.5*kr_scaled*sq_disp_rolling;
		}
	}
	return energy;
}

double Contact::get_f_normal_norm()
{
	return abs(dot(interaction->nvec, f_contact_total));
}

double Contact::get_normal_load()
{
	return normal_load;
}

vec3d Contact::get_f_tan()
{
	return f_contact_total - dot(interaction->nvec, f_contact_total)*interaction->nvec;
}

double Contact::get_f_tan_norm()
{
	return get_f_tan().norm();
}
