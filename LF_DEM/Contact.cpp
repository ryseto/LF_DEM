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
			cerr << lub_coeff << ' ' <<  log(lub_coeff) << endl;
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
	if (sys->friction) {
		/* disp_tan is orthogonal to the normal vector.
		 */
		disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
		if (staticfriction) {
			f_contact_tan = kt_scaled*disp_tan;
		} else {
			/* Dynamic friction:
			 * The spring force is set with the new value of mu*F_n.
			 * The total resistance force consists of lubrication force and dashpot.
			 * In the contact model, the dashpot part is considered to set the tangential spring.
			 * As approximation, the previously calculated resistance force is used here.
			 */
			disp_tan = (1/kt_scaled)*(supportable_tanforce*tvec-sys->get_ratio_dashpot_total()*lubforce_tan);
			f_contact_tan = kt_scaled*disp_tan;
		}
	}
}

void
Contact::frictionlaw(){
	interaction->lubrication.calcLubricationForce();
	double supportable_tanforce = mu*(f_contact_normal_norm+interaction->lubrication.get_lubforce_normal());
	double f_test;
	if (staticfriction){
		/* Friction force is the sum of spring force and dashpot force.
		 * The lubrication forces obtained by solving the force balance equations
		 * are considered as dashpot force.
		 */
		dashpot = interaction->lubrication.lubforce_p0;
		dashpot += interaction->lubrication.get_lubforce_normal()*interaction->nvec; // note: nvec directs from i to j.
		// @@@ The sign (+) is correct, but it is counter-intuitive, we should find a clearer way to do this
		f_test_vec = f_contact_tan+dashpot;
		f_test = f_test_vec.norm();
		if (f_test > supportable_tanforce){
			/* If f_est becomes larger than the supportable tangential force,
			 * the state is switched into the dynamics friction.
			 */
			staticfriction = false;
			sys->fout_dfric << previous_f_test << ' ' << previous_supportable_tanforce << endl;
			sys->fout_dfric << f_test << ' ' << supportable_tanforce << endl << endl;
			double log_lub_coeff = sys->get_log_lub_coeff_dynamicfriction();
			interaction->lubrication.setResistanceCoeffTang(log_lub_coeff);
			lubforce_tan = dashpot;
			tvec = f_test_vec/f_test;
		}
	} else {
		/* lubforce_tan is lubforce + dashpot
		 */ 
		lubforce_tan = interaction->lubrication.lubforce_p0;
		lubforce_tan += interaction->lubrication.get_lubforce_normal()*interaction->nvec;
		// @@@ The sign (+) is correct, but it is counter-intuitive, we should find a clearer way to do this
		f_test_vec = lubforce_tan+f_contact_tan;
		f_test = f_test_vec.norm();
		if (f_test < supportable_tanforce) {
			staticfriction = true;
			double log_lub_coeff = sys->get_log_lub_coeff_staticfriction();
			interaction->lubrication.setResistanceCoeffTang(log_lub_coeff);
			dashpot = lubforce_tan;
			if (previous_f_test != 0){
				sys->fout_sfric << previous_f_test << ' ' << previous_supportable_tanforce << endl;
				sys->fout_sfric << f_test << ' ' << supportable_tanforce << endl << endl;
			}
		} else {
			tvec = f_test_vec/f_test;
			//			static int cnt=0;
			//			if (sys->get_shear_strain() > 2) {
			//				if (cnt++ % 100 == 0) {
			//					cout << f_test << ' ' << supportable_tanforce << endl;
			//				}
			//			}
		}
	}
	previous_f_test = f_test;
	previous_supportable_tanforce = supportable_tanforce;
}

void
Contact::calcContactInteractionRelax(){
	f_contact_normal_norm = -kn_scaled*(interaction->get_gap_nondim()-0.02);
	f_contact_normal = -f_contact_normal_norm*interaction->nvec;
	if (sys->friction) {
		/* disp_tan is orthogonal to the normal vector.
		 */
		disp_tan -= dot(disp_tan, interaction->nvec)*interaction->nvec;
		f_contact_tan = kt_scaled*disp_tan;
		if (sys->frictionlaw == 1) {
			applyFrictionLaw_spring();
		} else {
			applyFrictionLaw_spring_dashpot();
		}
	}
}

/*
 * F_tangent < mu F_normal
 * F_normal = spring_force + dashpot
 * F_tangent = spring_force + dashpot
 *
 */

void
Contact::applyFrictionLaw_spring(){
	/* [NOTE]
	 * Test forces for the friction law include dashpot contributions in Luding model[1].
	 * However, in our overdamped formulation, we use only springs for the test forces.
	 * This approximated friction-law was used for our PRL2013 paper.
	 *
	 * [1] S. Luding. Cohesive, frictional powders: contact models for tension. Granular Matter, 10:235–246, 2008.
	 */
	double f_static = mu*f_contact_normal_norm;
	double sq_f_tan = f_contact_tan.sq_norm();
	if (sq_f_tan > f_static*f_static) {
		/*
		 * The static and dynamic friction coeffients are the same.
		 *
		 */
		disp_tan *= f_static/sqrt(sq_f_tan);
		f_contact_tan = kt_scaled*disp_tan;
		cnt_sliding++; // for output
	}
}

void
Contact::applyFrictionLaw_spring_dashpot(){
	/* [NOTE]
	 * Luding model[1] is modified to simulate static/dynamic friction of particles in fluid.
	 *
	 * [1] S. Luding. Cohesive, frictional powders: contact models for tension. Granular Matter, 10:235–246, 2008.
	 *
	 */
	interaction->lubrication.calcLubricationForce();
	supportable_tanforce = mu*(f_contact_normal_norm+interaction->lubrication.get_lubforce_normal());
	if (staticfriction == false){
		/* dynamic (sliding)
		 *
		 *  F = k*xi + F_lub
		 *  after update of xi,
		 *  If F > mu*Fn, keep the sliding state by updating xi = Fc_max / kt
		 *  If F < mu*Fn, switch into the no-sliding state
		 */
		vec3d lubforce_tan = interaction->lubrication.lubforce_p0;
		lubforce_tan += interaction->lubrication.get_lubforce_normal()*interaction->nvec;
		// @@@ The sign (+) is correct, but it is counter-intuitive, we should find a clearer way to do this
		double lubforce_tan_norm = lubforce_tan.norm();
		vec3d tvec = lubforce_tan/lubforce_tan_norm;
		double ftest = (lubforce_tan+f_contact_tan).norm();
		if (ftest < supportable_tanforce) {
			staticfriction = true;
		} else {
			disp_tan = (1/kt_scaled)*supportable_tanforce*tvec;
		}
	} else {
		/* static (non-sliding)
		 * In the non-sliding state of a contact point, the expected relative velocity is zero.
		 * But, it is finite in simulation.
		 * Thus, the total force is the sum of spring force and dashpot force.
		 * Test force for the friction law is here the sum of them.
		 * If tangential force is smaller than the supportable tangential force,
		 * it kepts the static friction.
		 */
		vec3d dashpot = interaction->lubrication.lubforce_p0;
		dashpot += interaction->lubrication.get_lubforce_normal()*interaction->nvec;
		// @@@ The sign (+) is correct, but it is counter-intuitive, we should find a clearer way to do this
		double f_test = (f_contact_tan+dashpot).norm();
		if (f_test > supportable_tanforce){
			/* If the tangential force becomes larger than the supportable tangential force,
			 * the state is switched into the dynamics friction.
			 * In the dynamic friction, the frictional force is given from only spring.
			 * So, we need to set disp_tan herew
			 */
			vec3d tvec = dashpot/dashpot.norm();
			disp_tan = (1/kt_scaled)*f_test*tvec;
			staticfriction = false;
		}
	}
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

