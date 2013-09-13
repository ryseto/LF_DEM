//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Interaction.h"

void
Interaction::init(System *sys_){
	sys = sys_;
	active = false;
	lubrication.init(sys);
	contact.init(sys, this);
}


void 
Interaction::setResistanceCoeff(double normal_rc, double tangent_rc){
	lubrication.setResistanceCoeff(normal_rc, tangent_rc);
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */
void
Interaction::calcNormalVectorDistanceGap(){
	rvec = sys->position[par_num[1]]-sys->position[par_num[0]];
	sys->periodize_diff(rvec, zshift);
	r = rvec.norm();
	nvec = rvec/r;
	nxnx = nvec.x*nvec.x;
	nxny = nvec.x*nvec.y;
	nxnz = nvec.x*nvec.z;
	nynz = nvec.y*nvec.z;
	nyny = nvec.y*nvec.y;
	nznz = nvec.z*nvec.z;
	gap_nondim = r/ro_12-2;
	if (contact.active) {
		double overlap_12 = 0.5*(a0+a1-r);
		a0_dash = a0-overlap_12;
		a1_dash = a1-overlap_12;
	}
	else{
		double lub_coeff = 1/(gap_nondim+sys->lub_reduce_parameter);
		lubrication.setResistanceCoeff(lub_coeff, log(lub_coeff));
		a0_dash = a0;
		a1_dash = a1;
	}


}

/* Activate interaction between particles i and j.
 * Always j>i is satisfied.
 */
void
Interaction::activate(int i, int j){
	active = true;
	if (j > i) {
		par_num[0] = i;
		par_num[1] = j;
	} else {
		par_num[0] = j;
		par_num[1] = i;
	}
	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);
	// tell them their new partner
	sys->interaction_partners[i].insert(j);
	sys->interaction_partners[j].insert(i);
	a0 = sys->radius[par_num[0]];
	a1 = sys->radius[par_num[1]];
	set_ro(a0+a1); // ro=a0+a1
	interaction_range_scaled = ro_12*sys->get_lub_max();
	/* NOTE:
	 * lub_coeff_contact includes kn.
	 * If the scaled kn is used there,
	 * particle size dependence appears in the simulation.
	 * I don't understand this point yet.
	 * lub_coeff_contact_scaled = 4*kn_scaled*sys->contact_relaxzation_time;
	 */
	/*
	 * The size dependence of colloidal force:
	 * a0*a1/(a1+a2)/2
	 * Is
	 */
	colloidalforce_amplitude = sys->get_colloidalforce_amplitude()*a0*a1/ro;
	colloidalforce_length = sys->get_colloidalforce_length();
	calcNormalVectorDistanceGap();

	// deal with contact
	contact.getInteractionData();
	if (gap_nondim <= 0) {
		contact.activate();
	}
	else{
		contact.deactivate();
	}
	contact.resetObservables();

	lubrication.getInteractionData();
	strain_lub_start = sys->get_shear_strain(); // for output
	lubrication.calcLubConstants();

void
Interaction::updateContactModel(){
	if (active) {
		kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
		kt_scaled = ro_12*sys->get_kt(); // F = kt_s
		if (contact) {
			lub_coeff = sys->get_lub_coeff_contact();
			log_lub_coeff = sys->get_log_lub_coeff_contact();
		}
	}
}

void
Interaction::deactivate(){
	// r > lub_max
#ifdef RECORD_HISTORY
	gap_history.clear();
#endif
	outputSummary();
	contact.deactivate();
	active = false;
	sys->interaction_list[par_num[0]].erase(this);
	sys->interaction_list[par_num[1]].erase(this);
	sys->interaction_partners[par_num[0]].erase(par_num[1]);
	sys->interaction_partners[par_num[1]].erase(par_num[0]);
}

void
Interaction::activate_contact(){
	// r < a0 + a1
	contact = true;
	staticfriction = false;
	disp_tan.reset();
	strain_contact_start = sys->get_shear_strain();
	lub_coeff = sys->get_lub_coeff_contact();
	log_lub_coeff = sys->get_log_lub_coeff_contact();

}

void
Interaction::deactivate_contact(){
	// r > a0 + a1
#ifdef RECORD_HISTORY
	outputHistory();
#endif
	contact = false;
	disp_tan.reset();
	f_contact_normal_norm = 0;
	f_contact_normal.reset();
	f_contact_tan.reset();
	duration_contact += sys->get_shear_strain()-strain_contact_start; // for output
}

void
Interaction::updateState(bool &deactivated){
	/* update tangential displacement: we do it before updating nvec
	 * as it should be along the tangential vector defined in the previous time step
	 */
	if (contact.active) {
		if (sys->friction) {
			calcRelativeVelocities();
			contact.incrementTangentialDisplacement();
		}
		calcNormalVectorDistanceGap();
		contact.calcContactInteraction();
		if (sys->colloidalforce) {
			/* For continuity, the colloidal force is kept as constant for h < 0.
			 * This force does not affect the friction law,
			 * i.e. it is separated from Fc_normal_norm.
			 */
			f_colloidal_norm = colloidalforce_amplitude;
			f_colloidal = -f_colloidal_norm*nvec;
		}
		if (sys->in_corrector) {
			/* Keep the contact state in predictor.
			 */
			if (gap_nondim > 0) {
				contact.deactivate();
			}
		}
	} else { /* separating */
		calcNormalVectorDistanceGap();
		if (sys->colloidalforce) {
			f_colloidal_norm = colloidalforce_amplitude*exp(-(r-ro)/colloidalforce_length);
			f_colloidal = -f_colloidal_norm*nvec;
		}
		if (sys->in_corrector) {
			/* If r > interaction_range_scaled, deactivate the interaction object.
			 */
			if (gap_nondim <= 0) {
				contact.activate();
			} else if (r > interaction_range_scaled) {
				deactivate();
				deactivated = true;
			}
		}
	}
#ifdef RECORD_HISTORY
	if (!sys->in_predictor) {
		gap_history.push_back(gap_nondim);
		if (contact) {
			disp_tan_sq_history.push_back(disp_tan.sq_norm());
			overlap_history.push_back(-gap_nondim);
		}
	}
#endif
}


void
Interaction::updateStateRelax(bool &deactivated){
	deactivated = false;
	if (active == false) {
		return;
	}
	/* update tangential displacement: we do it before updating nvec
	 * as it should be along the tangential vector defined in the previous time step
	 */
	calcNormalVectorDistanceGap();
	if (contact.active) {
		contact.calcContactInteractionRelax();
		if (gap_nondim > 0) {
			contact.deactivate();
		}
		f_colloidal_norm = colloidalforce_amplitude;
		f_colloidal = -f_colloidal_norm*nvec;
	} else {
		f_colloidal_norm = colloidalforce_amplitude*exp(-(r-ro)/colloidalforce_length);
		f_colloidal = -f_colloidal_norm*nvec;
		if (gap_nondim <= 0) {
			contact.activate();
		} else if (r > interaction_range_scaled) {
			deactivate();
			deactivated = true;
		}
	}
}


/*********************************
 *                                *
 *	   Contact Forces Methods    *
 *                                *
 *********************************/
/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * rvec = p[1] - p[0]
 * Fc_normal_norm is positive (for overlapping particles r < ro)
 */
void
Interaction::calcContactInteraction(){
	// gap_nondim < 0 (in contact)
	if (gap_nondim > 0) {
		return;
	}
	f_contact_normal_norm = -kn_scaled*gap_nondim; // gap_nondim is negative, therefore it is allways positive.
	f_contact_normal = -f_contact_normal_norm*nvec;
	if (sys->friction) {
		/* disp_tan is orthogonal to the normal vector.
		 */
		/*
		 * [possible issue]
		 * Is this projection of disp_tan ok for this place?
		 */
		disp_tan -= dot(disp_tan, nvec)*nvec;
		f_contact_tan = kt_scaled*disp_tan;
		if (sys->frictionlaw == 1) {
			applyFrictionLaw_spring();
		} else {
			applyFrictionLaw_spring_dashpot();
		}
	}
}

/*
 * We need to modify the friction law.
 * F_tangent < mu F_normal
 * The total forces should be used in the friction law.
 * F_normal = spring_force + dashpot
 * F_tangent = spring_force + dashpot
 *
 */

void
Interaction::applyFrictionLaw_spring(){
	/* [NOTE]
	 * Test forces for the friction law include dashpot contributions in Luding model[1].
	 * However, in our overdamped formulation, we use only springs for the test forces.
	 * This approximated friction-law was used for our PRL2013 paper.
	 *
	 * [1] S. Luding. Cohesive, frictional powders: contact models for tension. Granular Matter, 10:235–246, 2008.
	 */
	calcLubricationForce();
	double lubforce_norm = -dot(lubforce_i, nvec);
	double supportable_tangential_force = sys->get_mu_static()*(f_contact_normal_norm+lubforce_norm);
	
	//double f_static = sys->get_mu_static()*f_contact_normal_norm;
//	double sq_f_tan = f_contact_tan.sq_norm();
	double f_test = f_contact_tan.norm();

	if (f_test >= supportable_tangential_force) {
		/*
		 * The static and dynamic friction coeffients are the same.
		 *
		 */
		disp_tan *= supportable_tangential_force/f_test;
		f_contact_tan = kt_scaled*disp_tan;
		cnt_sliding++; // for output
	}
}

void
Interaction::applyFrictionLaw_spring_dashpot(){
	/* [NOTE]
	 * Luding model[1] is modified to simulate static/dynamic friction of particles in fluid.
	 *
	 * [1] S. Luding. Cohesive, frictional powders: contact models for tension. Granular Matter, 10:235–246, 2008.
	 *
	 */
	calcLubricationForce();
	double lubforce_norm = -dot(lubforce_i, nvec);
	double supportable_tangential_force = sys->get_mu_static()*(f_contact_normal_norm+lubforce_norm);
	if (staticfriction == false){
		/* dynamic (sliding)
		 *
		 *  F = Fc_max + F_lub
		 *  Fc_max = mu*Fn
		 *  after update of xi,
		 *  If kt*xi > Fc_max, keep the sliding state.
		 *  If kt*xi < Fc_max, switch into the no-sliding state.
		 */
		vec3d lubforce_tan = lubforce_i+lubforce_norm*nvec;
		double f_test = f_contact_tan.norm();
		double lubforce_tan_norm = lubforce_tan.norm();
		vec3d tvec = lubforce_tan/lubforce_tan_norm;
		if (f_test >= supportable_tangential_force){
			/* ft_max is not enough large to stop sliding.
			 * The force acting on the contact point can be estimated
			 * from the streathed spring. kt*disp_tan.
			 * Then, the sliding state is kept.
			 *
			 */
			disp_tan = (1/kt_scaled)*supportable_tangential_force*tvec;
		} else {
			/* ft_max becomes large or the tangential force acting on the contact point becomes smaller.
			 * So, the possibility becoming the static friction is considered.
			 *
			 */
			if ( f_test+(1-sys->ratio_dashpot_lubrication)*lubforce_tan_norm > 0 ){
				disp_tan = (1/kt_scaled)*(f_test+(1-sys->ratio_dashpot_lubrication)*lubforce_tan_norm)*tvec;
			} else {
				cerr << f_test << ' ' << 1-sys->ratio_dashpot_lubrication << ' ' << lubforce_tan_norm << endl;
				disp_tan = (1/kt_scaled)*supportable_tangential_force*tvec;
			}
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
		vec3d dashpot = lubforce_i+lubforce_norm*nvec;
		if (dot(f_contact_tan, dashpot)<0){
			cerr << "." ;
		}
		double f_test = (f_contact_tan+dashpot).norm();
		if (f_test > supportable_tangential_force){
			/* If the tangential force becomes larger than the supportable tangential force,
			 * the state is switched into the dynamics friction.
			 * In the dynamic friction, the frictional force is given from only spring.
			 * So, we need to set disp_tan here
			 */
			vec3d tvec = dashpot/dashpot.norm();
			disp_tan = (1/kt_scaled)*supportable_tangential_force*tvec;
			staticfriction = false;
		}
	}
}

void
Interaction::addUpContactForceTorque(){
	if (contact) {
		sys->contact_force[par_num[0]] += f_contact_normal;
		sys->contact_force[par_num[1]] -= f_contact_normal;
		if (sys->friction) {
			sys->contact_force[par_num[0]] += f_contact_tan;
			sys->contact_force[par_num[1]] -= f_contact_tan;
			vec3d t_ij = cross(nvec, f_contact_tan);
			sys->contact_torque[par_num[0]] += a0_dash*t_ij;
			sys->contact_torque[par_num[1]] += a1_dash*t_ij;
		}
	}
}

/*
 * Colloidal stabilizing force
 */
void
Interaction::addUpColloidalForce(){
	sys->colloidal_force[par_num[0]] += f_colloidal;
	sys->colloidal_force[par_num[1]] -= f_colloidal;
}

/* Relative velocity of particle 1 from particle 0.
 *
 * Use:
 *  sys->velocity and ang_velocity
 *
 */
void
Interaction::calcRelativeVelocities(){
	/* relative velocity particle 1 from particle 0.
	 */
	/*
	 * v1' = v1 - Lz = v1 - zshift*lz;
	 */
	/**** NOTE ********************************************
	 * In the Corrector, this relative_surface_velocity
	 * is also the correcting velocity.
	 * This correcting velocity should not involve the
	 * velocity diffrence due to crossing the z boundary.
	 * fix_interaction_status = true : in the Predictor
	 * fix_interaction_status = false : in the Corrector
	 *
	 * if p1 is upper, zshift = -1.
	 * zshift = -1; //  p1 (z ~ lz), p0 (z ~ 0)
	 *
	 ******************************************************/
	vec3d dv = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
	if (sys->in_predictor && zshift != 0) {
		dv.x += zshift*sys->vel_difference;
	}
	relative_surface_velocity = dv-cross(a0*sys->ang_velocity[par_num[0]]+a1*sys->ang_velocity[par_num[1]], nvec);
	relative_surface_velocity -= dot(relative_surface_velocity, nvec)*nvec; 
}


void
Interaction::addColloidalStress(){
	colloidal_stresslet_XF.set(rvec, f_colloidal);
}

#ifdef RECORD_HISTORY
void
Interaction::outputHistory(){
	cerr << "cnt_sliding =" << contact.cnt_sliding << ' '  << endl;
	if (sys->strain() > 1 && disp_tan_sq_history.size() > 10) {
		for (int i=0; i<disp_tan_sq_history.size(); i += 10) {
			cout << overlap_history[i] << ' ' << sqrt(disp_tan_sq_history[i]) << endl;
		}
		cout << endl;
		if (sys->cnt_monitored_data++ == 200) {
			exit(1);
		}
	}
	disp_tan_sq_history.clear();
	overlap_history.clear();
}
#endif

void
Interaction::outputSummary(){
	duration = sys->get_shear_strain()-strain_lub_start;
	sys->fout_int_data << strain_lub_start << ' '; // 1
	sys->fout_int_data << duration << ' '; // 2
	sys->fout_int_data << duration-contact.get_duration() << ' '; // 3
	sys->fout_int_data << contact.get_duration() << ' '; // 4
	sys->fout_int_data << contact.get_cnt_sliding() << ' '; //  7
	sys->fout_int_data << endl;
}

double
Interaction::getContactVelocity(){
	if (contact.active == false) {
		return 0;
	}
	return relative_surface_velocity.norm();
}

double
Interaction::getNormalVelocity(){
	sys->in_predictor = true;
	vec3d d_velocity = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
	if (zshift != 0) {
		d_velocity.x += zshift*sys->vel_difference;
	}
	return dot(d_velocity, nvec);
}

double
Interaction::getPotentialEnergy(){
	double energy;
	if (gap_nondim < 0) {
		energy = 0.5*sys->get_kn()*gap_nondim*gap_nondim;
		energy += -colloidalforce_amplitude*gap_nondim;
	} else {
		energy = colloidalforce_length*colloidalforce_amplitude*(exp(-(r-ro)/colloidalforce_length)-1);
	}
	return energy;
}
