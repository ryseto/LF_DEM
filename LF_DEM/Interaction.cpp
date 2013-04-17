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
	contact = false;
	active = false;
}

void
Interaction::r(const double &new_r){
	_r = new_r;
	_gap_nondim = _r/ro_2-2; // = h/ro_2
	if (_gap_nondim > 0) {
		lub_coeff = 1/(_gap_nondim+sys->lub_reduce_parameter);
	}
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */

void
Interaction::calcDistanceNormalVector(){
	r_vec = sys->position[par_num[1]]-sys->position[par_num[0]];
	sys->periodize_diff(r_vec, zshift);
	r(r_vec.norm());
	nr_vec = r_vec/_r;
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
	ro = a0+a1;
	ro_2 = ro/2;
	r_lub_max = ro_2*sys->lub_max;
	kn_scaled = ro_2*ro_2*sys->kn; // F = kn_scaled * _gap_nondim;  <-- gap is scaled
	kt_scaled = ro_2*sys->kt; // F = kt_scaled * disp_tan <-- disp is not scaled
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
	colloidal_force_amplitude = sys->cf_amp_dl*a0*a1/ro;
	lambda = a1/a0;
	invlambda = 1/lambda;
	calcDistanceNormalVector();
	if (_gap_nondim <= 0) {
		activate_contact();
	} else {
		contact = false;
	}
	cnt_sliding = 0;
	strain_lub_start = sys->strain(); // for output
	duration_contact = 0; // for output
}

void
Interaction::deactivate(){
	// r > lub_max
#ifdef RECORD_HISTORY
	gap_history.clear();
#endif
	outputSummary();
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
	disp_tan.reset();
	strain_contact_start = sys->strain();
	lub_coeff = sys->lub_coeff_contact;
}

void
Interaction::deactivate_contact(){
	// r > a0 + a1
#ifdef RECORD_HISTORY
	outputHistory();
#endif
	contact = false;
	disp_tan.reset();
	Fc_normal_norm = 0;
	Fc_normal.reset();
	Fc_tan.reset();
	duration_contact += sys->strain()-strain_contact_start; // for output
}

/*
 *
 */
void
Interaction::updateState(bool &deactivated){
	deactivated = false;
	if (active == false) {
		return;
	}
	/* update tangential displacement: we do it before updating nr_vec
	 * as it should be along the tangential vector defined in the previous time step
	 */
	if (contact) {
		if (sys->friction) {
			calcContactVelocity();
			if (sys->in_predictor) {
				disp_tan += contact_velocity*sys->dt;
				disp_tan_predictor = disp_tan;
			} else {
				disp_tan = disp_tan_predictor+contact_velocity*sys->dt;
			}
		}
		calcDistanceNormalVector();
		calcContactInteraction();
		if (sys->colloidalforce) {
			/* For continuity, the colloidal force is kept as constant for h < 0.
			 * This force does not affect the friction law,
			 * i.e. it is separated from Fc_normal_norm.
			 */
			F_colloidal_norm = colloidal_force_amplitude;
			F_colloidal = -F_colloidal_norm*nr_vec;
		}
		if (sys->in_corrector) {
			/* Keep the contact state in predictor.
			 */
			if (_gap_nondim > 0) {
				deactivate_contact();
			}
		}
	} else {
		calcDistanceNormalVector();
		if (sys->colloidalforce) {
			F_colloidal_norm = colloidal_force_amplitude*exp(-(_r-ro)/sys->cf_range_dl);
			F_colloidal = -F_colloidal_norm*nr_vec;
		}
		if (sys->in_corrector) {
			/* If r > r_lub_max, deactivate the interaction object.
			 */
			if (_gap_nondim <= 0) {
				activate_contact();
			} else if (_r > r_lub_max) {
				deactivate();
				deactivated = true;
			}
		}
	}
#ifdef RECORD_HISTORY
	if (!sys->in_predictor) {
		gap_history.push_back(_gap_nondim);
		if (contact) {
			disp_tan_sq_history.push_back(disp_tan.sq_norm());
			overlap_history.push_back(-_gap_nondim);
		}
	}
#endif
}

/*********************************
 *                                *
 *	   Contact Forces Methods    *
 *                                *
 *********************************/
/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * r_vec = p[1] - p[0]
 * Fc_normal_norm is positive (for overlapping particles r < ro)
 */
void
Interaction::calcContactInteraction(){
	Fc_normal_norm = -kn_scaled*_gap_nondim;
	Fc_normal = -Fc_normal_norm*nr_vec;
	if (sys->friction) {
		/* disp_tan is orthogonal to the normal vector.
		 */
		disp_tan -= dot(disp_tan, nr_vec)*nr_vec;
		Fc_tan = kt_scaled*disp_tan;
		checkBreakupStaticFriction();
	}
}

void
Interaction::checkBreakupStaticFriction(){
	double f_static = sys->mu_static*Fc_normal_norm;
	double sq_f_tan = Fc_tan.sq_norm();
	if (sq_f_tan > f_static*f_static) {
		/*
		 * The static and dynamic friction coeffients are the same.
		 *
		 */
		disp_tan *= f_static/sqrt(sq_f_tan);
		Fc_tan = kt_scaled*disp_tan;
		cnt_sliding++; // for output
	}
}

void
Interaction::addUpContactForceTorque(){
	if (contact) {
		sys->contact_force[par_num[0]] += Fc_normal;
		sys->contact_force[par_num[1]] -= Fc_normal;
		if (sys->friction) {
			sys->contact_force[par_num[0]] += Fc_tan;
			sys->contact_force[par_num[1]] -= Fc_tan;
			vec3d t_ij = cross(nr_vec, Fc_tan);
			sys->contact_torque[par_num[0]] += a0*t_ij;
			sys->contact_torque[par_num[1]] += a1*t_ij;
		}
	}
}

/*
 * Colloidal stabilizing force
 *
 *
 */
void
Interaction::addUpColloidalForce(){
	sys->colloidal_force[par_num[0]] += F_colloidal;
	sys->colloidal_force[par_num[1]] -= F_colloidal;
}

/* Relative velocity of particle 1 from particle 0.
 *
 * Use:
 *  sys->velocity and ang_velocity
 *
 */
void
Interaction::calcContactVelocity(){
	// relative velocity particle 1 from particle 0.
	/*
	 * v1' = v1 - Lz = v1 - zshift*lz;
	 */
	/**** NOTE ********************************************
	 * In the Corrector, this contact_velocity
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
	contact_velocity = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
	if (sys->in_predictor && zshift != 0) {
		contact_velocity.x += zshift*sys->vel_difference;
	}
	contact_velocity -=
	cross(a0*sys->ang_velocity[par_num[0]]+a1*sys->ang_velocity[par_num[1]], nr_vec);
}

/*********************************
 *                                *
 *  Lubrication Forces Methods    *
 *                                *
 *********************************/

// Resistance functions
void
Interaction::calcXA(){
	double l1 = 1+lambda;
	double l13 = l1*l1*l1;
	double g1_l = 2*lambda*lambda/l13;
	XA[0] = g1_l*lub_coeff;
	XA[1] = -2*XA[0]/l1;
	XA[2] = XA[1];
	XA[3] = XA[0]/lambda;
}

void
Interaction::calcXG(){
	double l1 = 1+lambda;// defined in XA
	double l13 = l1*l1*l1; // defined in XA
	double g1_l = 2*lambda*lambda/l13; // defined in XA
	double il1 = invlambda+1;
	XG[0] = 1.5*g1_l*lub_coeff;
	XG[1] = -4*XG[0]/l1/l1;
	XG[3] = -XG[0]/lambda;
	XG[2] = -4*XG[3]/il1/il1;
}

void
Interaction::calcXM(){
	double l1 = 1+lambda;
	double l13 = l1*l1*l1;
	double g1_l = 2*lambda*lambda/l13;
	XM[0] = 0.6*g1_l*lub_coeff;
	XM[1] = 8*lambda/l13*XM[0];
	XM[2] = XM[1];
	XM[3] = XM[0]/lambda;
}

void
Interaction::GE(double *GEi, double *GEj){
	calcXG();
	double nxnz_sr = nr_vec.x*nr_vec.z;
	double common_factor_1 = (4*a0*a0*XG[0]+ro*ro*XG[2])/6;
	GEi[0] = common_factor_1*nxnz_sr*nr_vec.x;
	GEi[1] = common_factor_1*nxnz_sr*nr_vec.y;
	GEi[2] = common_factor_1*nxnz_sr*nr_vec.z;
	double common_factor_2 = (4*a1*a1*XG[3]+ro*ro*XG[1])/6;
	GEj[0] = common_factor_2*nxnz_sr*nr_vec.x;
	GEj[1] = common_factor_2*nxnz_sr*nr_vec.y;
	GEj[2] = common_factor_2*nxnz_sr*nr_vec.z;
}

// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jeffrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
void
Interaction::pairVelocityStresslet(const vec3d &vi, const vec3d &vj,
								   stresslet &stresslet_i, stresslet &stresslet_j){
	calcXG();
	double n0n0_13 = nr_vec.x*nr_vec.x-1./3;
	double n1n1_13 = nr_vec.y*nr_vec.y-1./3;
	double n0n1 = nr_vec.x*nr_vec.y;
	double n0n2 = nr_vec.x*nr_vec.z;
	double n1n2 = nr_vec.y*nr_vec.z;
	double common_factor_i = -dot(nr_vec, (4*a0*a0*XG[0]*vi+ro*ro*XG[1]*vj)/6);
	double common_factor_j = -dot(nr_vec, (4*a1*a1*XG[3]*vj+ro*ro*XG[2]*vi)/6);
	stresslet_i.set(n0n0_13, n0n1, n0n2, n1n2, n1n1_13);
	stresslet_i *= common_factor_i;
	stresslet_j.set(n0n0_13, n0n1, n0n2, n1n2, n1n1_13);
	stresslet_j *= common_factor_j;
}

// convenient interface for pairVelocityStresslet(const vec3d &vi, const vec3d &vj, stresslet &stresslet_i, stresslet &stresslet_j)
void
Interaction::pairVelocityStresslet(double* &vel_array, stresslet &stresslet_i, stresslet &stresslet_j){
	vec3d vi, vj;
	int i3 = 3*par_num[0];
	int j3 = 3*par_num[1];
	vi.x = vel_array[i3];
	vi.y = vel_array[i3+1];
	vi.z = vel_array[i3+2];
	vj.x = vel_array[j3];
	vj.y = vel_array[j3+1];
	vj.z = vel_array[j3+2];
	pairVelocityStresslet(vi, vj, stresslet_i, stresslet_j);
}

void
Interaction::pairStrainStresslet(stresslet &stresslet_i, stresslet &stresslet_j){
	double n0n0_13 = nr_vec.x*nr_vec.x-1./3;
	double n1n1_13 = nr_vec.y*nr_vec.y-1./3;
	double n0n1 = nr_vec.x*nr_vec.y;
	double n0n2 = nr_vec.x*nr_vec.z;
	double n1n2 = nr_vec.y*nr_vec.z;
	double rororo = ro*ro*ro;
	calcXM();
	double common_factor_i = 5*(a0*a0*a0*XM[0]/3+rororo*XM[1]/24)*n0n2;
	stresslet_i.set(n0n0_13, n0n1, n0n2, n1n2, n1n1_13);
	stresslet_i *= common_factor_i;
	double common_factor_j = 5*(a1*a1*a1*XM[3]/3+rororo*XM[2]/24)*n0n2;
	stresslet_j.set(n0n0_13, n0n1, n0n2, n1n2, n1n1_13);
	stresslet_j *= common_factor_j;
}

void
Interaction::addHydroStress(){
	int i3 = 3*par_num[0];
	int j3 = 3*par_num[1];
	stresslet stresslet_GU_i;
	stresslet stresslet_GU_j;
	stresslet stresslet_ME_i;
	stresslet stresslet_ME_j;
	vec3d vi(sys->v_hydro[i3], sys->v_hydro[i3+1], sys->v_hydro[i3+2]);
	vec3d vj(sys->v_hydro[j3], sys->v_hydro[j3+1], sys->v_hydro[j3+2]);
	/*
	 *  First: -G*(U-Uinf) term
	 */
	pairVelocityStresslet(vi, vj, stresslet_GU_i, stresslet_GU_j);
	/*
	 *  Second: +M*Einf term
	 */
	pairStrainStresslet(stresslet_ME_i, stresslet_ME_j);
	sys->lubstress[par_num[0]] += stresslet_GU_i+stresslet_ME_i;
	sys->lubstress[par_num[1]] += stresslet_GU_j+stresslet_ME_j;
	//		lubstresslet.elm[u] = \
	stresslet_GU_i.elm[u]+stresslet_ME_i.elm[u]+stresslet_GU_j.elm[u]+stresslet_ME_j.elm[u];
}

/* Lubriction force between two particles is calculated.
 * Note that only the Brownian component of the velocity is NOT included here.
 * This part is used for ouput data.
 * lubforce_j = -lubforce_i
 */
void
Interaction::evaluateLubricationForce(){
	/*
	 *  First: -A*(U-Uinf) term
	 */
	vec3d vi = sys->velocity[par_num[0]];
	vec3d vj = sys->velocity[par_num[1]];
	if (sys->vel_difference > 0) {
		vi.x -= sys->position[par_num[0]].z;
		vj.x -= sys->position[par_num[1]].z;
	}
	calcXA();
	double cf_AU_i = -dot(a0*XA[0]*vi+0.5*ro*XA[1]*vj, nr_vec);
	/*
	 *  Second -tildeG*(-Einf)term
	 */
	calcXG();
	double cf_GE_i = nr_vec.x*nr_vec.z*(a0*a0*XG[0]*4+ro*ro*XG[2])/6;
	lubforce_i = (cf_AU_i+cf_GE_i)*nr_vec;
}

void
Interaction::addContactStress(){
	if (contact) {
		int i3 = 3*par_num[0];
		int j3 = 3*par_num[1];
		
		/*
		 * Fc_normal_norm = -kn_scaled*_gap_nondim; --> positive
		 * Fc_normal = -Fc_normal_norm*nr_vec;
		 * This force acts on particle 1.
		 * stress1 is a0*nr_vec[*]force.
		 * stress2 is (-a1*nr_vec)[*](-force) = a1*nr_vec[*]force
		 */
		vec3d contact_force = Fc_normal+Fc_tan;
		stresslet contactstressletXF(nr_vec, contact_force);
		sys->contactstressXF[par_num[0]] += a0*contactstressletXF;
		sys->contactstressXF[par_num[1]] += a1*contactstressletXF;
		// Add term G*V_cont
		stresslet stresslet_GU_i;
		stresslet stresslet_GU_j;
		vec3d vi(sys->v_cont[i3], sys->v_cont[i3+1], sys->v_cont[i3+2]);
		vec3d vj(sys->v_cont[j3], sys->v_cont[j3+1], sys->v_cont[j3+2]);
		pairVelocityStresslet(vi, vj, stresslet_GU_i, stresslet_GU_j);
		sys->contactstressGU[par_num[0]] += stresslet_GU_i;
		sys->contactstressGU[par_num[1]] += stresslet_GU_j;
	}
}

void
Interaction::addColloidalStress(){
	int i3 = 3*par_num[0];
	int j3 = 3*par_num[1];
	colloidalstressletXF.set(r_vec, F_colloidal);
	// Add term G*V_cont
	stresslet stresslet_colloid_GU_i;
	stresslet stresslet_colloid_GU_j;
	vec3d vi(sys->v_colloidal[i3], sys->v_colloidal[i3+1], sys->v_colloidal[i3+2]);
	vec3d vj(sys->v_colloidal[j3], sys->v_colloidal[j3+1], sys->v_colloidal[j3+2]);
	pairVelocityStresslet(vi, vj, stresslet_colloid_GU_i, stresslet_colloid_GU_j);
	sys->colloidalstressGU[par_num[0]] += stresslet_colloid_GU_i;
	sys->colloidalstressGU[par_num[1]] += stresslet_colloid_GU_j;
}

#ifdef RECORD_HISTORY
void
Interaction::outputHistory(){
	cerr << "cnt_sliding =" << cnt_sliding << ' '  << endl;
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
	duration = sys->strain()-strain_lub_start;
	sys->fout_int_data << strain_lub_start << ' '; // 1
	sys->fout_int_data << duration << ' '; // 2
	sys->fout_int_data << duration-duration_contact << ' '; // 3
	sys->fout_int_data << duration_contact << ' '; // 4
	sys->fout_int_data << cnt_sliding << ' '; //  7
	sys->fout_int_data << endl;
}

double
Interaction::getContactVelocity(){
	if (contact == false) {
		return 0;
	}
	sys->in_predictor = true;
	calcDistanceNormalVector();
	calcContactVelocity();
	return contact_velocity.norm();
}

double
Interaction::getNormalVelocity(){
	sys->in_predictor = true;
	calcDistanceNormalVector();
	vec3d d_velocity = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
	if (zshift != 0) {
		d_velocity.x += zshift*sys->vel_difference;
	}
	return dot(d_velocity, nr_vec);
}

double
Interaction::getPotentialEnergy(){
	double energy;
	//	double h = _r - ro;
	if (_gap_nondim < 0) {
		energy = 0.5*sys->kn*_gap_nondim*_gap_nondim;
		energy += -colloidal_force_amplitude*_gap_nondim;
	} else {
		energy = sys->cf_range_dl*colloidal_force_amplitude*(exp(-(_r-ro)/sys->cf_range_dl)-1);
	}
	return energy;
}
