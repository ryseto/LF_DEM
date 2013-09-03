//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Interaction.h"
#include "LubricationFunctions.h"

void
Interaction::init(System *sys_){
	sys = sys_;
	contact = false;
	active = false;
	c13 = 0;
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */
void
Interaction::calcNormalVectorDistanceGap(){
	r_vec = sys->position[par_num[1]]-sys->position[par_num[0]];
	sys->periodize_diff(r_vec, zshift);
	r = r_vec.norm();
	nvec = r_vec/r;
	nxnx = nvec.x*nvec.x;
	nxny = nvec.x*nvec.y;
	nxnz = nvec.x*nvec.z;
	nynz = nvec.y*nvec.z;
	nyny = nvec.y*nvec.y;
	nznz = nvec.z*nvec.z;
	gap_nondim = r/ro_12-2;
	if (!contact) {
		lub_coeff = 1/(gap_nondim+sys->lub_reduce_parameter);
		log_lub_coeff = log(lub_coeff);
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
	i6 = 6*par_num[0];
	j6 = 6*par_num[1];
	
	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);
	// tell them their new partner
	sys->interaction_partners[i].insert(j);
	sys->interaction_partners[j].insert(i);
	a0 = sys->radius[par_num[0]];
	a1 = sys->radius[par_num[1]];
	set_ro(a0+a1); // ro=a0+a1
	lub_max_scaled = ro_12*sys->get_lub_max();
	kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
	kt_scaled = ro_12*sys->get_kt(); // F = kt_scaled * disp_tan <-- disp is not scaled
	tangential_dashpot_coeff = ro_12*sys->get_tang_coeff_contact(); // **** to be checked ****
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
	if (gap_nondim <= 0) {
		activate_contact();
	} else {
		contact = false;
	}
	cnt_sliding = 0;
	strain_lub_start = sys->get_shear_strain(); // for output
	duration_contact = 0; // for output
	calcLubConstants();
}

void
Interaction::updateContactModel(){
	if (active) {
		kn_scaled = ro_12*ro_12*sys->get_kn(); // F = kn_scaled * _gap_nondim;  <-- gap is scaled
		kt_scaled = ro_12*sys->get_kt(); // F = kt_s
		if (contact) {
			lub_coeff = sys->get_lub_coeff_contact();
			log_lub_coeff = log(lub_coeff);
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
	active = false;
	contact = false;
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
	strain_contact_start = sys->get_shear_strain();
	lub_coeff = sys->get_lub_coeff_contact();
	log_lub_coeff = log(1/(sys->lub_reduce_parameter));
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
	if (contact) {
		if (sys->friction) {
			calcRelativeVelocities();
			if (sys->in_predictor) {
				disp_tan += contact_velocity*sys->get_dt();
				disp_tan_predictor = disp_tan;
			} else {
				disp_tan = disp_tan_predictor+contact_velocity*sys->get_dt();
			}
		}
		calcNormalVectorDistanceGap();
		calcContactInteraction();
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
				deactivate_contact();
			}
		}
	} else { /* separating */
		calcNormalVectorDistanceGap();
		if (sys->colloidalforce) {
			f_colloidal_norm = colloidalforce_amplitude*exp(-(r-ro)/colloidalforce_length);
			f_colloidal = -f_colloidal_norm*nvec;
		}
		if (sys->in_corrector) {
			/* If r > r_lub_max, deactivate the interaction object.
			 */
			if (gap_nondim <= 0) {
				activate_contact();
			} else if (r > lub_max_scaled) {
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
	if (contact) {
		f_contact_normal_norm = -kn_scaled*(gap_nondim-0.02);
		f_contact_normal = -f_contact_normal_norm*nvec;
		if (gap_nondim > 0) {
			deactivate_contact();
		}
		if (sys->friction) {
			/* disp_tan is orthogonal to the normal vector.
			 */
			disp_tan -= dot(disp_tan, nvec)*nvec;
			f_contact_tan = kt_scaled*disp_tan;
			checkBreakupStaticFriction();
		}
		f_colloidal_norm = colloidalforce_amplitude;
		f_colloidal = -f_colloidal_norm*nvec;
	} else {
		f_colloidal_norm = colloidalforce_amplitude*exp(-(r-ro)/colloidalforce_length);
		f_colloidal = -f_colloidal_norm*nvec;
		if (gap_nondim <= 0) {
			activate_contact();
		} else if (r > lub_max_scaled) {
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
 * r_vec = p[1] - p[0]
 * Fc_normal_norm is positive (for overlapping particles r < ro)
 */
void
Interaction::calcContactInteraction(){
	// gap_nondim < 0 (in contact)
	//
	f_contact_normal_norm = -kn_scaled*gap_nondim; // gap_nondim is negative, therefore it is allways positive.
	
	f_contact_normal = -f_contact_normal_norm*nvec;
	if (sys->friction) {
		/* disp_tan is orthogonal to the normal vector.
		 */
		disp_tan -= dot(disp_tan, nvec)*nvec;
		f_contact_tan = kt_scaled*disp_tan;
		checkBreakupStaticFriction();
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
Interaction::checkBreakupStaticFriction(){
	/* [NOTE]
	 * Test forces for the friction law include dashpot contributions in Luding model[1].
	 * However, in our overdamped formulation, we use only springs for the test forces.
	 * This approximated friction-law was used for our PRL2013 paper.
	 *
	 * [1] S. Luding. Cohesive, frictional powders: contact models for tension. Granular Matter, 10:235–246, 2008.
	 */
	double f_static = sys->get_mu_static()*f_contact_normal_norm;
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
Interaction::addUpContactForceTorque(){
	if (contact) {
		
		//		cerr << " I am interaction " <<endl;
		//		info();
		
		sys->contact_force[par_num[0]] += f_contact_normal;
		/*		if(f_contact_normal < 0){
			cout << "neg force " << f_contact_normal<< " " << nvec <<" " << sys->position[par_num[0]] << " " << sys->position[par_num[1]] << endl;
			getchar();
			}*/
		sys->contact_force[par_num[1]] -= f_contact_normal;
		if (sys->friction) {
			sys->contact_force[par_num[0]] += f_contact_tan;
			sys->contact_force[par_num[1]] -= f_contact_tan;
			vec3d t_ij = cross(nvec, f_contact_tan);
			sys->contact_torque[par_num[0]] += a0*t_ij;
			sys->contact_torque[par_num[1]] += a1*t_ij;
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
	vec3d dv = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
	if (sys->in_predictor && zshift != 0) {
		dv.x += zshift*sys->vel_difference;
	}
	contact_velocity = dv-cross(a0*sys->ang_velocity[par_num[0]]+a1*sys->ang_velocity[par_num[1]], nvec);
	normal_relative_velocity = dot(dv, nvec); // Negative when approaching.

}

/*********************************
 *                                *
 *  Lubrication Forces Methods    *
 *                                *
 *********************************/

void
Interaction::calcLubConstants(){
	lambda = a1/a0;
	invlambda = 1/lambda;
	lambda_square = lambda*lambda;
	lambda_cubic = lambda*lambda*lambda;
	lambda_p_1 = 1+lambda;
	lambda_p_1_square = lambda_p_1*lambda_p_1;
	lambda_p_1_cubic = lambda_p_1_square*lambda_p_1;
	a0a0_23 = a0*a0*(2./3);
	a1a1_23 = a1*a1*(2./3);
	roro_16 = ro*ro*(1./6);
	double a0a0a0 = a0*a0*a0;
	double a1a1a1 = a1*a1*a1;
	double rororo = ro*ro*ro;
	a0a0a0_43 = a0a0a0*(4./3);
	a1a1a1_43 = a1a1a1*(4./3);
	rororo_16 = rororo*(1./6);
	a0a0a0_109 = a0a0a0*(10./9);
	a1a1a1_109 = a1a1a1*(10./9);
	rororo_536 = rororo*(5./36);
	/* XA
	 * X_{a,b}(l) = X_{b,a}(l) = X_{3-a,3-b}(1/l)
	 * X21(l) = X12(l)
	 * X22(l) = X11(1/l)
	 */
	g1_XA = func_g1_XA(lambda);
	g1_inv_XA = func_g1_XA(invlambda);
	cXA[0] = g1_XA;
	cXA[1] = (-2/lambda_p_1)*g1_XA;
	cXA[2] = cXA[1];
	cXA[3] = g1_inv_XA;
	/* YA
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YA = func_g2_YA(lambda);
	g2_inv_YA = func_g2_YA(invlambda);
	cYA[0] = g2_YA;
	cYA[1] = (-2/lambda_p_1)*g2_YA;
	cYA[2] = cYA[1];
	cYA[3] = g2_inv_YA;
	/* YB
	 * Y_{a,b}(l) = -Y_{3-a,3-b}(1/l)
	 * Y21(l) = -Y12(1/l)
	 * Y22(l) = -Y11(1/l)
	 */
	g2_YB = func_g2_YB(lambda);
	g2_inv_YB = func_g2_YB(invlambda);
	cYB[0] = g2_YB;
	cYB[1] = -4/lambda_p_1_square*g2_YB;
	cYB[2] = 4*lambda_square/lambda_p_1_square*g2_inv_YB;
	cYB[3] = -g2_inv_YB;
	/* YC
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l})
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YC = func_g2_YC(lambda);
	g2_inv_YC = func_g2_YC(invlambda);
	g4_YC = func_g4_YC(lambda);
	cYC[0] = g2_YC;
	cYC[1] = g4_YC;
	cYC[2] = cYC[1];
	cYC[3] = g2_inv_YC;
	/* XG
	 * X_{a,b}(l) = -X_{3-a,3-b}(1/l)
	 * X21(l) = -X12(1/l)
	 * X22(l) = -X11(1/l)
	 */
	g1_XG = func_g1_XG(lambda);
	g1_inv_XG = func_g1_XG(invlambda);
	cXG[0] = g1_XG;
	cXG[1] = -4/lambda_p_1_square*g1_XG;
	cXG[2] = 4*lambda_square/lambda_p_1_square*g1_inv_XG;
	cXG[3] = -g1_inv_XG;
	/* YG
	 * Y_{a,b}(l) = -Y_{3-a,3-b}(1/l)
	 * Y21(l) = -Y12(1/l)
	 * Y22(l) = -Y11(1/l)
	 */
	g2_YG = func_g2_YG(lambda);
	g2_inv_YG = func_g2_YG(invlambda);
	cYG[0] = g2_YG;
	cYG[1] = -(4/lambda_p_1_square)*g2_YG;
	cYG[2] = (4*lambda_square/lambda_p_1_square)*g2_inv_YG;
	cYG[3] = -g2_inv_YG;
	/* YH
	 * Y_{a,b}(l) = Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(1/l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YH = func_g2_YH(lambda);
	g2_inv_YH = func_g2_YH(invlambda);
	g5_YH = func_g5_YH(lambda);
	g5_inv_YH = func_g5_YH(invlambda);
	cYH[0] = g2_YH;
	cYH[1] = (8/lambda_p_1_cubic)*g5_YH;
	cYH[2] = (8*lambda_cubic/lambda_p_1_cubic)*g5_inv_YH;
	cYH[3] = g2_inv_YH;
	/* XM
	 * X_{a,b}(l) = X_{b,a}(l)= X_{3-a,3-b}(1/l)
	 * X21(l) = X12(l)
	 * X22(l) = X11(1/l)
	 */
	g1_XM = func_g1_XM(lambda);
	g1_inv_XM = func_g1_XM(invlambda);
	g4_XM = func_g4_XM(lambda);
	cXM[0] = g1_XM;
	cXM[1] = (8/lambda_p_1_cubic)*g4_XM;
	cXM[2] = cXM[1];
	cXM[3] = g1_inv_XM;
	/* YM
	 * Y_{a,b}(l) = Y_{b,a}(l)= Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YM = func_g2_YM(lambda);
	g2_inv_YM = func_g2_YM(invlambda);
	g5_YM = func_g5_YM(lambda);
	cYM[0] = g2_YM;
	cYM[1] = (8/lambda_p_1_cubic)*g5_YM;
	cYM[2] = cYM[1];
	cYM[3] = g2_inv_YM;
}

// Resistance functions
void
Interaction::calcXFunctions(){
	for (int j=0; j<4; j++) {
		XA[j] = cXA[j]*lub_coeff;
		XG[j] = cXG[j]*lub_coeff;
	}
}

void
Interaction::calcXFunctionsStress(){
	for (int j=0; j<4; j++) {
		XA[j] = cXA[j]*lub_coeff;
		XG[j] = cXG[j]*lub_coeff;
		XM[j] = cXM[j]*lub_coeff;
	}
}

void
Interaction::calcXYFunctions(){
	for (int j=0; j<4; j++) {
		XA[j] = cXA[j]*lub_coeff;
		YA[j] = cYA[j]*log_lub_coeff;
		YB[j] = cYB[j]*log_lub_coeff;
		YC[j] = cYC[j]*log_lub_coeff;
		XG[j] = cXG[j]*lub_coeff;
		YG[j] = cYG[j]*log_lub_coeff;
		YH[j] = cYH[j]*log_lub_coeff;
	}
}

void
Interaction::calcXYFunctionsStress(){
	for (int j=0; j<4; j++) {
		XA[j] = cXA[j]*lub_coeff;
		YA[j] = cYA[j]*log_lub_coeff;
		YB[j] = cYB[j]*log_lub_coeff;
		YC[j] = cYC[j]*log_lub_coeff;
		XG[j] = cXG[j]*lub_coeff;
		YG[j] = cYG[j]*log_lub_coeff;
		YH[j] = cYH[j]*log_lub_coeff;
		XM[j] = cXM[j]*lub_coeff;
		YM[j] = cYM[j]*log_lub_coeff;
	}
}

void
Interaction::calcGE(double *GEi, double *GEj){
	/* NOTE:
	 * Calculation of XG and YG needs to be done before that.
	 *
	 * lubrication_model = 1
	 * 1/xi level
	 *
	 * GE1 = nx*nz*(XG11+XG21)*nvec
	 * GE2 = nx*nz*(XG12+XG22)*nvec
	 */
	double cGE_i = (scaledXG0()+scaledXG2())*nxnz;
	double cGE_j = (scaledXG1()+scaledXG3())*nxnz;
	GEi[0] = cGE_i*nvec.x;
	GEi[1] = cGE_i*nvec.y;
	GEi[2] = cGE_i*nvec.z;
	GEj[0] = cGE_j*nvec.x;
	GEj[1] = cGE_j*nvec.y;
	GEj[2] = cGE_j*nvec.z;
}

void
Interaction::calcGEHE(double *GEi, double *GEj, double *HEi, double *HEj){
	/*
	 * lubrication_model = 1
	 * upto log(1/xi) level
	 *
	 * GE1 = nx*nz*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(nz,0,nx);
	 * GE2 = nx*nz*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(nz,0,nx);
	 */
	double nxnx_nznz = nxnx-nznz;
	double YG0_YG2 = scaledYG0()+scaledYG2();
	double YG1_YG3 = scaledYG1()+scaledYG3();
	double cGE_i = (scaledXG0()+scaledXG2()-2*YG0_YG2)*nxnz;
	double cGE_j = (scaledXG1()+scaledXG3()-2*YG1_YG3)*nxnz;
	double cHE_i = scaledYH0()+scaledYH2();
	double cHE_j = scaledYH3()+scaledYH1();
	GEi[0] = cGE_i*nvec.x+YG0_YG2*nvec.z;
	GEi[1] = cGE_i*nvec.y;
	GEi[2] = cGE_i*nvec.z+YG0_YG2*nvec.x;
	GEj[0] = cGE_j*nvec.x+YG1_YG3*nvec.z;
	GEj[1] = cGE_j*nvec.y;
	GEj[2] = cGE_j*nvec.z+YG1_YG3*nvec.x;
	HEi[0] =  cHE_i*nxny;
	HEi[1] = -cHE_i*nxnx_nznz;
	HEi[2] = -cHE_i*nynz;
	HEj[0] =  cHE_j*nxny;
	HEj[1] = -cHE_j*nxnx_nznz;
	HEj[2] = -cHE_j*nynz;
}

// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jeffrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
void
Interaction::pairVelocityStresslet(const vec3d &vi, const vec3d &vj,
								   const vec3d &oi, const vec3d &oj,
								   StressTensor &stresslet_i, StressTensor &stresslet_j){
	/*
	 * (xx, xy, xz, yz, yy, zz)
	 *
	 * S_{11}^{GUX}+S_{12}^{GUX}
	 *         = - vec{n}.(XG11*v1+XG12*v2)*(ninj-(1/3)*delta_{ij})
	 *
	 */
	double cXG_i = -dot(nvec, scaledXG0()*vi+scaledXG1()*vj);
	double cXG_j = -dot(nvec, scaledXG2()*vi+scaledXG3()*vj);
	StressTensor XGU_i(nxnx-c13, nxny, nxnz, nynz, nyny-c13, nznz-c13);
	StressTensor XGU_j = XGU_i;
	XGU_i *= cXG_i;
	XGU_j *= cXG_j;
	stresslet_i = XGU_i;
	stresslet_j = XGU_j;
	if (sys->lubrication_model == 1){
		return;
	}
	StressTensor YGU_i;
	StressTensor YGU_j;
	double cYGUi_xx = nvec.x*vi.x+nvec.x*vi.x-2*nvec.x*nvec.x*dot(nvec, vi);
	double cYGUj_xx = nvec.x*vj.x+nvec.x*vj.x-2*nvec.x*nvec.x*dot(nvec, vj);
	YGU_i.elm[0]  = -scaledYG0()*cYGUi_xx-scaledYG1()*cYGUj_xx;
	YGU_j.elm[0]  = -scaledYG2()*cYGUi_xx-scaledYG3()*cYGUj_xx;
	
	double cYGUi_xy = nvec.x*vi.y+nvec.y*vi.x-2*nvec.x*nvec.y*dot(nvec, vi);
	double cYGUj_xy = nvec.x*vj.y+nvec.y*vj.x-2*nvec.x*nvec.y*dot(nvec, vj);
	YGU_i.elm[1]  = -scaledYG0()*cYGUi_xy-scaledYG1()*cYGUj_xy;
	YGU_j.elm[1]  = -scaledYG2()*cYGUi_xy-scaledYG3()*cYGUj_xy;
	
	double cYGUi_xz = nvec.x*vi.z+nvec.z*vi.x-2*nvec.x*nvec.z*dot(nvec, vi);
	double cYGUj_xz = nvec.x*vj.z+nvec.z*vj.x-2*nvec.x*nvec.z*dot(nvec, vj);
	YGU_i.elm[2]  = -scaledYG0()*cYGUi_xz-scaledYG1()*cYGUj_xz;
	YGU_j.elm[2]  = -scaledYG2()*cYGUi_xz-scaledYG3()*cYGUj_xz;
	
	double cYGUi_yz = nvec.y*vi.z+nvec.z*vi.y-2*nvec.y*nvec.z*dot(nvec, vi);
	double cYGUj_yz = nvec.y*vj.z+nvec.z*vj.y-2*nvec.y*nvec.z*dot(nvec, vj);
	YGU_i.elm[3]  = -scaledYG0()*cYGUi_yz-scaledYG1()*cYGUj_yz;
	YGU_j.elm[3]  = -scaledYG2()*cYGUi_yz-scaledYG3()*cYGUj_yz;

	double cYGUi_yy = nvec.y*vi.y+nvec.y*vi.y-2*nvec.y*nvec.y*dot(nvec, vi);
	double cYGUj_yy = nvec.y*vj.y+nvec.y*vj.y-2*nvec.y*nvec.y*dot(nvec, vj);
	
	YGU_i.elm[4]  = -scaledYG0()*cYGUi_yy-scaledYG1()*cYGUj_yy;
	YGU_j.elm[4]  = -scaledYG2()*cYGUi_yy-scaledYG3()*cYGUj_yy;

	double cYGUi_zz = nvec.z*vi.z+nvec.z*vi.z-2*nvec.z*nvec.z*dot(nvec, vi);
	double cYGUj_zz = nvec.z*vj.z+nvec.z*vj.z-2*nvec.z*nvec.z*dot(nvec, vj);
	
	YGU_i.elm[5]  = -scaledYG0()*cYGUi_zz-scaledYG1()*cYGUj_zz;
	YGU_j.elm[5]  = -scaledYG2()*cYGUi_zz-scaledYG3()*cYGUj_zz;
	
	stresslet_i += YGU_i;
	stresslet_j += YGU_j;
	
	StressTensor YHO_i;
	StressTensor YHO_j;
	double cYHOi_xx = nxnz*oi.y-nxny*oi.z;
	double cYHOj_xx = nxnz*oj.y-nxny*oj.z;
	YHO_i.elm[0]  = -2*(scaledYM0()*cYHOi_xx+scaledYM1()*cYHOj_xx);
	YHO_j.elm[0]  = -2*(scaledYM2()*cYHOi_xx+scaledYM3()*cYHOj_xx);

	double cYHOi_xy = nxnx*oi.z-nxnz*oi.x+nynz*oi.y-nyny*oi.z;
	double cYHOj_xy = nxnx*oj.z-nxnz*oj.x+nynz*oj.y-nyny*oj.z;
	YHO_i.elm[1]  = -scaledYM0()*cYHOi_xy-scaledYM1()*cYHOj_xy;
	YHO_j.elm[1]  = -scaledYM2()*cYHOi_xy-scaledYM3()*cYHOj_xy;
	
	double cYHOi_xz = nxny*oi.x-nxnx*oi.y+nznz*oi.y-nynz*oi.z;
	double cYHOj_xz = nxny*oj.x-nxnx*oj.y+nznz*oj.y-nynz*oj.z;
	YHO_i.elm[2]  = -scaledYM0()*cYHOi_xz-scaledYM1()*cYHOj_xz;
	YHO_j.elm[2]  = -scaledYM2()*cYHOi_xz-scaledYM3()*cYHOj_xz;
	
	double cYHOi_yz = nyny*oi.x-nynz*oi.y+nxnz*oi.z-nynz*oi.x;
	double cYHOj_yz = nyny*oj.x-nynz*oj.y+nxnz*oj.z-nynz*oj.x;
	YHO_i.elm[3]  = -scaledYM0()*cYHOi_yz-scaledYM1()*cYHOj_yz;
	YHO_j.elm[3]  = -scaledYM2()*cYHOi_yz-scaledYM3()*cYHOj_yz;
	
	double cYHOi_yy = nxny*oi.z-nynz*oi.x;
	double cYHOj_yy = nxny*oj.z-nynz*oj.x;
	YHO_i.elm[4]  = -2*(scaledYM0()*cYHOi_yy+scaledYM1()*cYHOj_yy);
	YHO_j.elm[4]  = -2*(scaledYM2()*cYHOi_yy+scaledYM3()*cYHOj_yy);
	
	double cYHOi_zz = nynz*oi.x-nxnz*oi.y;
	double cYHOj_zz = nynz*oj.x-nxnz*oj.y;
	YHO_i.elm[5]  = -2*(scaledYM0()*cYHOi_zz+scaledYM1()*cYHOj_zz);
	YHO_j.elm[5]  = -2*(scaledYM2()*cYHOi_zz+scaledYM3()*cYHOj_zz);
	
	stresslet_i += YHO_i;
	stresslet_j += YHO_j;
}

// convenient interface for pairVelocityStresslet(const vec3d &vi, const vec3d &vj, stresslet &stresslet_i, stresslet &stresslet_j)
//void
//Interaction::pairVelocityStresslet(double* &vel_array, StressTensor &stresslet_i, StressTensor &stresslet_j){
//	vec3d vi, vj;
//	int i3 = 3*par_num[0];
//	int j3 = 3*par_num[1];
//	vi.x = vel_array[i3];
//	vi.y = vel_array[i3+1];
//	vi.z = vel_array[i3+2];
//	vj.x = vel_array[j3];
//	vj.y = vel_array[j3+1];
//	vj.z = vel_array[j3+2];
//	pairVelocityStresslet(vi, vj, stresslet_i, stresslet_j);
//}

void
Interaction::pairStrainStresslet(StressTensor &stresslet_i, StressTensor &stresslet_j){
	/*
	 * S_{ab} = M_{ijxz}E_{xz}+M_{ijzx}E_{zx}
	 *        = (rate/2)*(M_{ijxz}+M_{ijzx})
	 *        = M_{ijxz}
	 * rate = 1 and M_{ijxz} = M_{ijzx}
	 *   M_{ijxz} = XM_{ab}(3/2)(ninj-(1/3)delta_ij)*(nznx)
	 *   M_{ijzx} = XM_{ab}(3/2)(ninj-(1/3)delta_ij)*(nznx)
	 *
	 * (S_{11}+S_{12})^{ME,X}
	 *   = M11_{ijxz} + M12_{ijxz}
     *   = XM_{11}(3/2)(ninj-(1/3)delta_ij)*(nznx) + XM_{12}(3/2)(ninj-(1/3)delta_ij)*(nznx)
	 *   = [(3/2)(XM_{11}+XM_{12})*nxnz]*(ninj-(1/3)delta_ij)
	 */
	double cXM_i = (3./2)*(scaledXM0()+scaledXM1())*nxnz;
	double cXM_j = (3./2)*(scaledXM2()+scaledXM3())*nxnz;
	StressTensor XME_i(nxnx-c13, nxny, nxnz, nynz, nyny-c13, nznz-c13);
	StressTensor XME_j = XME_i;
	XME_i *= cXM_i;
	XME_j *= cXM_j;
	stresslet_i = XME_i;
	stresslet_j = XME_j;
	if (sys->lubrication_model == 1){
		return;
	}
	double cYM_i = (1./2)*(scaledYM0()+scaledYM1());
	double cYM_j = (1./2)*(scaledYM2()+scaledYM3());
	StressTensor YME_i(2*nxnz     -4*nxnx*nxnz,
					     nynz     -4*nxny*nxnz,
						 nxnx+nznz-4*nxnz*nxnz,
					     nxny     -4*nynz*nxnz,
								  -4*nyny*nxnz,
					   2*nxnz     -4*nznz*nxnz);
	StressTensor YME_j = YME_i;
	YME_i *= cYM_i;
	YME_j *= cYM_j;
	stresslet_i += YME_i;
	stresslet_j += YME_j;
}


void
Interaction::calcTestStress(){
	// @@@@@ FOR TEST
	StressTensor stresslet_GU_HO_i;
	StressTensor stresslet_GU_HO_j;
	StressTensor stresslet_ME_i;
	StressTensor stresslet_ME_j;
	vec3d vi(sys->v_total[i6], sys->v_total[i6+1], sys->v_total[i6+2]);
	vec3d vj(sys->v_total[j6], sys->v_total[j6+1], sys->v_total[j6+2]);
	vec3d oi(sys->v_total[i6+3], sys->v_total[i6+4], sys->v_total[i6+5]);
	vec3d oj(sys->v_total[j6+3], sys->v_total[j6+4], sys->v_total[j6+5]);
	/*
	 *  First: -G*(U-Uinf) term
	 */
	pairVelocityStresslet(vi, vj, oi, oj, stresslet_GU_HO_i, stresslet_GU_HO_j);
	/*
	 *  Second: +M*Einf term
	 */
	pairStrainStresslet(stresslet_ME_i, stresslet_ME_j);
	sys->test_totalstress[par_num[0]] += stresslet_GU_HO_i+stresslet_ME_i;
	sys->test_totalstress[par_num[1]] += stresslet_GU_HO_j+stresslet_ME_j;
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
	//	calcXA();
	if (sys->lubrication_model == 1){
		calcXFunctions();
	} else if (sys->lubrication_model == 2){
		calcXYFunctions();
	}
	double cf_AU_i = -dot(a0*XA[0]*vi+ro_12*XA[1]*vj, nvec);
	/*
	 *  Second -tildeG*(-Einf)term
	 */
	//calcXG();
	double cf_GE_i = nxnz*(a0a0_23*XG[0]+roro_16*XG[2]);
	lubforce_i = (cf_AU_i+cf_GE_i)*nvec;
}

void
Interaction::addHydroStress(){
	StressTensor stresslet_GU_i;
	StressTensor stresslet_GU_j;
	StressTensor stresslet_ME_i;
	StressTensor stresslet_ME_j;
	vec3d vi(sys->v_hydro[i6], sys->v_hydro[i6+1], sys->v_hydro[i6+2]);
	vec3d vj(sys->v_hydro[j6], sys->v_hydro[j6+1], sys->v_hydro[j6+2]);
	vec3d oi(sys->v_hydro[i6+3], sys->v_hydro[i6+4], sys->v_hydro[i6+5]);
	vec3d oj(sys->v_hydro[j6+3], sys->v_hydro[j6+4], sys->v_hydro[j6+5]);
	/*
	 *  First: -G*(U-Uinf) term
	 */
	pairVelocityStresslet(vi, vj, oi, oj, stresslet_GU_i, stresslet_GU_j);
	/*
	 *  Second: +M*Einf term
	 */
	pairStrainStresslet(stresslet_ME_i, stresslet_ME_j);
	sys->lubstress[par_num[0]] += stresslet_GU_i+stresslet_ME_i;
	sys->lubstress[par_num[1]] += stresslet_GU_j+stresslet_ME_j;
}

void
Interaction::addContactStress(){
	if (contact) {
		/*
		 * Fc_normal_norm = -kn_scaled*gap_nondim; --> positive
		 * Fc_normal = -Fc_normal_norm*nvec;
		 * This force acts on particle 1.
		 * stress1 is a0*nvec[*]force.
		 * stress2 is (-a1*nvec)[*](-force) = a1*nvec[*]force
		 */
		contact_stresslet_XF_normal.set(r_vec, f_contact_normal);
		contact_stresslet_XF_tan.set(r_vec, f_contact_tan);
	} else {
		contact_stresslet_XF_normal.reset();
		contact_stresslet_XF_tan.reset();
	}

	// Add term G*V_cont
	StressTensor stresslet_GU_i;
	StressTensor stresslet_GU_j;
	vec3d vi(sys->v_cont[i6], sys->v_cont[i6+1], sys->v_cont[i6+2]);
	vec3d vj(sys->v_cont[j6], sys->v_cont[j6+1], sys->v_cont[j6+2]);
	vec3d oi(sys->v_cont[i6+3], sys->v_cont[i6+4], sys->v_cont[i6+5]);
	vec3d oj(sys->v_cont[j6+3], sys->v_cont[j6+4], sys->v_cont[j6+5]);
	pairVelocityStresslet(vi, vj, oi, oj, stresslet_GU_i, stresslet_GU_j);
	sys->contactstressGU[par_num[0]] += stresslet_GU_i;
	sys->contactstressGU[par_num[1]] += stresslet_GU_j;

}

void
Interaction::addColloidalStress(){
	colloidal_stresslet_XF.set(r_vec, f_colloidal);
	// Add term G*V_cont
	StressTensor stresslet_colloid_GU_i;
	StressTensor stresslet_colloid_GU_j;
	vec3d vi(sys->v_colloidal[i6], sys->v_colloidal[i6+1], sys->v_colloidal[i6+2]);
	vec3d vj(sys->v_colloidal[j6], sys->v_colloidal[j6+1], sys->v_colloidal[j6+2]);
	vec3d oi(sys->v_colloidal[i6+3], sys->v_colloidal[i6+4], sys->v_colloidal[i6+5]);
	vec3d oj(sys->v_colloidal[j6+3], sys->v_colloidal[j6+4], sys->v_colloidal[j6+5]);

	pairVelocityStresslet(vi, vj, oi, oj, stresslet_colloid_GU_i, stresslet_colloid_GU_j);
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
	duration = sys->get_shear_strain()-strain_lub_start;
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
	return contact_velocity.norm();
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
