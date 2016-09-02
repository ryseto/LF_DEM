//
//  Lubrication.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <tuple>
#include "ContactDashpot.h"
#include "LubricationFunctions.h"
#include "Interaction.h"
#include "System.h"

ContactDashpot::ContactDashpot():
_active(false),
p0(0),
p1(0),
p0_6(0),
p1_6(0),
a0(0),
a1(0),
ro(0),
ro_12(0),
normal_coeff(0),
tangential_coeff(0)
{
	for (int i=0; i<4; i++) {
		XA[i] = 0;
		YA[i] = 0;
		YB[i] = 0;
		YC[i] = 0;
	}
}

void ContactDashpot::init(System *sys_, Interaction* int_)
{
	sys = sys_;
	interaction = int_;
	nvec = &(interaction->nvec);
}

void ContactDashpot::setParticleData()
{
	std::tie(p0, p1) = interaction->get_par_num();
	p0_6 = 6*p0;
	p1_6 = 6*p1;
	a0 = sys->radius[p0];
	a1 = sys->radius[p1];
	ro = a0+a1;
	ro_12 = ro/2;
}

void ContactDashpot::activate()
{
	_active = true;
	sys->updateNumberOfPairwiseResistances(p0, p1, +1);
}

void ContactDashpot::deactivate()
{
	_active = false;
	sys->updateNumberOfPairwiseResistances(p0, p1, -1);
}

void ContactDashpot::setDashpotResistanceCoeffs(double kn, double kt,
												double rtime_normal, double rtime_tan)
{
	if (rtime_normal > 0) {
		/* t = beta/kn
		 *  beta = t*kn
		 * normal_coeff = 4*beta = 4*kn*rtime_normal
		 */
		normal_coeff = 4*kn*rtime_normal;
	} else {
		if (sys->lubrication) {// take the same resistance as lubrication
			// 1/(h+c) --> 1/c
			normal_coeff = 1/sys->p.lub_reduce_parameter;
		} else {
			throw std::runtime_error(" ContactDashpot:: Error: normal relaxation time set negative, but no lubrication.");
		}
	}

	if (rtime_tan > 0) {
		if (sys->lubrication) { // the contact can get unstable if the tangential resistance difference is too big between with and wihout contact
			throw std::runtime_error(" ContactDashpot:: Error: with lubrication, tangential relaxation time cannot be set positive.");
		}
		tangential_coeff = 6*kt*rtime_tan;
	} else {
		if (sys->lubrication) {// take the same resistance as lubrication
			// 1/(h+c) --> 1/c
			if (sys->p.lub_reduce_parameter < 1) {
				tangential_coeff = log(1/sys->p.lub_reduce_parameter);
			} else {
				tangential_coeff = 0;
			}
		} else {
			throw std::runtime_error(" ContactDashpot:: Error: tangential relaxation time set negative, but no lubrication.");
		}
	}
	calcDashpotResistances();
}

void ContactDashpot::calcDashpotResistances()
{
	double lambda = a1/a0;
	double invlambda = 1/lambda;
	double lambda_square = lambda*lambda;
	double lambda_p_1 = 1+lambda;
	double lambda_p_1_square = lambda_p_1*lambda_p_1;
	double a0a0_23 = a0*a0*(2.0/3);
	double a1a1_23 = a1*a1*(2.0/3);
	double roro_16 = ro*ro*(1.0/6);
	double a0a0a0_43 = a0*a0*a0*(4.0/3);
	double a1a1a1_43 = a1*a1*a1*(4.0/3);
	double rororo_16 = ro*ro*ro*(1.0/6);
	/* XA
	 * X_{a,b}(l) = X_{b,a}(l) = X_{3-a,3-b}(1/l)
	 * X21(l) = X12(l)
	 * X22(l) = X11(1/l)
	 */
	double g1_XA     = func_g1_XA(lambda);
	double g1_inv_XA = func_g1_XA(invlambda);
	XA[0] = a0*g1_XA;
	XA[1] = ro_12*(-2/lambda_p_1)*g1_XA;
	XA[2] = ro_12*XA[1];
	XA[3] = a1*g1_inv_XA;
	for (int j=0; j<4; j++) {
		XA[j] *= normal_coeff;
	}

	/* YA
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	double g2_YA     = func_g2_YA(lambda);
	double g2_inv_YA = func_g2_YA(invlambda);
	YA[0] = a0*g2_YA;
	YA[1] = ro_12*(-2/lambda_p_1)*g2_YA;
	YA[2] = ro_12*YA[1];
	YA[3] = a1*g2_inv_YA;
	for (int j=0; j<4; j++) {
		YA[j] *= tangential_coeff;
	}

	/* YB
	 * Y_{a,b}(l) = -Y_{3-a,3-b}(1/l)
	 * Y21(l) = -Y12(1/l)
	 * Y22(l) = -Y11(1/l)
	 */
	double g2_YB     = func_g2_YB(lambda);
	double g2_inv_YB = func_g2_YB(invlambda);
	YB[0] = a0a0_23*g2_YB;
	YB[1] = roro_16*(-4/lambda_p_1_square)*g2_YB;
	YB[2] = roro_16*(4*lambda_square/lambda_p_1_square)*g2_inv_YB;
	YB[3] = -a1a1_23*g2_inv_YB;
	for (int j=0; j<4; j++) {
		YB[j] *= tangential_coeff;
	}

	/* YC
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l})
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	double g2_YC     = func_g2_YC(lambda);
	double g2_inv_YC = func_g2_YC(invlambda);
	double g4_YC     = func_g4_YC(lambda);
	YC[0] = a0a0a0_43*g2_YC;
	YC[1] = rororo_16*g4_YC;
	YC[2] = rororo_16*YC[1];
	YC[3] = a1a1a1_43*g2_inv_YC;
	for (int j=0; j<4; j++) {
		YC[j] *= tangential_coeff;
	}
}

// for FT/UW version
struct ODBlock ContactDashpot::RFU_ODBlock() const
{
	double n0n0 = nvec->x*nvec->x;
	double n0n1 = nvec->x*nvec->y;
	double n0n2 = nvec->x*nvec->z;
	double n1n1 = nvec->y*nvec->y;
	double n1n2 = nvec->y*nvec->z;
	double n2n2 = nvec->z*nvec->z;
	double one_n0n0 = 1-n0n0;
	double one_n1n1 = 1-n1n1;
	double one_n2n2 = 1-n2n2;
	double XA1mYA1 = XA[1] - YA[1];
	struct ODBlock block;
	// column 0
	block.col0[0] =  XA[1]*n0n0 + YA[1]*one_n0n0;
	block.col0[1] = XA1mYA1*n0n1;
	block.col0[2] = XA1mYA1*n0n2;
	block.col0[3] = -YB[2]*nvec->z;
	block.col0[4] =  YB[2]*nvec->y;
	// column 1
	block.col1[0] =  XA[1]*n1n1 + YA[1]*one_n1n1;
	block.col1[1] = XA1mYA1*n1n2;
	block.col1[2] = -YB[2]*nvec->x;
	// column 2
	block.col2[0] =  XA[1]*n2n2 + YA[1]*one_n2n2;
	// column 3
	block.col3[0] =  YB[1]*nvec->z;
	block.col3[1] = -YB[1]*nvec->y;
	block.col3[2] =  YC[1]*one_n0n0;
	block.col3[3] = -YC[1]*n0n1;
	block.col3[4] = -YC[1]*n0n2;
	// column 4
	block.col4[0] =  YB[1]*nvec->x;
	block.col4[1] =  YC[1]*one_n1n1;
	block.col4[2] = -YC[1]*n1n2;
	// column 5
	block.col5[0] =  YC[1]*one_n2n2;
	return block;
}

// Diagonal Blocks Terms, FT/UW version
std::pair<struct DBlock, struct DBlock> ContactDashpot::RFU_DBlocks() const
{
	struct DBlock b0;
	struct DBlock b1;

	double n0n0 = nvec->x*nvec->x;
	double n0n1 = nvec->x*nvec->y;
	double n0n2 = nvec->x*nvec->z;
	double n1n1 = nvec->y*nvec->y;
	double n1n2 = nvec->y*nvec->z;
	double n2n2 = nvec->z*nvec->z;
	double one_n0n0 = 1-n0n0;
	double one_n1n1 = 1-n1n1;
	double one_n2n2 = 1-n2n2;

	double XA0mYA0 = XA[0] - YA[0];

	// (*,0)
	b0.col0[0] =  XA[0]*n0n0 + YA[0]*one_n0n0; // 00 element of the dblock
	b0.col0[1] = XA0mYA0*n0n1;           // 10
	b0.col0[2] = XA0mYA0*n0n2;           // 20
	b0.col0[3] = -YB[0]*nvec->z;                   // 40
	b0.col0[4] =  YB[0]*nvec->y;                   // 50
	// (*,1)
	b0.col1[0] =  XA[0]*n1n1 + YA[0]*one_n1n1; // 11
	b0.col1[1] = XA0mYA0*n1n2;           // 21
	b0.col1[2] = -YB[0]*nvec->x;                   // 51
	// (*,2)
	b0.col2[0] =  XA[0]*n2n2 + YA[0]*one_n2n2; // 22
	// (*,3)
	b0.col3[0] =  YC[0]*one_n0n0;                 // 33
	b0.col3[1] = -YC[0]*n0n1;                     // 43
	b0.col3[2] = -YC[0]*n0n2;                     // 53
	// (*,4)
	b0.col4[0] =  YC[0]*one_n1n1;                 // 44
	b0.col4[1] = -YC[0]*n1n2;                     // 54
	// (*,5)
	b0.col5[0] =  YC[0]*one_n2n2;                 // 55

	double XA3mYA3 = XA[3] - YA[3];
// (*,0)
	b1.col0[0] =  XA[3]*n0n0 + YA[3]*one_n0n0; // 00 element of the dblock
	b1.col0[1] = XA3mYA3*n0n1;           // 10
	b1.col0[2] = XA3mYA3*n0n2;           // 20
	b1.col0[3] = -YB[3]*nvec->z;                   // 40
	b1.col0[4] =  YB[3]*nvec->y;                   // 50
	// (*,1)
	b1.col1[0] =  XA[3]*n1n1 + YA[3]*one_n1n1; // 11
	b1.col1[1] = XA3mYA3*n1n2;           // 21
	b1.col1[2] = -YB[3]*nvec->x;                   // 51
	// (*,2)
	b1.col2[0] =  XA[3]*n2n2 + YA[3]*one_n2n2; // 22
	// (*,3)
	b1.col3[0] =  YC[3]*one_n0n0;                 // 33
	b1.col3[1] = -YC[3]*n0n1;                     // 43
	b1.col3[2] = -YC[3]*n0n2;                     // 53
	// (*,4)
	b1.col4[0] =  YC[3]*one_n1n1;                 // 44
	b1.col4[1] = -YC[3]*n1n2;                     // 54
	// (*,5)
	b1.col5[0] =  YC[3]*one_n2n2;                 // 55

	return std::make_pair(b0, b1);
}

vec3d ContactDashpot::getForceOnP0(const vec3d &vel_p0,
                                   const vec3d &vel_p1,
                                   const vec3d &ang_vel_p0,
                                   const vec3d &ang_vel_p1) const
{
	/** \brief Resistance force acting on particle p0.

	 Used by the dynamics to get the R_FU*U_inf force term and by the output data.

	 vel_p0 (ang_vel_p0) is the TOTAL (angular) velocity of p0 (not the non-affine part),
	 vel_p1 (ang_vel_p1) is the total (angular) velocity of p1.
	 This method deals with Lees-Edwards PBC by itself, so vel_p0 and vel_p1
	 are the actual velocities in System.
	 */

	/*
	 *  First: -A*(U-Uinf) term
	 * Eq. (1.6a) in Jeffrey&Onishi 1984
	 * A_{ij}^{ab} = XA_{ab}ni*nj + YA_{ab}(del_{ij}-ni*nj)
	 * B~_{ji}^{ab} = YB_{ba}epsilon_{jik} nk
	 *
	 */
	if (is_active()) {
		vec3d vi(vel_p0);
		vec3d vj(vel_p1);
		vj += interaction->z_offset*sys->get_vel_difference();
		/* XAU_i */
		vec3d force_p0 = -dot(XA[0]*vi+XA[1]*vj, nvec)*(*nvec);

		/* YAU_i */
		force_p0 += -YA[0]*(vi-(*nvec)*dot(nvec, vi)) - YA[1]*(vj-(*nvec)*dot(nvec, vj));
		/* YBO_i */
		force_p0 += -YB[0]*cross(nvec, ang_vel_p0)    - YB[2]*cross(nvec, ang_vel_p1);
		return force_p0;
	} else {
		return vec3d();
	}
}


vec3d ContactDashpot::getForceOnP0_nonaffine(const vec3d &na_vel_p0,
					                                   const vec3d &na_vel_p1,
					                                   const vec3d &na_ang_vel_p0,
					                                   const vec3d &na_ang_vel_p1) const
{
	/** \brief Resistance force acting on particle p0.

	 Used by the dynamics to get the R_FU*U_inf force term and by the output data.

	 vel_p0 (ang_vel_p0) is the TOTAL (angular) velocity of p0 (not the non-affine part),
	 vel_p1 (ang_vel_p1) is the total (angular) velocity of p1.
	 This method deals with Lees-Edwards PBC by itself, so vel_p0 and vel_p1
	 are the actual velocities in System.
	 */

	/*
	 *  First: -A*(U-Uinf) term
	 * Eq. (1.6a) in Jeffrey&Onishi 1984
	 * A_{ij}^{ab} = XA_{ab}ni*nj + YA_{ab}(del_{ij}-ni*nj)
	 * B~_{ji}^{ab} = YB_{ba}epsilon_{jik} nk
	 *
	 */
	if (is_active()) {
		/* XAU_i */
		vec3d force_p0 = -dot(XA[0]*na_vel_p0+XA[1]*na_vel_p1, nvec)*(*nvec);

		/* YAU_i */
		force_p0 += - YA[0]*(na_vel_p0-(*nvec)*dot(nvec, na_vel_p0))
		            - YA[1]*(na_vel_p1-(*nvec)*dot(nvec, na_vel_p1));
		/* YBO_i */
		force_p0 += - YB[0]*cross(nvec, na_ang_vel_p0)
		            - YB[2]*cross(nvec, na_ang_vel_p1);
		return force_p0;
	} else {
		return vec3d();
	}
}


std::tuple<vec3d, vec3d, vec3d, vec3d> ContactDashpot::getRFU_Uinf(const vec3d &u_inf_p0,
																   const vec3d &u_inf_p1,
																   const vec3d &omega_inf) const
{
	/** \brief */
	if (is_active()) {
		vec3d force_p0 = getForceOnP0(u_inf_p0, u_inf_p1, omega_inf, omega_inf);
		vec3d torque_p0 = cross((*nvec)*a0, force_p0);
		return std::make_tuple(force_p0, -force_p0, torque_p0, (a1/a0)*torque_p0);
	} else {
		return std::make_tuple(vec3d(), vec3d(), vec3d(), vec3d());
	}
}
