//
//  Lubrication.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2016 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <tuple>
#include <stdexcept>
#include "Lubrication.h"
#include "LubricationFunctions.h"
#include "Interaction.h"
#include "System.h"

Lubrication::Lubrication():
ro_12(0), // = ro/2
a0a0_23(0),
a1a1_23(0),
roro_16(0),
a0a0a0_43(0),
a1a1a1_43(0),
rororo_16(0),
a0a0a0_109(0),
a1a1a1_109(0),
rororo_536(0),
g1_XA(0),
g1_inv_XA(0),
g2_YA(0),
g2_inv_YA(0),
g2_YB(0),
g2_inv_YB(0),
g2_YC(0),
g2_inv_YC(0),
g4_YC(0),
g1_XG(0),
g1_inv_XG(0),
g2_YG(0),
g2_inv_YG(0),
g2_YH(0),
g2_inv_YH(0),
g5_YH(0),
g5_inv_YH(0),
g1_XM(0),
g1_inv_XM(0),
g4_XM(0),
g2_YM(0),
g2_inv_YM(0),
g5_YM(0),
lub_coeff(0),
log_lub_coeff(0),
a0(0),
a1(0),
ro(0),
lambda(0),
invlambda(0),
lambda_square(0),
lambda_cubic(0),
lambda_p_1(0),
lambda_p_1_square(0),
lambda_p_1_cubic(0),
_active(false),
range(0),
p0(0),
p1(0),
p0_6(0),
p1_6(0)
{
	for (int i=0; i<4; i++) {
		XA[i] = 0;
		YA[i] = 0;
		YB[i] = 0;
		YC[i] = 0;
		XG[i] = 0;
		YG[i] = 0;
		YH[i] = 0;
		XM[i] = 0;
		YM[i] = 0;
		cXA[i] = 0;
		cYA[i] = 0;
		cYB[i] = 0;
		cYC[i] = 0;
		cXG[i] = 0;
		cYG[i] = 0;
		cYH[i] = 0;
		cXM[i] = 0;
		cYM[i] = 0;
	}
}

void Lubrication::init(System *sys_, Interaction* int_)
{
	interaction = int_;
	sys = sys_;
	nvec = &(interaction->nvec);
}

void Lubrication::setParticleData()
{
	/**
	 	Set things having to do with particle index and radii.
		No position, velocity or force data is used here.
		**/
	p0 = interaction->p0;
	p1 = interaction->p1;
	p0_6 = 6*p0;
	p1_6 = 6*p1;
	calcLubConstants();
	range = sys->calcLubricationRange(p0, p1);
}

void Lubrication::activate()
{
	_active = true;
	sys->updateNumberOfPairwiseResistances(p0, p1, +1);
}

void Lubrication::deactivate()
{
	_active = false;
	sys->updateNumberOfPairwiseResistances(p0, p1, -1);
}

void Lubrication::updateActivationState()
{
	if (sys->p.lubrication_model > 0) {
		bool in_range = interaction->r < range && interaction->r > ro;
		if (!is_active() && in_range) {
			activate();
		}
		if (is_active() && !in_range) {
			deactivate();
		}
	}
}

void Lubrication::setResistanceCoeff(double lub_coeff_, double log_lub_coeff_)
{
	lub_coeff = lub_coeff_; // normal
	log_lub_coeff = log_lub_coeff_; // tangential
	if (sys->p.lubrication_model == 1){
		calcXFunctions();
	}
	if (sys->p.lubrication_model == 2){
		calcXYFunctions();
	}
}

/*********************************
 *                                *
 *  Lubrication Forces Methods    *
 *                                *
 *********************************/

void Lubrication::calcLubConstants()
{
	a0 = interaction->a0;
	a1 = interaction->a1;
	ro = interaction->ro;
	ro_12 = interaction->ro_12;
	lambda = a1/a0;
	invlambda = 1/lambda;
	lambda_square = lambda*lambda;
	lambda_cubic = lambda*lambda*lambda;
	lambda_p_1 = 1+lambda;
	lambda_p_1_square = lambda_p_1*lambda_p_1;
	lambda_p_1_cubic = lambda_p_1_square*lambda_p_1;
	a0a0_23 = a0*a0*(2.0/3);
	a1a1_23 = a1*a1*(2.0/3);
	roro_16 = ro*ro*(1.0/6);
	double a0a0a0 = a0*a0*a0;
	double a1a1a1 = a1*a1*a1;
	double rororo = ro*ro*ro;
	a0a0a0_43 = a0a0a0*(4.0/3);
	a1a1a1_43 = a1a1a1*(4.0/3);
	rororo_16 = rororo*(1.0/6);
	a0a0a0_109 = a0a0a0*(10.0/9);
	a1a1a1_109 = a1a1a1*(10.0/9);
	rororo_536 = rororo*(5.0/36);
	/* XA
	 * X_{a,b}(l) = X_{b,a}(l) = X_{3-a,3-b}(1/l)
	 * X21(l) = X12(l)
	 * X22(l) = X11(1/l)
	 */
	g1_XA     = func_g1_XA(lambda);
	g1_inv_XA = func_g1_XA(invlambda);
	cXA[0] = a0*g1_XA;
	cXA[1] = ro_12*(-2/lambda_p_1)*g1_XA;
	cXA[2] = cXA[1];
	cXA[3] = a1*g1_inv_XA;
	/* YA
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YA     = func_g2_YA(lambda);
	g2_inv_YA = func_g2_YA(invlambda);
	cYA[0] = a0*g2_YA;
	cYA[1] = ro_12*(-2/lambda_p_1)*g2_YA;
	cYA[2] = cYA[1];
	cYA[3] = a1*g2_inv_YA;
	/* YB
	 * Y_{a,b}(l) = -Y_{3-a,3-b}(1/l)
	 * Y21(l) = -Y12(1/l)
	 * Y22(l) = -Y11(1/l)
	 */
	g2_YB     = func_g2_YB(lambda);
	g2_inv_YB = func_g2_YB(invlambda);
	cYB[0] = a0a0_23*g2_YB;
	cYB[1] = roro_16*(-4/lambda_p_1_square)*g2_YB;
	cYB[2] = roro_16*(4*lambda_square/lambda_p_1_square)*g2_inv_YB;
	cYB[3] = -a1a1_23*g2_inv_YB;
	/* YC
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l})
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YC     = func_g2_YC(lambda);
	g2_inv_YC = func_g2_YC(invlambda);
	g4_YC     = func_g4_YC(lambda);
	cYC[0] = a0a0a0_43*g2_YC;
	cYC[1] = rororo_16*g4_YC;
	cYC[2] = cYC[1];
	cYC[3] = a1a1a1_43*g2_inv_YC;
	/* XG
	 * X_{a,b}(l) = -X_{3-a,3-b}(1/l)
	 * X21(l) = -X12(1/l)
	 * X22(l) = -X11(1/l)
	 */
	g1_XG     = func_g1_XG(lambda);
	g1_inv_XG = func_g1_XG(invlambda);
	cXG[0] = a0a0_23*g1_XG;
	cXG[1] = roro_16*(-4/lambda_p_1_square)*g1_XG;
	cXG[2] = roro_16*(4*lambda_square/lambda_p_1_square)*g1_inv_XG;
	cXG[3] = -a1a1_23*g1_inv_XG;
	/* YG
	 * Y_{a,b}(l) = -Y_{3-a,3-b}(1/l)
	 * Y21(l) = -Y12(1/l)
	 * Y22(l) = -Y11(1/l)
	 */
	g2_YG     = func_g2_YG(lambda);
	g2_inv_YG = func_g2_YG(invlambda);
	cYG[0] = a0a0_23*g2_YG;
	cYG[1] = roro_16*(-4/lambda_p_1_square)*g2_YG;
	cYG[2] = roro_16*(4*lambda_square/lambda_p_1_square)*g2_inv_YG;
	cYG[3] = -a1a1_23*g2_inv_YG;
	/* YH
	 * Y_{a,b}(l) = Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(1/l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YH     = func_g2_YH(lambda);
	g2_inv_YH = func_g2_YH(invlambda);
	g5_YH     = func_g5_YH(lambda);
	g5_inv_YH = func_g5_YH(invlambda);
	cYH[0] = a0a0a0_43*g2_YH;
	cYH[1] = rororo_16*(8/lambda_p_1_cubic)*g5_YH;
	cYH[2] = rororo_16*(8*lambda_cubic/lambda_p_1_cubic)*g5_inv_YH;
	cYH[3] = a1a1a1_43*g2_inv_YH;
	/* XM
	 * X_{a,b}(l) = X_{b,a}(l)= X_{3-a,3-b}(1/l)
	 * X21(l) = X12(l)
	 * X22(l) = X11(1/l)
	 */
	g1_XM     = func_g1_XM(lambda);
	g1_inv_XM = func_g1_XM(invlambda);
	g4_XM     = func_g4_XM(lambda);
	cXM[0] = a0a0a0_109*g1_XM;
	cXM[1] = rororo_536*(8/lambda_p_1_cubic)*g4_XM;
	cXM[2] = cXM[1];
	cXM[3] = a1a1a1_109*g1_inv_XM;
	/* YM
	 * Y_{a,b}(l) = Y_{b,a}(l)= Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	g2_YM     = func_g2_YM(lambda);
	g2_inv_YM = func_g2_YM(invlambda);
	g5_YM     = func_g5_YM(lambda);
	cYM[0] = a0a0a0_109*g2_YM;
	cYM[1] = rororo_536*(8/lambda_p_1_cubic)*g5_YM;
	cYM[2] = cYM[1];
	cYM[3] = a1a1a1_109*g2_inv_YM;
}

// Resistance functions
void Lubrication::calcXFunctions()
{
	for (int j=0; j<4; j++) {
		XA[j] = cXA[j]*lub_coeff;
		XG[j] = cXG[j]*lub_coeff;
	}
}

void Lubrication::calcXFunctionsStress()
{
	calcXFunctions();
	for (int j=0; j<4; j++) {
		XM[j] = cXM[j]*lub_coeff;
	}
}

void Lubrication::calcXYFunctions()
{
	calcXFunctions();

	for (int j=0; j<4; j++) {
		YA[j] = cYA[j]*log_lub_coeff;
		YB[j] = cYB[j]*log_lub_coeff;
		YC[j] = cYC[j]*log_lub_coeff;
		YG[j] = cYG[j]*log_lub_coeff;
		YH[j] = cYH[j]*log_lub_coeff;
	}
}

void Lubrication::calcXYFunctionsStress()
{
	calcXYFunctions();
	for (int j=0; j<4; j++) {
		XM[j] = cXM[j]*lub_coeff;
		YM[j] = cYM[j]*log_lub_coeff;
	}
}

std::tuple<vec3d, vec3d> Lubrication::calcGE_squeeze()
{
	/* NOTE:
	 * Calculation of XG and YG needs to be done before that.
	 *
	 * lubrication_model = 1 or 3
	 * 1/xi level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21)*nvec
	 * GE2 = (nvecnvec:E)*(XG12+XG22)*nvec
	 */
	double nnE;
	if (sys->p.cross_shear) {
		double costheta, sintheta;
		std::tie(costheta, sintheta) = sys->getCosSinShearAngle();
		// nnE = nynz;
		nnE = costheta*nvec->x*nvec->z+sintheta*nvec->y*nvec->z;
	} else {
		nnE = nvec->x*nvec->z;
	}
	double cGE_p0 = (XG[0]+XG[2])*nnE;
	double cGE_p1 = (XG[1]+XG[3])*nnE;
	vec3d GEi, GEj;
	GEi.x = cGE_p0*nvec->x;
	GEi.y = cGE_p0*nvec->y;
	GEi.z = cGE_p0*nvec->z;
	GEj.x = cGE_p1*nvec->x;
	GEj.y = cGE_p1*nvec->y;
	GEj.z = cGE_p1*nvec->z;
	return std::make_tuple(GEi, GEj);
}

std::tuple<vec3d, vec3d> Lubrication::calcGE_squeeze_tangential()
{
	/*
	 * lubrication_model = 2
	 * upto log(1/xi) level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(E+tE).nvec;
	 * GE2 = (nvecnvec:E)*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(E+tE).nvec;
	 */
	double nnE;
	if (sys->p.cross_shear) {
		double costheta, sintheta;
		std::tie(costheta, sintheta) = sys->getCosSinShearAngle();
		// nnE = nynz;
		nnE = costheta*nvec->x*nvec->z+sintheta*nvec->y*nvec->z;
	} else {
		nnE = nvec->x*nvec->z;
	}
	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE;
	vec3d GEi, GEj;
	if (!sys->p.cross_shear) {
		GEi.x =  cGE_i*nvec->x+YG0_YG2*nvec->z;
		GEi.y =  cGE_i*nvec->y;
		GEi.z =  cGE_i*nvec->z+YG0_YG2*nvec->x;
		GEj.x =  cGE_j*nvec->x+YG1_YG3*nvec->z;
		GEj.y =  cGE_j*nvec->y;
		GEj.z =  cGE_j*nvec->z+YG1_YG3*nvec->x;
	} else {
		double costheta, sintheta;
		std::tie(costheta, sintheta) = sys->getCosSinShearAngle();
		double costheta_nx_sintheta_ny = costheta*nvec->x+sintheta*nvec->y;
		GEi.x =  cGE_i*nvec->x+YG0_YG2*costheta*nvec->z;
		GEi.y =  cGE_i*nvec->y+YG0_YG2*sintheta*nvec->z;
		GEi.z =  cGE_i*nvec->z+YG0_YG2*costheta_nx_sintheta_ny;
		GEj.x =  cGE_j*nvec->x+YG1_YG3*costheta*nvec->z;
		GEj.y =  cGE_j*nvec->y+YG1_YG3*sintheta*nvec->z;
		GEj.z =  cGE_j*nvec->z+YG1_YG3*costheta_nx_sintheta_ny;
	}
	return std::make_tuple(GEi, GEj);
}

std::tuple<vec3d, vec3d, vec3d, vec3d> Lubrication::calcGEHE_squeeze_tangential()
{
	/*
	 * lubrication_model = 2
	 * upto log(1/xi) level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(E+tE).nvec;
	 * GE2 = (nvecnvec:E)*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(E+tE).nvec;
	 */

	double nxnx = nvec->x*nvec->x;
 	double nxny = nvec->x*nvec->y;
 	double nxnz = nvec->x*nvec->z;
 	double nynz = nvec->y*nvec->z;
 	double nyny = nvec->y*nvec->y;
 	double nznz = nvec->z*nvec->z;
	double nnE;
 	if (sys->p.cross_shear) {
 		double costheta, sintheta;
 		std::tie(costheta, sintheta) = sys->getCosSinShearAngle();
 		// nnE = nynz;
 		nnE = costheta*nxnz+sintheta*nynz;
 	} else {
 		nnE = nxnz;
 	}

	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE;
	double cHE_i =  YH[0]+YH[2];
	double cHE_j =  YH[1]+YH[3];
	vec3d GEi, GEj, HEi, HEj;
	if (!sys->p.cross_shear) {
		GEi.x =  cGE_i*nvec->x+YG0_YG2*nvec->z;
		GEi.y =  cGE_i*nvec->y;
		GEi.z =  cGE_i*nvec->z+YG0_YG2*nvec->x;
		GEj.x =  cGE_j*nvec->x+YG1_YG3*nvec->z;
		GEj.y =  cGE_j*nvec->y;
		GEj.z =  cGE_j*nvec->z+YG1_YG3*nvec->x;
		double nxnx_nznz = nxnx-nznz;
		HEi.x =  cHE_i*nxny;
		HEi.y = -cHE_i*nxnx_nznz;
		HEi.z = -cHE_i*nynz;
		HEj.x =  cHE_j*nxny;
		HEj.y = -cHE_j*nxnx_nznz;
		HEj.z = -cHE_j*nynz;
	} else {
		double costheta, sintheta;
		std::tie(costheta, sintheta) = sys->getCosSinShearAngle();
		double costheta_nx_sintheta_ny = costheta*nvec->x+sintheta*nvec->y;
		GEi.x =  cGE_i*nvec->x+YG0_YG2*costheta*nvec->z;
		GEi.y =  cGE_i*nvec->y+YG0_YG2*sintheta*nvec->z;
		GEi.z =  cGE_i*nvec->z+YG0_YG2*costheta_nx_sintheta_ny;
		GEj.x =  cGE_j*nvec->x+YG1_YG3*costheta*nvec->z;
		GEj.y =  cGE_j*nvec->y+YG1_YG3*sintheta*nvec->z;
		GEj.z =  cGE_j*nvec->z+YG1_YG3*costheta_nx_sintheta_ny;
		double nyny_nznz = nyny-nznz;
		double nxnx_nznz = nxnx-nznz;
		HEi.x =  cHE_i*( costheta*nxny      + sintheta*nyny_nznz);
		HEi.y = -cHE_i*( costheta*nxnx_nznz + sintheta*nxny);
		HEi.z =  cHE_i*(-costheta*nynz      + sintheta*nxnz);
		HEj.x =  cHE_j*( costheta*nxny      + sintheta*nyny_nznz);
		HEj.y = -cHE_j*( costheta*nxnx_nznz + sintheta*nxny);
		HEj.z =  cHE_j*(-costheta*nynz      + sintheta*nxnz);
	}
	return std::make_tuple(GEi, GEj, HEi, HEj);
}

// for FT/UW version
struct ODBlock Lubrication::RFU_ODBlock_squeeze()
{
	double n0n0 = nvec->x*nvec->x;
	double n0n1 = nvec->x*nvec->y;
	double n0n2 = nvec->x*nvec->z;
	double n1n1 = nvec->y*nvec->y;
	double n1n2 = nvec->y*nvec->z;
	double n2n2 = nvec->z*nvec->z;

	struct ODBlock block;
	// column 0
	block.col0[0] = XA[1]*n0n0;
	block.col0[1] = XA[1]*n0n1;
	block.col0[2] = XA[1]*n0n2;
	block.col0[3] = 0;
	block.col0[4] =  0;
	// column 1
	block.col1[0] =  XA[1]*n1n1;
	block.col1[1] = XA[1]*n1n2;
	block.col1[2] = 0;
	// column 2
	block.col2[0] =  XA[1]*n2n2;
	// column 3
	block.col3[0] =  0;
	block.col3[1] =  0;
	block.col3[2] =  0;
	block.col3[3] =  0;
	block.col3[4] =  0;
	// column 4
	block.col4[0] =  0;
	block.col4[1] =  0;
	block.col4[2] = -0;
	// column 5
	block.col5[0] =  0;
	return block;
}

// for FT/UW version
struct ODBlock Lubrication::RFU_ODBlock_squeeze_tangential()
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

// Diagonal Blocks Terms, F/U version
std::pair<struct DBlock, struct DBlock> Lubrication::RFU_DBlocks_squeeze()
{
	struct DBlock b0;
	struct DBlock b1;

	double n0n0 = nvec->x*nvec->x;
	double n0n1 = nvec->x*nvec->y;
	double n0n2 = nvec->x*nvec->z;
	double n1n1 = nvec->y*nvec->y;
	double n1n2 = nvec->y*nvec->z;
	double n2n2 = nvec->z*nvec->z;

	// (*,0)
	b0.col0[0] = XA[0]*n0n0; // 00 element of the dblock
	b0.col0[1] = XA[0]*n0n1;           // 10
	b0.col0[2] = XA[0]*n0n2;           // 20
	b0.col0[3] = 0;                   // 40
	b0.col0[4] = 0;                   // 50
	// (*,1)
	b0.col1[0] = XA[0]*n1n1; // 11
	b0.col1[1] = XA[0]*n1n2;           // 21
	b0.col1[2] = 0;                   // 51
	// (*,2)
	b0.col2[0] = XA[0]*n2n2 ; // 22
	// (*,3)
	b0.col3[0] = 0;                 // 33
	b0.col3[1] = 0;                     // 43
	b0.col3[2] = 0;                     // 53
	// (*,4)
	b0.col4[0] = 0;                 // 44
	b0.col4[1] = 0;                     // 54
	// (*,5)
	b0.col5[0] = 0;                 // 55

	 // (*,0)
 	b1.col0[0] = XA[3]*n0n0; // 00 element of the dblock
 	b1.col0[1] = XA[3]*n0n1;           // 10
 	b1.col0[2] = XA[3]*n0n2;           // 20
 	b1.col0[3] = 0;                   // 40
 	b1.col0[4] = 0;                   // 50
 	// (*,1)
 	b1.col1[0] = XA[3]*n1n1; // 11
 	b1.col1[1] = XA[3]*n1n2;           // 21
 	b1.col1[2] = 0;                   // 51
 	// (*,2)
 	b1.col2[0] = XA[3]*n2n2; // 22
 	// (*,3)
 	b1.col3[0] = 0;                 // 33
 	b1.col3[1] = 0;                     // 43
 	b1.col3[2] = 0;                     // 53
 	// (*,4)
 	b1.col4[0] = 0;                 // 44
 	b1.col4[1] = 0;                     // 54
 	// (*,5)
 	b1.col5[0] = 0;                 // 55

	return std::make_pair(b0, b1);
}

// Diagonal Blocks Terms, FT/UW version
std::pair<struct DBlock, struct DBlock> Lubrication::RFU_DBlocks_squeeze_tangential()
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



// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jeffrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
void Lubrication::pairVelocityStresslet(const vec3d& vi, const vec3d& vj,
										const vec3d& oi, const vec3d& oj,
										StressTensor& stresslet_i,
										StressTensor& stresslet_j)
{
	/*
	 * (xx, xy, xz, yz, yy, zz)
	 *
	 * S_{11}^{GUX}+S_{12}^{GUX}
	 *         = - vec{n}.(XG11*v1+XG12*v2)*(ninj-(1/3)*delta_{ij})
	 *
	 */
	double nxnx = nvec->x*nvec->x;
	double nxny = nvec->x*nvec->y;
	double nxnz = nvec->x*nvec->z;
	double nynz = nvec->y*nvec->z;
	double nyny = nvec->y*nvec->y;
	double nznz = nvec->z*nvec->z;

	double cXG_i = -dot(nvec, XG[0]*vi+XG[1]*vj);
	double cXG_j = -dot(nvec, XG[2]*vi+XG[3]*vj);
	StressTensor XGU_i(nxnx, nxny, nxnz, nynz, nyny, nznz);
	StressTensor XGU_j = XGU_i;
	XGU_i *= cXG_i;
	XGU_j *= cXG_j;
	stresslet_i = XGU_i;
	stresslet_j = XGU_j;
	if (sys->p.lubrication_model == 1 || sys->p.lubrication_model == 3) {
		return;
	}
	StressTensor YGU_i;
	StressTensor YGU_j;
	double cYGUi_xx = nvec->x*vi.x+nvec->x*vi.x-2*nxnx*dot(nvec, vi);
	double cYGUj_xx = nvec->x*vj.x+nvec->x*vj.x-2*nxnx*dot(nvec, vj);
	YGU_i.elm[0] = -YG[0]*cYGUi_xx-YG[1]*cYGUj_xx;
	YGU_j.elm[0] = -YG[2]*cYGUi_xx-YG[3]*cYGUj_xx;

	double cYGUi_xy = nvec->x*vi.y+nvec->y*vi.x-2*nxny*dot(nvec, vi);
	double cYGUj_xy = nvec->x*vj.y+nvec->y*vj.x-2*nxny*dot(nvec, vj);
	YGU_i.elm[1] = -YG[0]*cYGUi_xy-YG[1]*cYGUj_xy;
	YGU_j.elm[1] = -YG[2]*cYGUi_xy-YG[3]*cYGUj_xy;

	double cYGUi_xz = nvec->x*vi.z+nvec->z*vi.x-2*nxnz*dot(nvec, vi);
	double cYGUj_xz = nvec->x*vj.z+nvec->z*vj.x-2*nxnz*dot(nvec, vj);
	YGU_i.elm[2] = -YG[0]*cYGUi_xz-YG[1]*cYGUj_xz;
	YGU_j.elm[2] = -YG[2]*cYGUi_xz-YG[3]*cYGUj_xz;

	double cYGUi_yz = nvec->y*vi.z+nvec->z*vi.y-2*nynz*dot(nvec, vi);
	double cYGUj_yz = nvec->y*vj.z+nvec->z*vj.y-2*nynz*dot(nvec, vj);
	YGU_i.elm[3] = -YG[0]*cYGUi_yz-YG[1]*cYGUj_yz;
	YGU_j.elm[3] = -YG[2]*cYGUi_yz-YG[3]*cYGUj_yz;

	double cYGUi_yy = nvec->y*vi.y+nvec->y*vi.y-2*nyny*dot(nvec, vi);
	double cYGUj_yy = nvec->y*vj.y+nvec->y*vj.y-2*nyny*dot(nvec, vj);

	YGU_i.elm[4] = -YG[0]*cYGUi_yy-YG[1]*cYGUj_yy;
	YGU_j.elm[4] = -YG[2]*cYGUi_yy-YG[3]*cYGUj_yy;

	double cYGUi_zz = nvec->z*vi.z+nvec->z*vi.z-2*nznz*dot(nvec, vi);
	double cYGUj_zz = nvec->z*vj.z+nvec->z*vj.z-2*nznz*dot(nvec, vj);

	YGU_i.elm[5] = -YG[0]*cYGUi_zz-YG[1]*cYGUj_zz;
	YGU_j.elm[5] = -YG[2]*cYGUi_zz-YG[3]*cYGUj_zz;

	stresslet_i += YGU_i;
	stresslet_j += YGU_j;

	StressTensor YHO_i;
	StressTensor YHO_j;
	double cYHOi_xx = nxnz*oi.y-nxny*oi.z;
	double cYHOj_xx = nxnz*oj.y-nxny*oj.z;
	YHO_i.elm[0] = -2*(YM[0]*cYHOi_xx+YM[1]*cYHOj_xx);
	YHO_j.elm[0] = -2*(YM[2]*cYHOi_xx+YM[3]*cYHOj_xx);

	double cYHOi_xy = nxnx*oi.z-nxnz*oi.x+nynz*oi.y-nyny*oi.z;
	double cYHOj_xy = nxnx*oj.z-nxnz*oj.x+nynz*oj.y-nyny*oj.z;
	YHO_i.elm[1] = -YM[0]*cYHOi_xy-YM[1]*cYHOj_xy;
	YHO_j.elm[1] = -YM[2]*cYHOi_xy-YM[3]*cYHOj_xy;

	double cYHOi_xz = nxny*oi.x-nxnx*oi.y+nznz*oi.y-nynz*oi.z;
	double cYHOj_xz = nxny*oj.x-nxnx*oj.y+nznz*oj.y-nynz*oj.z;
	YHO_i.elm[2] = -YM[0]*cYHOi_xz-YM[1]*cYHOj_xz;
	YHO_j.elm[2] = -YM[2]*cYHOi_xz-YM[3]*cYHOj_xz;

	double cYHOi_yz = nyny*oi.x-nxny*oi.y+nxnz*oi.z-nznz*oi.x;
	double cYHOj_yz = nyny*oj.x-nxny*oj.y+nxnz*oj.z-nznz*oj.x;
	YHO_i.elm[3] = -YM[0]*cYHOi_yz-YM[1]*cYHOj_yz;
	YHO_j.elm[3] = -YM[2]*cYHOi_yz-YM[3]*cYHOj_yz;

	double cYHOi_yy = nxny*oi.z-nynz*oi.x;
	double cYHOj_yy = nxny*oj.z-nynz*oj.x;
	YHO_i.elm[4] = -2*(YM[0]*cYHOi_yy+YM[1]*cYHOj_yy);
	YHO_j.elm[4] = -2*(YM[2]*cYHOi_yy+YM[3]*cYHOj_yy);

	double cYHOi_zz = nynz*oi.x-nxnz*oi.y;
	double cYHOj_zz = nynz*oj.x-nxnz*oj.y;
	YHO_i.elm[5] = -2*(YM[0]*cYHOi_zz+YM[1]*cYHOj_zz);
	YHO_j.elm[5] = -2*(YM[2]*cYHOi_zz+YM[3]*cYHOj_zz);

	stresslet_i += YHO_i;
	stresslet_j += YHO_j;
}

void Lubrication::pairStrainStresslet(StressTensor& stresslet_i,
									  StressTensor& stresslet_j)
{
	/**
		\brief The \f$ M:\hat{E}^{\infty} \f$ component of the stress.

		(Jeffrey 1992 notations)
		For a non-dimensionalized strain rate tensor \f$ \hat{E}^{\infty}_{ij} = \frac{1}{2}(\delta_{ia}\delta_{jb} + \delta_{ib}\delta_{ja}) \f$:

		\f$ (M_{\alpha\beta}:\hat{E}^{\infty})_{ij} = \frac{3}{2} X_{\alpha\beta} n_i n_j n_a n_b \\ \qquad \qquad \qquad + \frac{1}{2} Y_{\alpha\beta} ( -4 n_i n_j n_a n_b + n_j\delta_{ib}n_a + n_j\delta_{ia}n_b + n_i \delta_{jb}n_a + n_i \delta_{ja}n_b )\f$

		\b Note that we remove the explicit trace removal terms from Jeffrey's expressions. The stress tensors this method computes are \b NOT traceless.

		On return the stress on particle 1 is
		\f$ S_1 = (M_{11}+M_{12}):\hat{E}^{\infty} \f$

		On particle 2 it is
		\f$ S_2 = (M_{21}+M_{22}):\hat{E}^{\infty} \f$

	 */
	double nxnx = nvec->x*nvec->x;
 	double nxny = nvec->x*nvec->y;
 	double nxnz = nvec->x*nvec->z;
 	double nynz = nvec->y*nvec->z;
 	double nyny = nvec->y*nvec->y;
 	double nznz = nvec->z*nvec->z;
	double nnE;
 	if (sys->p.cross_shear) {
 		double costheta, sintheta;
 		std::tie(costheta, sintheta) = sys->getCosSinShearAngle();
 		// nnE = nynz;
 		nnE = costheta*nxnz+sintheta*nynz;
 	} else {
 		nnE = nxnz;
 	}
	double cXM_i = (3.0/2)*(XM[0]+XM[1])*nnE;
	double cXM_j = (3.0/2)*(XM[2]+XM[3])*nnE;
	StressTensor XME_i(nxnx, nxny, nxnz, nynz, nyny, nznz);
	StressTensor XME_j = XME_i;
	XME_i *= cXM_i;
	XME_j *= cXM_j;
	stresslet_i = XME_i;
	stresslet_j = XME_j;
	if (sys->p.lubrication_model == 1 || sys->p.lubrication_model == 3) {
		return;
	}
	double cYM_i = (1.0/2)*(YM[0]+YM[1]);
	double cYM_j = (1.0/2)*(YM[2]+YM[3]);

	StressTensor YME_i;
	if (!sys->p.cross_shear) {
		YME_i.elm[0] = 2*nxnz     -4*nxnx*nnE;
		YME_i.elm[1] =   nynz     -4*nxny*nnE;
		YME_i.elm[2] =   nxnx+nznz-4*nxnz*nnE;
		YME_i.elm[3] =   nxny     -4*nynz*nnE;
		YME_i.elm[4] =            -4*nyny*nnE;
		YME_i.elm[5] = 2*nxnz     -4*nznz*nnE;
	}
	if (sys->p.cross_shear) {
		double costheta, sintheta;
		std::tie( costheta, sintheta ) = sys->getCosSinShearAngle();
		YME_i.elm[0] = 2*costheta*nxnz                             -4*nxnx*nnE;
		YME_i.elm[1] =   costheta*nynz        +sintheta*nxnz       -4*nxny*nnE;
		YME_i.elm[2] =   costheta*(nxnx+nznz) +sintheta*nxny       -4*nxnz*nnE;
		YME_i.elm[3] =   costheta*nxny        +sintheta*(nyny+nznz)-4*nynz*nnE;
		YME_i.elm[4] = 2*sintheta*nynz                             -4*nyny*nnE;
		YME_i.elm[5] = 2*costheta*nxnz        +2*sintheta*nynz     -4*nznz*nnE;
	}

	StressTensor YME_j = YME_i;
	YME_i *= cYM_i;
	YME_j *= cYM_j;
	stresslet_i += YME_i;
	stresslet_j += YME_j;
}

/* Lubriction force between two particles is calculated.
 * Note that only the Brownian component of the velocity is NOT included here (IS THAT TRUE?).
 *     @@@@ na_velocity includes Brownian component(??)
 * This part is used for ouput data.
 * lubforce_p1 = -lubforce_p0
 *
 */
void Lubrication::calcPairwiseForce()
{
	/*
	 *  First: -A*(U-Uinf) term
	 */
	/* Eq. (1.6a) in Jeffrey&Onishi 1984
	 * A_{ij}^{ab} = XA_{ab}ni*nj + YA_{ab}(del_{ij}-ni*nj)
	 * B~_{ji}^{ab} = YB_{ba}epsilon_{jik} nk
	 *
	 */
	double sr = sys->get_shear_rate();
	vec3d vi(sys->na_velocity[p0]);
	vec3d vj(sys->na_velocity[p1]);
	if (sys->p.lubrication_model == 1 || sys->p.lubrication_model == 3) {
		calcXFunctions();
	} else if (sys->p.lubrication_model == 2) {
		calcXYFunctions();
	}
	/* XAU_i */
	lubforce_p0 = -dot(XA[0]*vi+XA[1]*vj, nvec)*(*nvec);
	if (sys->p.lubrication_model == 2) {
		vec3d oi(sys->na_ang_velocity[p0]);
		vec3d oj(sys->na_ang_velocity[p1]);
		/* YAU_i */
		lubforce_p0 += -YA[0]*(vi-(*nvec)*dot(nvec, vi)) - YA[1]*(vj-(*nvec)*dot(nvec, vj));
		/* YBO_i */
		lubforce_p0 += -YB[0]*cross(nvec, oi)            - YB[2]*cross(nvec, oj);
	}
	if (!sys->zero_shear) {
		vec3d GEi, GEj;
		if (sys->p.lubrication_model == 1) {
			std::tie(GEi, GEj) = calcGE_squeeze();
		} else {
			std::tie(GEi, GEj) = calcGE_squeeze_tangential();
		}
		/* XGE_i */
		lubforce_p0 += sr*GEi;
	}
	return;
}

void Lubrication::addStressME()
{
	/** Stress from the M*Einf term (see Jeffrey 1991)
	*/
	StressTensor stresslet_ME_i;
	StressTensor stresslet_ME_j;
	if (!sys->zero_shear) {
		pairStrainStresslet(stresslet_ME_i, stresslet_ME_j);
		double sr = sys->get_shear_rate();
		stresslet_ME_i *= sr;
		stresslet_ME_j *= sr;
	}
	sys->lubstress[p0] += stresslet_ME_i;
	sys->lubstress[p1] += stresslet_ME_j;
}

void Lubrication::addStressesGU()
{
	/*
	 *  First: -G*(U-Uinf) term
	 */
	StressTensor stresslet_hydro_GU_i;
	StressTensor stresslet_hydro_GU_j;
	pairVelocityStresslet(sys->vel_hydro[p0], sys->vel_hydro[p1],
						  sys->ang_vel_hydro[p0], sys->ang_vel_hydro[p1],
						  stresslet_hydro_GU_i, stresslet_hydro_GU_j);

	// Add term G*V_cont
	/* [note]
	 * We must not check for interaction->contact.is_active() condition,
	 * because the contact velocities of p0 and p1 can be non-zero
	 * even when they are not in contact.
	 */
	StressTensor stresslet_contact_GU_i;
	StressTensor stresslet_contact_GU_j;
	pairVelocityStresslet(sys->vel_contact[p0], sys->vel_contact[p1],
						  sys->ang_vel_contact[p0], sys->ang_vel_contact[p1],
						  stresslet_contact_GU_i, stresslet_contact_GU_j);
	sys->contactstressGU[p0] += stresslet_contact_GU_i;
	sys->contactstressGU[p1] += stresslet_contact_GU_j;
	// Add term G*V_repulsive
	if (sys->repulsiveforce) {
		StressTensor stresslet_repulsive_GU_i;
		StressTensor stresslet_repulsive_GU_j;
		pairVelocityStresslet(sys->vel_repulsive[p0], sys->vel_repulsive[p1],
							  sys->ang_vel_repulsive[p0], sys->ang_vel_repulsive[p1],
							  stresslet_repulsive_GU_i, stresslet_repulsive_GU_j);
		sys->repulsivestressGU[p0] += stresslet_repulsive_GU_i;
		sys->repulsivestressGU[p1] += stresslet_repulsive_GU_j;
	}
	// Add term G*V_brownian
	if (sys->brownian) {
		StressTensor stresslet_brownian_GU_i;
		StressTensor stresslet_brownian_GU_j;
		pairVelocityStresslet(sys->vel_brownian[p0], sys->vel_brownian[p1],
							  sys->ang_vel_brownian[p0], sys->ang_vel_brownian[p1],
							  stresslet_brownian_GU_i, stresslet_brownian_GU_j);
		sys->brownianstressGU[p0] += stresslet_brownian_GU_i;
		sys->brownianstressGU[p1] += stresslet_brownian_GU_j;
	}
	if (sys->vel_hydro_from_fixed.size()) {
		StressTensor stresslet_fixed_GU_i;
		StressTensor stresslet_fixed_GU_j;
		if (p1 < sys->np_mobile) { // p0 and p1 mobile
			pairVelocityStresslet(sys->vel_hydro_from_fixed[p0], sys->vel_hydro_from_fixed[p1],
								  sys->ang_vel_hydro_from_fixed[p0], sys->ang_vel_hydro_from_fixed[p1],
								  stresslet_fixed_GU_i, stresslet_fixed_GU_j);
		} else if (p0 >= sys->np_mobile) { // p0 and p1 fixed
			pairVelocityStresslet(sys->na_velocity[p0], sys->na_velocity[p1],
								  sys->na_ang_velocity[p0], sys->na_ang_velocity[p1],
								  stresslet_fixed_GU_i, stresslet_fixed_GU_j);
		} else { // p0 mobile and p1 fixed
			pairVelocityStresslet(sys->vel_hydro_from_fixed[p0], sys->na_velocity[p1],
								  sys->ang_vel_hydro_from_fixed[p0], sys->na_ang_velocity[p1],
								  stresslet_fixed_GU_i, stresslet_fixed_GU_j);
		}
		sys->hydrofromfixedstressGU[p0] += stresslet_fixed_GU_i;
		sys->hydrofromfixedstressGU[p1] += stresslet_fixed_GU_j;
	}
}

void Lubrication::updateResistanceCoeff()
{
	if (interaction->reduced_gap > 0) {
		double coeff = 1/(interaction->reduced_gap+sys->p.lub_reduce_parameter);
		setResistanceCoeff(coeff, log(coeff));
	} else {
		setResistanceCoeff(sys->lub_coeff_contact, sys->log_lub_coeff_contact_tan_total);
	}
}
