//
//  Lubrication.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2016 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <tuple>
#include <stdexcept>
#include <sstream>
#include "Lubrication.h"
#include "LubricationFunctions.h"
#include "Interaction.h"
#include "System.h"
using namespace std;

Lubrication::Lubrication():
_active(false),
range(0),
p0(0),
p1(0),
p0_6(0),
p1_6(0),
lub_coeff(0),
log_lub_coeff(0),
a0(0),
a1(0),
ro(0),
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
lambda(0),
invlambda(0),
lambda_square(0),
lambda_cubic(0),
lambda_p_1(0),
lambda_p_1_square(0),
lambda_p_1_cubic(0),
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
tangential(true)
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
	if (!sys->lubrication) {
		throw runtime_error("Lubrication::init called but sys->lubrication is set to false.");
	}
	if (sys->p.lubrication_model == "normal") {
		tangential = false;
	} else if (sys->p.lubrication_model == "tangential") {
		tangential = true;
	} else {
		ostringstream error_str;
		error_str << "Unknown lubrication_model " << sys->p.lubrication_model << endl;
		throw runtime_error(error_str.str());
	}
}

void Lubrication::setParticleData()
{
	/**
	 	Set things having to do with particle index and radii.
		No position, velocity or force data is used here.
		**/
	tie(p0, p1) = interaction->get_par_num();
	p0_6 = 6*p0;
	p1_6 = 6*p1;
	calcLubConstants();
	range = sys->calcLubricationRange(p0, p1);
}

void Lubrication::activate()
{
	_active = true;
	sys->declareResistance(p0, p1);
}

void Lubrication::deactivate()
{
	_active = false;
	sys->eraseResistance(p0, p1);
}

void Lubrication::updateActivationState()
{
	bool in_range = interaction->separation_distance() < range && interaction->separation_distance() > ro;
	if (!is_active() && in_range) {
		activate();
	}
	if (is_active() && !in_range) {
		deactivate();
	}
}

void Lubrication::setResistanceCoeff(double lub_coeff_, double log_lub_coeff_)
{
	lub_coeff = lub_coeff_; // normal
	log_lub_coeff = log_lub_coeff_; // tangential
	if (!tangential){
		calcXFunctions();
	} else {
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
	a0 = sys->radius[p0];
	a1 = sys->radius[p1];
	ro = a0+a1;
	ro_12 = ro/2;
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
	for (int j=0; j<4; j++) {
		XM[j] = cXM[j]*lub_coeff;
		YM[j] = cYM[j]*log_lub_coeff;
	}
}

std::tuple<vec3d, vec3d> Lubrication::calcGE_squeeze(const Sym2Tensor& E_inf) const
{
	/* NOTE:
	 * Calculation of XG and YG needs to be done before that.
	 *
	 * mode normal
	 * 1/xi level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21)*nvec
	 * GE2 = (nvecnvec:E)*(XG12+XG22)*nvec
	 */
	double nnE;
	nnE = dot(nvec, dot(E_inf, *nvec));
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

std::tuple<vec3d, vec3d> Lubrication::calcGE_squeeze_tangential(const Sym2Tensor& E_inf)  const
{
	/*
	* mode normal+tangential
	 * upto log(1/xi) level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(E+tE).nvec;
	 * GE2 = (nvecnvec:E)*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(E+tE).nvec;
	 */
	double nnE = dot(nvec, dot(E_inf, *nvec));
	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE;
	vec3d Einf_nvec = dot(E_inf, *nvec);
	vec3d GEi = cGE_i*(*nvec) + (YG0_YG2*2)*Einf_nvec;
	vec3d GEj = cGE_j*(*nvec) + (YG1_YG3*2)*Einf_nvec;

	return std::make_tuple(GEi, GEj);
}

std::tuple<vec3d, vec3d, vec3d, vec3d> Lubrication::calcGEHE_squeeze_tangential(const Sym2Tensor& E_inf) const
{
	/*
	 * mode normal+tangential
	 * upto log(1/xi) level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(E+tE).nvec;
	 * GE2 = (nvecnvec:E)*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(E+tE).nvec;
	 */
	double nnE = dot(nvec, dot(E_inf, *nvec));
	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE;
	double cHE_i = YH[0]+YH[2];
	double cHE_j = YH[1]+YH[3];
	vec3d Einf_nvec = dot(E_inf, *nvec);
	vec3d GEi = cGE_i*(*nvec) + (YG0_YG2*2)*Einf_nvec;
	vec3d GEj = cGE_j*(*nvec) + (YG1_YG3*2)*Einf_nvec;
	vec3d nvec_Einf = dot(*nvec, E_inf);
	vec3d nvec_Einf_x_nvec = cross(nvec_Einf, (*nvec));
	vec3d HEi = -2*cHE_i*nvec_Einf_x_nvec;
	vec3d HEj = -2*cHE_j*nvec_Einf_x_nvec;

	return std::make_tuple(GEi, GEj, HEi, HEj);
}

// for FT/UW version
struct ODBlock Lubrication::RFU_ODBlock_squeeze() const
{
	struct ODBlock block;
	// column 0
	block.col0[0] = XA[1]*nvec->x*nvec->x;
	block.col0[1] = XA[1]*nvec->x*nvec->y;
	block.col0[2] = XA[1]*nvec->x*nvec->z;
	block.col0[3] = 0;
	block.col0[4] = 0;
	// column 1
	block.col1[0] = XA[1]*nvec->y*nvec->y;
	block.col1[1] = XA[1]*nvec->y*nvec->z;
	block.col1[2] = 0;
	// column 2
	block.col2[0] = XA[1]*nvec->z*nvec->z;
	// column 3
	block.col3[0] = 0;
	block.col3[1] = 0;
	block.col3[2] = 0;
	block.col3[3] = 0;
	block.col3[4] = 0;
	// column 4
	block.col4[0] = 0;
	block.col4[1] = 0;
	block.col4[2] = -0;
	// column 5
	block.col5[0] = 0;
	return block;
}

// for FT/UW version
struct ODBlock Lubrication::RFU_ODBlock_squeeze_tangential() const
{
	double XA1_YA1 = XA[1]-YA[1];
	struct ODBlock block;
	// column 0
	block.col0[0] =  YA[1] + XA1_YA1*nvec->x*nvec->x;
	block.col0[1] =  XA1_YA1*nvec->x*nvec->y;
	block.col0[2] =  XA1_YA1*nvec->x*nvec->z;
	block.col0[3] = -YB[2]*nvec->z;
	block.col0[4] =  YB[2]*nvec->y;
	// column 1
	block.col1[0] =  YA[1] + XA1_YA1*nvec->y*nvec->y;
	block.col1[1] =  XA1_YA1*nvec->y*nvec->z;
	block.col1[2] = -YB[2]*nvec->x;
	// column 2
	block.col2[0] =  YA[1] + XA1_YA1*nvec->z*nvec->z;
	// column 3
	block.col3[0] =  YB[1]*nvec->z;
	block.col3[1] = -YB[1]*nvec->y;
	block.col3[2] =  YC[1]*(1-nvec->x*nvec->x);

	block.col3[3] = -YC[1]*nvec->x*nvec->y;
	block.col3[4] = -YC[1]*nvec->x*nvec->z;
	// column 4
	block.col4[0] =  YB[1]*nvec->x;
	block.col4[1] =  YC[1]*(1-nvec->y*nvec->y);
	block.col4[2] = -YC[1]*nvec->y*nvec->z;
	// column 5
	block.col5[0] =  YC[1]*(1-nvec->z*nvec->z);
	return block;
}

// Diagonal Blocks Terms, F/U version
std::pair<struct DBlock, struct DBlock> Lubrication::RFU_DBlocks_squeeze() const
{
	struct DBlock b0;
	struct DBlock b1;
	// (*,0)
	b0.col0[0] = XA[0]*nvec->x*nvec->x; // 00 element of the dblock
	b0.col0[1] = XA[0]*nvec->x*nvec->y;           // 10
	b0.col0[2] = XA[0]*nvec->x*nvec->z;           // 20
	b0.col0[3] = 0;                   // 40
	b0.col0[4] = 0;                   // 50
	// (*,1)
	b0.col1[0] = XA[0]*nvec->y*nvec->y; // 11
	b0.col1[1] = XA[0]*nvec->y*nvec->z; // 21
	b0.col1[2] = 0;                   // 51
	// (*,2)
	b0.col2[0] = XA[0]*nvec->z*nvec->z ; // 22
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
 	b1.col0[0] = XA[3]*nvec->x*nvec->x; // 00 element of the dblock
 	b1.col0[1] = XA[3]*nvec->x*nvec->y;           // 10
 	b1.col0[2] = XA[3]*nvec->x*nvec->z;           // 20
 	b1.col0[3] = 0;                   // 40
 	b1.col0[4] = 0;                   // 50
 	// (*,1)
 	b1.col1[0] = XA[3]*nvec->y*nvec->y; // 11
	b1.col1[1] = XA[3]*nvec->y*nvec->z; // 21
 	b1.col1[2] = 0;                   // 51
 	// (*,2)
	b1.col2[0] = XA[3]*nvec->z*nvec->z; // 22
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
std::pair<struct DBlock, struct DBlock> Lubrication::RFU_DBlocks_squeeze_tangential() const
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
	double XA0mYA0 = XA[0]-YA[0];
	// (*,0)
	b0.col0[0] =  XA[0]*n0n0+YA[0]*one_n0n0; // 00 element of the dblock
	b0.col0[1] = XA0mYA0*n0n1;           // 10
	b0.col0[2] = XA0mYA0*n0n2;           // 20
	b0.col0[3] = -YB[0]*nvec->z;                   // 40
	b0.col0[4] =  YB[0]*nvec->y;                   // 50
	// (*,1)
	b0.col1[0] =  XA[0]*n1n1+YA[0]*one_n1n1; // 11
	b0.col1[1] = XA0mYA0*n1n2;           // 21
	b0.col1[2] = -YB[0]*nvec->x;                   // 51
	// (*,2)
	b0.col2[0] =  XA[0]*n2n2+YA[0]*one_n2n2; // 22
	// (*,3)
	b0.col3[0] =  YC[0]*one_n0n0;                 // 33
	b0.col3[1] = -YC[0]*n0n1;                     // 43
	b0.col3[2] = -YC[0]*n0n2;                     // 53
	// (*,4)
	b0.col4[0] =  YC[0]*one_n1n1;                 // 44
	b0.col4[1] = -YC[0]*n1n2;                     // 54
	// (*,5)
	b0.col5[0] =  YC[0]*one_n2n2;                 // 55

	double XA3mYA3 = XA[3]-YA[3];
	 // (*,0)
 	b1.col0[0] =  XA[3]*n0n0+YA[3]*one_n0n0; // 00 element of the dblock
 	b1.col0[1] =  XA3mYA3*n0n1;           // 10
 	b1.col0[2] =  XA3mYA3*n0n2;           // 20
 	b1.col0[3] = -YB[3]*nvec->z;                   // 40
 	b1.col0[4] =  YB[3]*nvec->y;                   // 50
 	// (*,1)
 	b1.col1[0] =  XA[3]*n1n1+YA[3]*one_n1n1; // 11
 	b1.col1[1] =  XA3mYA3*n1n2;           // 21
 	b1.col1[2] = -YB[3]*nvec->x;                   // 51
 	// (*,2)
 	b1.col2[0] =  XA[3]*n2n2+YA[3]*one_n2n2; // 22
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
void Lubrication::addGUStresslet(const vec3d& vi, const vec3d& vj,
                                 const vec3d& oi, const vec3d& oj,
                                 Sym2Tensor& stresslet_i,
                                 Sym2Tensor& stresslet_j) const
{
	/*
	 * (xx, xy, xz, yz, yy, zz)
	 *
	 * S_{11}^{GUX}+S_{12}^{GUX}
	 *         = - vec{n}.(XG11*v1+XG12*v2)*(ninj-(1/3)*delta_{ij})
	 *
	 */
	double nvec_vi = dot(nvec, vi);
	double nvec_vj = dot(nvec, vj);
	Sym2Tensor nvec_nvec = outer(*nvec);
	stresslet_i += -(XG[0]*nvec_vi+XG[1]*nvec_vj)*nvec_nvec; // XGU_i
	stresslet_j += -(XG[2]*nvec_vi+XG[3]*nvec_vj)*nvec_nvec; // XGU_j
	if (!tangential) {
		return;
	}
	Sym2Tensor tmp_i = outer_sym(*nvec, vi) - nvec_vi*nvec_nvec;
	Sym2Tensor tmp_j = outer_sym(*nvec, vj) - nvec_vj*nvec_nvec;
	stresslet_i += -2*(YG[0]*tmp_i+YG[1]*tmp_j); // YGU_i
	stresslet_j += -2*(YG[2]*tmp_i+YG[3]*tmp_j); // YGU_j
	Sym2Tensor tmp2_i = outer_sym(*nvec, cross(oi, *nvec));
	Sym2Tensor tmp2_j = outer_sym(*nvec, cross(oj, *nvec));
	stresslet_i += -2*(YM[0]*tmp2_i+YM[1]*tmp2_j); // YHO_i
	stresslet_j += -2*(YM[2]*tmp2_i+YM[3]*tmp2_j); // YHO_j
}

void Lubrication::addMEStresslet(const Sym2Tensor& E_inf,
                                 Sym2Tensor& stresslet_i,
                                 Sym2Tensor& stresslet_j) const
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
	double nnE = dot(nvec, dot(E_inf, *nvec));
	Sym2Tensor nvec_nvec = outer(*nvec);
	double coeff = 1.5*nnE;
	stresslet_i += coeff*(XM[0]+XM[1])*nvec_nvec; // XME_i
	stresslet_j += coeff*(XM[2]+XM[3])*nvec_nvec; // XME_j
	if (tangential) {
		vec3d Einf_nvec = dot(E_inf, *nvec);
		Sym2Tensor nvec_Einf_nvec = outer_sym(*nvec, Einf_nvec);
		Sym2Tensor nvec_Einf_nvec__nnE_nvec_nvec = 2*(nvec_Einf_nvec-nnE*nvec_nvec);
		stresslet_i += (YM[0]+YM[1])*nvec_Einf_nvec__nnE_nvec_nvec; // YME_i
		stresslet_j += (YM[2]+YM[3])*nvec_Einf_nvec__nnE_nvec_nvec; // YME_j
	}
}

/* Lubriction force between two particles is calculated.
 * Note that only the Brownian component of the velocity is NOT included here (IS THAT TRUE?).
 *     @@@@ na_velocity includes Brownian component(??)
 * This part is used for ouput data.
 * lubforce_p1 = -lubforce_p0
 *
 */
vec3d Lubrication::getTotalForce() const
{
	/**
	* \brief The total lubrication force acting on p0.
	* NOTE: the velocities must be computed first by the System class.
	*/

	/* Eq. (1.6a) in Jeffrey&Onishi 1984
	 * A_{ij}^{ab} = XA_{ab}ni*nj + YA_{ab}(del_{ij}-ni*nj)
	 * B~_{ji}^{ab} = YB_{ba}epsilon_{jik} nk
	 *
	 */
	vec3d vi(sys->na_velocity[p0]);
	vec3d vj(sys->na_velocity[p1]);
	/* XAU_i */
	vec3d lubforce_p0 = -dot(XA[0]*vi+XA[1]*vj, nvec)*(*nvec);
	if (tangential) {
		vec3d oi(sys->na_ang_velocity[p0]);
		vec3d oj(sys->na_ang_velocity[p1]);
		/* YAU_i */
		lubforce_p0 += -YA[0]*(vi-(*nvec)*dot(nvec, vi)) - YA[1]*(vj-(*nvec)*dot(nvec, vj));
		/* YBO_i */
		lubforce_p0 += -YB[0]*cross(nvec, oi)            - YB[2]*cross(nvec, oj);
	}
	if (!sys->zero_shear) {
		vec3d GEi, GEj;
		if (!tangential) {
			std::tie(GEi, GEj) = calcGE_squeeze(sys->getEinfty());
		} else {
			std::tie(GEi, GEj) = calcGE_squeeze_tangential(sys->getEinfty());
		}
		/* XGE_i */
		lubforce_p0 += GEi;
	}
	return lubforce_p0;
}

void Lubrication::updateResistanceCoeff()
{
	if (interaction->get_reduced_gap() > 0) {
		double coeff = 1/(interaction->get_reduced_gap()+sys->p.lub_reduce_parameter);
		if (coeff > 1.) {
			setResistanceCoeff(coeff, log(coeff));
		} else {
			setResistanceCoeff(coeff, 0.);
		}
	}
}
