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
#include "PairwiseInteraction.h"

namespace Interactions {

namespace Lub {

Lubrication::Lubrication(PairwiseInteraction* interaction, const LubParams &params):
inter(interaction),
regularization_length(params.regularization_length),
coeffs(inter->a0, inter->a1)
{
	if (params.model == "normal") {
		tangential = false;
	} else if (params.model == "tangential") {
		tangential = true;
	} else {
		std::ostringstream error_str;
		error_str << "Unknown lubrication_model " << params.model << std::endl;
		throw std::runtime_error(error_str.str());
	}
	lub_coeff_min = 1/(params.max_gap+params.regularization_length);

	updateResistanceCoeff();
}

void Lubrication::setResistanceCoeff(double lub_coeff_, double log_lub_coeff_)
{
	lub_coeff = lub_coeff_; // normal
	log_lub_coeff = log_lub_coeff_; // tangential
	if (!tangential){
		calcXFunctions();
		calcXFunctionsStress();
	} else {
		calcXYFunctions();
		calcXYFunctionsStress();
	}
}

/*********************************
 *                                *
 *  Lubrication Forces Methods    *
 *                                *
 *********************************/

// Resistance functions
void Lubrication::calcXFunctions()
{
	for (int j=0; j<4; j++) {
		XA[j] = coeffs.cXA[j]*lub_coeff;
		XG[j] = coeffs.cXG[j]*lub_coeff;
	}
}

void Lubrication::calcXFunctionsStress()
{
	for (int j=0; j<4; j++) {
		XM[j] = coeffs.cXM[j]*lub_coeff;
	}
}

void Lubrication::calcXYFunctions()
{
	calcXFunctions();
	for (int j=0; j<4; j++) {
		YA[j] = coeffs.cYA[j]*log_lub_coeff;
		YB[j] = coeffs.cYB[j]*log_lub_coeff;
		YC[j] = coeffs.cYC[j]*log_lub_coeff;
		YG[j] = coeffs.cYG[j]*log_lub_coeff;
		YH[j] = coeffs.cYH[j]*log_lub_coeff;
	}
}

void Lubrication::calcXYFunctionsStress()
{
	for (int j=0; j<4; j++) {
		XM[j] = coeffs.cXM[j]*lub_coeff;
		YM[j] = coeffs.cYM[j]*log_lub_coeff;
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
	nnE = dot(inter->nvec, dot(E_inf, inter->nvec));
	double cGE_p0 = (XG[0]+XG[2])*nnE;
	double cGE_p1 = (XG[1]+XG[3])*nnE;
	vec3d GEi, GEj;
	GEi.x = cGE_p0*inter->nvec.x;
	GEi.y = cGE_p0*inter->nvec.y;
	GEi.z = cGE_p0*inter->nvec.z;
	GEj.x = cGE_p1*inter->nvec.x;
	GEj.y = cGE_p1*inter->nvec.y;
	GEj.z = cGE_p1*inter->nvec.z;
	return std::make_tuple(GEi, GEj);
}

std::tuple<vec3d, vec3d> Lubrication::calcGE_squeeze(const Sym2Tensor& E0, const Sym2Tensor& E1) const
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
	double nnE0, nnE1;
	nnE0 = dot(inter->nvec, dot(E0, inter->nvec));
	nnE1 = dot(inter->nvec, dot(E1, inter->nvec));
	double cGE_p0 = (XG[0]+XG[2])*nnE0;
	double cGE_p1 = (XG[1]+XG[3])*nnE1;
	vec3d GEi, GEj;
	GEi.x = cGE_p0*(inter->nvec).x;
	GEi.y = cGE_p0*(inter->nvec).y;
	GEi.z = cGE_p0*(inter->nvec).z;
	GEj.x = cGE_p1*(inter->nvec).x;
	GEj.y = cGE_p1*(inter->nvec).y;
	GEj.z = cGE_p1*(inter->nvec).z;
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
	double nnE = dot(inter->nvec, dot(E_inf, inter->nvec));
	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE;
	vec3d Einf_nvec = dot(E_inf, inter->nvec);
	vec3d GEi = cGE_i*(inter->nvec) + (YG0_YG2*2)*Einf_nvec;
	vec3d GEj = cGE_j*(inter->nvec) + (YG1_YG3*2)*Einf_nvec;

	return std::make_tuple(GEi, GEj);
}

std::tuple<vec3d, vec3d> Lubrication::calcGE_squeeze_tangential(const Sym2Tensor& E0, const Sym2Tensor& E1)  const
{
	/*
	 * mode normal+tangential
	 * upto log(1/xi) level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(E+tE).nvec;
	 * GE2 = (nvecnvec:E)*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(E+tE).nvec;
	 */
	double nnE0 = dot(inter->nvec, dot(E0, inter->nvec));
	double nnE1 = dot(inter->nvec, dot(E1, inter->nvec));
	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE0;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE1;
	vec3d E0_nvec = dot(E0, inter->nvec);
	vec3d E1_nvec = dot(E1, inter->nvec);
	vec3d GEi = cGE_i*inter->nvec + (YG0_YG2*2)*E0_nvec;
	vec3d GEj = cGE_j*inter->nvec + (YG1_YG3*2)*E1_nvec;
	
	return std::make_tuple(GEi, GEj);
}

std::tuple<vec3d, vec3d> Lubrication::calcGE(const Sym2Tensor& E_inf)  const
{
	if (!tangential) {
		return calcGE_squeeze(E_inf);
	} else {
		return calcGE_squeeze_tangential(E_inf);
	}
}

std::tuple<vec3d, vec3d> Lubrication::calcGE(const Sym2Tensor& E0, const Sym2Tensor& E1)  const
{
	if (!tangential) {
		return calcGE_squeeze(E0, E1);
	} else {
		return calcGE_squeeze_tangential(E0, E1);
	}
}

std::tuple<vec3d, vec3d, vec3d, vec3d> Lubrication::calcGEHE(const Sym2Tensor& E_inf) const
{
	/*
	 * mode normal+tangential
	 * upto log(1/xi) level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(E+tE).nvec;
	 * GE2 = (nvecnvec:E)*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(E+tE).nvec;
	 */
	double nnE = dot(inter->nvec, dot(E_inf, inter->nvec));
	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE;
	double cHE_i = YH[0]+YH[2];
	double cHE_j = YH[1]+YH[3];
	vec3d Einf_nvec = dot(E_inf, inter->nvec);
	vec3d GEi = cGE_i*inter->nvec + (YG0_YG2*2)*Einf_nvec;
	vec3d GEj = cGE_j*inter->nvec + (YG1_YG3*2)*Einf_nvec;
	vec3d nvec_Einf = dot(inter->nvec, E_inf);
	vec3d nvec_Einf_x_nvec = cross(nvec_Einf, inter->nvec);
	vec3d HEi = -2*cHE_i*nvec_Einf_x_nvec;
	vec3d HEj = -2*cHE_j*nvec_Einf_x_nvec;

	return std::make_tuple(GEi, GEj, HEi, HEj);
}

std::tuple<vec3d, vec3d, vec3d, vec3d> Lubrication::calcGEHE(const Sym2Tensor& E0, const Sym2Tensor& E1) const
{
	/*
	 * mode normal+tangential
	 * upto log(1/xi) level
	 *
	 * GE1 = (nvecnvec:E)*(XG11+XG21-2*(YG11+YG21))*nvec+(YG11+YG21)*(E+tE).nvec;
	 * GE2 = (nvecnvec:E)*(XG12+XG22-2*(YG12+YG22))*nvec+(YG12+YG22)*(E+tE).nvec;
	 */
	double nnE0 = dot(inter->nvec, dot(E0, inter->nvec));
	double nnE1 = dot(inter->nvec, dot(E1, inter->nvec));
	double YG0_YG2 = YG[0]+YG[2];
	double YG1_YG3 = YG[1]+YG[3];
	double cGE_i = (XG[0]+XG[2]-2*YG0_YG2)*nnE0;
	double cGE_j = (XG[1]+XG[3]-2*YG1_YG3)*nnE1;
	double cHE_i = YH[0]+YH[2];
	double cHE_j = YH[1]+YH[3];
	vec3d E0_nvec = dot(E0, inter->nvec);
	vec3d E1_nvec = dot(E1, inter->nvec);
	vec3d GEi = cGE_i*inter->nvec + (YG0_YG2*2)*E0_nvec;
	vec3d GEj = cGE_j*inter->nvec + (YG1_YG3*2)*E1_nvec;
	vec3d nvec_E0_x_nvec = cross(dot(inter->nvec, E0), inter->nvec);
	vec3d nvec_E1_x_nvec = cross(dot(inter->nvec, E1), inter->nvec);
	vec3d HEi = -2*cHE_i*nvec_E0_x_nvec;
	vec3d HEj = -2*cHE_j*nvec_E1_x_nvec;
	return std::make_tuple(GEi, GEj, HEi, HEj);
}

// for FT/UW version
struct ODBlock Lubrication::RFU_ODBlock_squeeze() const
{
	struct ODBlock block;
	const auto &nvec = inter->nvec;
	// column 0
	block.col0[0] = XA[1]*nvec.x*nvec.x;
	block.col0[1] = XA[1]*nvec.x*nvec.y;
	block.col0[2] = XA[1]*nvec.x*nvec.z;
	block.col0[3] = 0;
	block.col0[4] = 0;
	// column 1
	block.col1[0] = XA[1]*nvec.y*nvec.y;
	block.col1[1] = XA[1]*nvec.y*nvec.z;
	block.col1[2] = 0;
	// column 2
	block.col2[0] = XA[1]*nvec.z*nvec.z;
	// column 3
	block.col3[0] = 0;
	block.col3[1] = 0;
	block.col3[2] = 0;
	block.col3[3] = 0;
	block.col3[4] = 0;
	// column 4
	block.col4[0] = 0;
	block.col4[1] = 0;
	block.col4[2] = 0;
	// column 5
	block.col5[0] = 0;
	return block;
}

// for FT/UW version
struct ODBlock Lubrication::RFU_ODBlock_squeeze_tangential() const
{
	double XA1_YA1 = XA[1]-YA[1];
	struct ODBlock block;
	const auto &nvec = inter->nvec;

	// column 0
	block.col0[0] =  YA[1] + XA1_YA1*nvec.x*nvec.x;
	block.col0[1] =  XA1_YA1*nvec.x*nvec.y;
	block.col0[2] =  XA1_YA1*nvec.x*nvec.z;
	block.col0[3] = -YB[2]*nvec.z;
	block.col0[4] =  YB[2]*nvec.y;
	// column 1
	block.col1[0] =  YA[1] + XA1_YA1*nvec.y*nvec.y;
	block.col1[1] =  XA1_YA1*nvec.y*nvec.z;
	block.col1[2] = -YB[2]*nvec.x;
	// column 2
	block.col2[0] =  YA[1] + XA1_YA1*nvec.z*nvec.z;
	// column 3
	block.col3[0] =  YB[1]*nvec.z;
	block.col3[1] = -YB[1]*nvec.y;
	block.col3[2] =  YC[1]*(1-nvec.x*nvec.x);

	block.col3[3] = -YC[1]*nvec.x*nvec.y;
	block.col3[4] = -YC[1]*nvec.x*nvec.z;
	// column 4
	block.col4[0] =  YB[1]*nvec.x;
	block.col4[1] =  YC[1]*(1-nvec.y*nvec.y);
	block.col4[2] = -YC[1]*nvec.y*nvec.z;
	// column 5
	block.col5[0] =  YC[1]*(1-nvec.z*nvec.z);
	return block;
}


struct ODBlock Lubrication::RFU_ODBlock() const
{
	if (tangential) {
		return RFU_ODBlock_squeeze_tangential();
	} else {
		return RFU_ODBlock_squeeze();
	}
}


// Diagonal Blocks Terms, F/U version
std::pair<struct DBlock, struct DBlock> Lubrication::RFU_DBlocks_squeeze() const
{
	struct DBlock b0;
	struct DBlock b1;
	const auto &nvec = inter->nvec;

	// (*,0)
	b0.col0[0] = XA[0]*nvec.x*nvec.x; // 00 element of the dblock
	b0.col0[1] = XA[0]*nvec.x*nvec.y;           // 10
	b0.col0[2] = XA[0]*nvec.x*nvec.z;           // 20
	b0.col0[3] = 0;                   // 40
	b0.col0[4] = 0;                   // 50
	// (*,1)
	b0.col1[0] = XA[0]*nvec.y*nvec.y; // 11
	b0.col1[1] = XA[0]*nvec.y*nvec.z; // 21
	b0.col1[2] = 0;                   // 51
	// (*,2)
	b0.col2[0] = XA[0]*nvec.z*nvec.z ; // 22
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
 	b1.col0[0] = XA[3]*nvec.x*nvec.x; // 00 element of the dblock
 	b1.col0[1] = XA[3]*nvec.x*nvec.y;           // 10
 	b1.col0[2] = XA[3]*nvec.x*nvec.z;           // 20
 	b1.col0[3] = 0;                   // 40
 	b1.col0[4] = 0;                   // 50
 	// (*,1)
 	b1.col1[0] = XA[3]*nvec.y*nvec.y; // 11
	b1.col1[1] = XA[3]*nvec.y*nvec.z; // 21
 	b1.col1[2] = 0;                   // 51
 	// (*,2)
	b1.col2[0] = XA[3]*nvec.z*nvec.z; // 22
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
	const auto &nvec = inter->nvec;

	double n0n0 = nvec.x*nvec.x;
	double n0n1 = nvec.x*nvec.y;
	double n0n2 = nvec.x*nvec.z;
	double n1n1 = nvec.y*nvec.y;
	double n1n2 = nvec.y*nvec.z;
	double n2n2 = nvec.z*nvec.z;
	double one_n0n0 = 1-n0n0;
	double one_n1n1 = 1-n1n1;
	double one_n2n2 = 1-n2n2;
	double XA0mYA0 = XA[0]-YA[0];
	// (*,0)
	b0.col0[0] =  XA[0]*n0n0+YA[0]*one_n0n0; // 00 element of the dblock
	b0.col0[1] = XA0mYA0*n0n1;           // 10
	b0.col0[2] = XA0mYA0*n0n2;           // 20
	b0.col0[3] = -YB[0]*nvec.z;                   // 40
	b0.col0[4] =  YB[0]*nvec.y;                   // 50
	// (*,1)
	b0.col1[0] =  XA[0]*n1n1+YA[0]*one_n1n1; // 11
	b0.col1[1] = XA0mYA0*n1n2;           // 21
	b0.col1[2] = -YB[0]*nvec.x;                   // 51
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
 	b1.col0[3] = -YB[3]*nvec.z;                   // 40
 	b1.col0[4] =  YB[3]*nvec.y;                   // 50
 	// (*,1)
 	b1.col1[0] =  XA[3]*n1n1+YA[3]*one_n1n1; // 11
 	b1.col1[1] =  XA3mYA3*n1n2;           // 21
 	b1.col1[2] = -YB[3]*nvec.x;                   // 51
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

std::pair<struct DBlock, struct DBlock> Lubrication::RFU_DBlocks() const
{
	if (tangential) {
		return RFU_DBlocks_squeeze_tangential();
	} else {
		return RFU_DBlocks_squeeze();
	}
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
	const auto &nvec = inter->nvec;
	double nvec_vi = dot(nvec, vi);
	double nvec_vj = dot(nvec, vj);
	Sym2Tensor nvec_nvec = outer(nvec);
	stresslet_i += -(XG[0]*nvec_vi+XG[1]*nvec_vj)*nvec_nvec; // XGU_i
	stresslet_j += -(XG[2]*nvec_vi+XG[3]*nvec_vj)*nvec_nvec; // XGU_j
	if (!tangential) {
		return;
	}
	Sym2Tensor tmp_i = outer_sym(nvec, vi) - nvec_vi*nvec_nvec;
	Sym2Tensor tmp_j = outer_sym(nvec, vj) - nvec_vj*nvec_nvec;
	stresslet_i += -2*(YG[0]*tmp_i+YG[1]*tmp_j); // YGU_i
	stresslet_j += -2*(YG[2]*tmp_i+YG[3]*tmp_j); // YGU_j
	Sym2Tensor tmp2_i = outer_sym(nvec, cross(oi, nvec));
	Sym2Tensor tmp2_j = outer_sym(nvec, cross(oj, nvec));
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
	const auto &nvec = inter->nvec;
	double nnE = dot(nvec, dot(E_inf, nvec));
	Sym2Tensor nvec_nvec = outer(nvec);
	double coeff = 1.5*nnE;
	stresslet_i += coeff*(XM[0]+XM[1])*nvec_nvec; // XME_i
	stresslet_j += coeff*(XM[2]+XM[3])*nvec_nvec; // XME_j
	if (tangential) {
		vec3d Einf_nvec = dot(E_inf, nvec);
		Sym2Tensor nvec_Einf_nvec = outer_sym(nvec, Einf_nvec);
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
vec3d Lubrication::getTotalForce(const vec3d &na_v0, 
								 const vec3d &na_v1, 
								 const vec3d &na_ang_v0, 
								 const vec3d &na_ang_v1,
								 const Sym2Tensor &E_inf) const
{
	/**
	* \brief The total lubrication force acting on p0.
	*/

	/* Eq. (1.6a) in Jeffrey&Onishi 1984
	 * A_{ij}^{ab} = XA_{ab}ni*nj + YA_{ab}(del_{ij}-ni*nj)
	 * B~_{ji}^{ab} = YB_{ba}epsilon_{jik} nk
	 *
	 */
	const auto &nvec = inter->nvec;

	/* XAU_i */
	vec3d lubforce_p0 = -dot(XA[0]*na_v0+XA[1]*na_v1, nvec)*nvec;
	if (tangential) {
		/* YAU_i */
		lubforce_p0 += -YA[0]*(na_v0-nvec*dot(nvec, na_v0)) - YA[1]*(na_v1-nvec*dot(nvec, na_v1));
		/* YBO_i */
		lubforce_p0 += -YB[0]*cross(nvec, na_ang_v0)        - YB[2]*cross(nvec, na_ang_v1);
	}
	vec3d GEi, GEj;
	if (!tangential) {
		std::tie(GEi, GEj) = calcGE_squeeze(E_inf);
	} else {
		std::tie(GEi, GEj) = calcGE_squeeze_tangential(E_inf);
	}
	lubforce_p0 += GEi;
	return lubforce_p0;
}

void Lubrication::updateResistanceCoeff()
{
	if (inter->getReducedGap() > 0) {
		double coeff = 1/(inter->getReducedGap()+regularization_length);
		double coeff_norm = coeff-lub_coeff_min;
		if (coeff > 1.) {
			setResistanceCoeff(coeff_norm, log(coeff));
		} else {
			setResistanceCoeff(coeff_norm, 0.);
		}
	}
}

} // namespace Lub

} // namespace Interactions
