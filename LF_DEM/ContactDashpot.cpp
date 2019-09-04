//
//  Lubrication.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <tuple>
#include "ContactDashpot.h"
#include "LubricationCoefficients.h"

namespace Interactions {

ContactDashpot::ContactDashpot(PairwiseInteraction *inter, double norm_coeff, double tan_coeff) :
interaction(inter)
{
	for (int i=0; i<4; i++) {
		XA[i] = 0;
		YA[i] = 0;
		YB[i] = 0;
		YC[i] = 0;
	}

	ro = inter->a0+inter->a1;
	ro_12 = ro/2;
	setDashpotResistanceCoeffs(norm_coeff, tan_coeff);
}

void ContactDashpot::setDashpotResistanceCoeffs(double normal, double tangential)
{
	normal_coeff = normal;
	tangential_coeff = tangential;
	calcDashpotResistances();
}

void ContactDashpot::calcDashpotResistances()
{
	Lub::LubricationCoefficients coeffs (interaction->a0, interaction->a1);
	for (unsigned j=0; j<4; j++) {
		XA[j] = normal_coeff*coeffs.cXA[j];
	}

	if (tangential_coeff == 0) {
		return;
	}

	for (unsigned j=0; j<4; j++) {
		YA[j] = tangential_coeff*coeffs.cYA[j];
	}
		for (unsigned j=0; j<4; j++) {
		YB[j] = tangential_coeff*coeffs.cYB[j];
	}
	for (unsigned j=0; j<4; j++) {
		YC[j] = tangential_coeff*coeffs.cYC[j];
	}
}

// for FT/UW version
struct ODBlock ContactDashpot::RFU_ODBlock() const
{
	struct ODBlock block;
	const auto &nvec = interaction->nvec;

	if (tangential_coeff > 0) {
		double XA1_YA1 = XA[1]-YA[1];
		// column 0
		block.col0[0] = YA[1] + XA1_YA1*nvec.x*nvec.x;
		block.col0[1] = XA1_YA1*nvec.x*nvec.y;
		block.col0[2] = XA1_YA1*nvec.x*nvec.z;
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
	} else {
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
	}
	return block;
}

// Diagonal Blocks Terms, FT/UW version
std::pair<struct DBlock, struct DBlock> ContactDashpot::RFU_DBlocks() const
{
	struct DBlock b0;
	struct DBlock b1;
	const auto &nvec = interaction->nvec;

	if (tangential_coeff > 0) {
		double XA0_YA0 = XA[0] - YA[0];
		// (*,0)
		b0.col0[0] = YA[0] + XA0_YA0*nvec.x*nvec.x; // 00 element of the dblock
		b0.col0[1] = XA0_YA0*nvec.x*nvec.y;         // 10
		b0.col0[2] = XA0_YA0*nvec.x*nvec.z;           // 20
		b0.col0[3] = -YB[0]*nvec.z;                   // 40
		b0.col0[4] =  YB[0]*nvec.y;                   // 50
		// (*,1)
		b0.col1[0] = YA[0] + XA0_YA0*nvec.y*nvec.y; // 11
		b0.col1[1] = XA0_YA0*nvec.y*nvec.z;           // 21
		b0.col1[2] = -YB[0]*nvec.x;                   // 51
		// (*,2)
		b0.col2[0] = YA[0] + XA0_YA0*nvec.z*nvec.z; // 22
		// (*,3)
		b0.col3[0] =  YC[0]*(1-nvec.x*nvec.x);                 // 33
		b0.col3[1] = -YC[0]*nvec.x*nvec.y;                     // 43
		b0.col3[2] = -YC[0]*nvec.x*nvec.z;                     // 53
		// (*,4)
		b0.col4[0] =  YC[0]*(1-nvec.y*nvec.y);                 // 44
		b0.col4[1] = -YC[0]*nvec.y*nvec.z;                    // 54
		// (*,5)
		b0.col5[0] =  YC[0]*(1-nvec.z*nvec.z);                 // 55
		
		double XA3_YA3 = XA[3] - YA[3];
		// (*,0)
		b1.col0[0] = YA[3] + XA3_YA3*nvec.x*nvec.x; // 00 element of the dblock
		b1.col0[1] = XA3_YA3*nvec.x*nvec.y;           // 10
		b1.col0[2] = XA3_YA3*nvec.x*nvec.z;          // 20
		b1.col0[3] = -YB[3]*nvec.z;                   // 40
		b1.col0[4] =  YB[3]*nvec.y;                   // 50
		// (*,1)
		b1.col1[0] = YA[3] + XA3_YA3*nvec.y*nvec.y; // 11
		b1.col1[1] = XA3_YA3*nvec.y*nvec.z;           // 21
		b1.col1[2] = -YB[3]*nvec.x;                   // 51
		// (*,2)
		b1.col2[0] = YA[3] + XA3_YA3*nvec.z*nvec.z; // 22
		// (*,3)
		b1.col3[0] =  YC[3]*(1-nvec.x*nvec.x);                 // 33
		b1.col3[1] = -YC[3]*nvec.x*nvec.y;                     // 43
		b1.col3[2] = -YC[3]*nvec.x*nvec.z;                     // 53
		// (*,4)
		b1.col4[0] =  YC[3]*(1-nvec.y*nvec.y);                 // 44
		b1.col4[1] = -YC[3]*nvec.y*nvec.z;                     // 54
		// (*,5)
		b1.col5[0] =  YC[3]*(1-nvec.z*nvec.z);;                 // 55
	} else {
		b0.col0[0] = XA[0]*nvec.x*nvec.x; // 00 element of the dblock
		b0.col0[1] = XA[0]*nvec.x*nvec.y;           // 10
		b0.col0[2] = XA[0]*nvec.x*nvec.z;           // 20
		b0.col0[3] = 0;                   // 40
		b0.col0[4] = 0;                   // 50
		// (*,1)
		b0.col1[0] = XA[0]*nvec.y*nvec.y; // 11
		b0.col1[1] = XA[0]*nvec.y*nvec.z;           // 21
		b0.col1[2] = 0;                   // 51
		// (*,2)
		b0.col2[0] = XA[0]*nvec.z*nvec.z; // 22
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
		b1.col1[1] = XA[3]*nvec.y*nvec.z;           // 21
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
	}
	return std::make_pair(b0, b1);
}


vec3d ContactDashpot::getForceOnP0(const struct PairVelocity &pvel) const
{
	/** \brief Resistance force acting on particle p0.

	 Used by the dynamics to get the R_FU*U_inf force term and by the output data.

	 */

	const auto &nvec = interaction->nvec;
	/* XAU_i */
	vec3d force_p0 = -dot(XA[0]*pvel.U[0]+XA[1]*pvel.U[1], nvec)*nvec;
	if (tangential_coeff > 0) {
		/* YAU_i */
		force_p0 += -YA[0]*(pvel.U[0]-nvec*dot(nvec, pvel.U[0])) \
					-YA[1]*(pvel.U[1]-nvec*dot(nvec, pvel.U[1]));
		/* YBO_i */
		force_p0 += -YB[0]*cross(nvec, pvel.O[0]) \
					-YB[2]*cross(nvec, pvel.O[1]);
	}
	return force_p0;
}

std::tuple<vec3d, vec3d, vec3d, vec3d> ContactDashpot::getForcesTorques(const struct PairVelocity &pvel) const
{
	/** \brief */
	vec3d force_p0 = getForceOnP0(pvel);
	vec3d torque_p0;
	if (tangential_coeff > 0) {
		torque_p0 = interaction->a0*cross(interaction->nvec, force_p0);
	}
	return std::make_tuple(force_p0, -force_p0, torque_p0, (interaction->a1/interaction->a0)*torque_p0);
}

} // namespace Interactions