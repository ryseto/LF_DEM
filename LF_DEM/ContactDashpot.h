//
//  ContactDashpot.h
//  LF_DEM
//
//  Copyright (c) 2016 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class ContactDashpot
 \brief ContactDashpot interaction, to be called from a Contact object

	Contains both normal and tangential dashpots.
	As for the force, it formally acts as a distance independant Lubrication object.
	The stresses, however, have to be computed as -xF, from the Contact object.

 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__ContactDashpot__
#define __LF_DEM__ContactDashpot__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "MatrixBlocks.h"
#include "PairwiseInteraction.h"
#include "PairVelocity.h"


namespace Interactions
{

class ContactDashpot{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	PairwiseInteraction *interaction;
	//======= particles data  ====================//
	double ro;
	double ro_12; // = ro/2
	double normal_coeff;
	double tangential_coeff;
	double XA[4]; // ii ij ji jj
	double YA[4]; // ii ij ji jj
	double YB[4]; // ii ij ji jj
	double YC[4]; // ii ij ji jj
	void calcDashpotResistances();
	void setDashpotResistanceCoeffs(double normal, double tangential);
 public:
	ContactDashpot(PairwiseInteraction *inter, double norm_coeff, double tan_coeff);

	//===== forces/stresses  ==========================
	std::tuple<vec3d, vec3d, vec3d, vec3d> getForcesTorques(const struct PairVelocity &pvel) const;	
	vec3d getForceOnP0(const struct PairVelocity &pvel) const;
	
	//=============  Resistance Matrices ====================/
	struct ODBlock RFU_ODBlock() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks() const;
};

} // namespace Interactions
#endif /* defined(__LF_DEM__ContactDashpot__) */
