//
//  Lubrication.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2016 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Lubrication
 \brief Lubrication interaction, to be called from an Interaction object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__Lubrication__
#define __LF_DEM__Lubrication__
#include <string>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "MatrixBlocks.h"
#include "LubricationCoefficients.h"
#include "LubricationParams.h"
#include "PairVelocity.h"


namespace Interactions {

class PairwiseInteraction;

namespace Lub {

class Lubrication{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	PairwiseInteraction *inter;
	double regularization_length;
	double lub_coeff;
	double log_lub_coeff;
	double lub_coeff_min;
	double XA[4]; // ii ij ji jj
	double YA[4]; // ii ij ji jj
	double YB[4]; // ii ij ji jj
	double YC[4]; // ii ij ji jj
	double XG[4]; // ii ij ji jj
	double YG[4]; // ii ij ji jj
	double YH[4]; // ii ij ji jj
	double XM[4]; // ii ij ji jj
	double YM[4]; // ii ij ji jj

	LubricationCoefficients coeffs;

	void setResistanceCoeff(double normal_rc, double tangent_rc);

	//=============  Resistance Matrices ====================/
	void calcXFunctionsStress();
	void calcXYFunctionsStress();
	void calcXFunctions();
	void calcXYFunctions();

	struct ODBlock RFU_ODBlock_squeeze_tangential() const;
	struct ODBlock RFU_ODBlock_squeeze() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks_squeeze_tangential() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks_squeeze() const;


	std::tuple<vec3d,vec3d> calcGE_squeeze(const Sym2Tensor& E_inf) const;
	std::tuple<vec3d,vec3d> calcGE_squeeze(const Sym2Tensor& E0, const Sym2Tensor& E1) const;
	std::tuple<vec3d,vec3d> calcGE_squeeze_tangential(const Sym2Tensor& E_inf) const;
	std::tuple<vec3d,vec3d> calcGE_squeeze_tangential(const Sym2Tensor& E0, const Sym2Tensor& E1) const;
	
	bool tangential;

 public:
	Lubrication(PairwiseInteraction* interaction, const LubParams &params);
	void updateResistanceCoeff();

	//===== forces/stresses  ========================== //
	vec3d getTotalForce(const struct PairVelocity &vel, 
						const Sym2Tensor &E_inf) const;
	void addMEStresslet(const Sym2Tensor& E_inf,
						Sym2Tensor& stresslet_i,
						Sym2Tensor& stresslet_j) const;
	void addGUStresslet(const struct PairVelocity &vel, 
						Sym2Tensor& stresslet_i, Sym2Tensor& stresslet_j) const;

	std::tuple<vec3d,vec3d> calcGE(const Sym2Tensor& E_inf) const;
	std::tuple<vec3d,vec3d> calcGE(const Sym2Tensor& E0, const Sym2Tensor& E1) const;
	std::tuple<vec3d,vec3d,vec3d,vec3d> calcGEHE(const Sym2Tensor& E_inf) const;
	std::tuple<vec3d,vec3d,vec3d,vec3d> calcGEHE(const Sym2Tensor& E0, const Sym2Tensor& E1) const;

	struct ODBlock RFU_ODBlock() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks() const;
};

inline double calcLubricationRange(double input_max_gap, double a0, double a1)
{
	double max_gap;
	double rad_ratio = a0/a1;
	if (rad_ratio > 2 || rad_ratio < 0.5) {
		double minradius = (a0<a1 ? a0 : a1);
		max_gap = input_max_gap*2*minradius/(a0+a1);   // max gap proportional to min radius when large size diff
	} else {
		max_gap = input_max_gap;
	}
	return (a0 + a1)*(1 + 0.5*max_gap);
}

} // namespace Lub

} // namespace Interactions

#endif /* defined(__LF_DEM__Lubrication__) */
