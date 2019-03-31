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
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "MatrixBlocks.h"

class System;
class Interaction;

class Lubrication{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	Interaction *interaction;
	System *sys;
	bool _active;
	double range;
	//======= particles data  ====================//
	unsigned int p0;
	unsigned int p1;
	unsigned int p0_6;
	unsigned int p1_6;
	vec3d *nvec;
	double lub_coeff;
	double log_lub_coeff;
	double lub_coeff_min;
	double a0;
	double a1;
	double ro;
	double ro_12; // = ro/2
	double a0a0_23;
	double a1a1_23;
	double roro_16;
	double a0a0a0_43;
	double a1a1a1_43;
	double rororo_16;
	double a0a0a0_109;
	double a1a1a1_109;
	double rororo_536;
	double lambda;
	double invlambda;
	double lambda_square;
	double lambda_cubic;
	double lambda_p_1;
	double lambda_p_1_square;
	double lambda_p_1_cubic;
	double XA[4]; // ii ij ji jj
	double YA[4]; // ii ij ji jj
	double YB[4]; // ii ij ji jj
	double YC[4]; // ii ij ji jj
	double XG[4]; // ii ij ji jj
	double YG[4]; // ii ij ji jj
	double YH[4]; // ii ij ji jj
	double XM[4]; // ii ij ji jj
	double YM[4]; // ii ij ji jj
	double cXA[4];
	double cYA[4];
	double cYB[4];
	double cYC[4];
	double cXG[4];
	double cYG[4];
	double cYH[4];
	double cXM[4];
	double cYM[4];
	double g1_XA;
	double g1_inv_XA;
	double g2_YA;
	double g2_inv_YA;
	double g2_YB;
	double g2_inv_YB;
	double g2_YC;
	double g2_inv_YC;
	double g4_YC;
	double g1_XG;
	double g1_inv_XG;
	double g2_YG;
	double g2_inv_YG;
	double g2_YH;
	double g2_inv_YH;
	double g5_YH;
	double g5_inv_YH;
	double g1_XM;
	double g1_inv_XM;
	double g4_XM;
	double g2_YM;
	double g2_inv_YM;
	double g5_YM;

 public:
	Lubrication();
	bool tangential; // 1: normal; 2: normal+tangential
	void init(System *sys_, Interaction *int_);
	void setParent(Interaction *int_){interaction = int_;};
	inline bool is_active() const {return _active;};
	void activate();
	void deactivate();
	void updateActivationState(bool contact_active);
	void setParticleData();
	void calcLubConstants();
	//===== forces/stresses  ========================== //
	vec3d getTangentialForce() const
	{
		vec3d lubforce_p0 = getTotalForce();
		return lubforce_p0-dot(lubforce_p0, nvec)*(*nvec);
	}
	vec3d getTotalForce() const;
	double getNormalForceValue() const
	{
		return dot(getTotalForce(), nvec);
	}
	void addMEStresslet(const Sym2Tensor& E_inf,
						Sym2Tensor& stresslet_i,
						Sym2Tensor& stresslet_j) const;
	void addGUStresslet(const vec3d& vi, const vec3d& vj,
						const vec3d& oi, const vec3d& oj,
						Sym2Tensor& stresslet_i, Sym2Tensor& stresslet_j) const;
	void updateResistanceCoeff();
	void setResistanceCoeff(double normal_rc, double tangent_rc);
	//void setResistanceCoeffTang(double tangent_rc);
	//=============  Resistance Matrices ====================/
	void calcXFunctionsStress();
	void calcXYFunctionsStress();
	std::tuple<vec3d,vec3d> calcGE_squeeze(const Sym2Tensor& E_inf) const;
	std::tuple<vec3d, vec3d> calcGE_squeeze(const Sym2Tensor& E_local_0,
											const Sym2Tensor& E_local_1) const;
	std::tuple<vec3d,vec3d> calcGE_squeeze_tangential(const Sym2Tensor& E_inf) const;
	std::tuple<vec3d,vec3d> calcGE_squeeze_tangential(const Sym2Tensor& E_local_0,
													  const Sym2Tensor& E_local_1) const;
	
	std::tuple<vec3d,vec3d,vec3d,vec3d> calcGEHE_squeeze_tangential(const Sym2Tensor& E_inf) const;

	std::tuple<vec3d,vec3d,vec3d,vec3d> calcGEHE_squeeze_tangential(const Sym2Tensor& E_local_0,
																	const Sym2Tensor& E_local_1) const;
	
	
	struct ODBlock RFU_ODBlock_squeeze_tangential() const;
	struct ODBlock RFU_ODBlock_squeeze() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks_squeeze_tangential() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks_squeeze() const;
	void calcXFunctions();
	void calcXYFunctions();
	void calcLubricationForce();
	void resetLubricationForce()
	{
		force = 0;
	}
	double force;
	
};

#endif /* defined(__LF_DEM__Lubrication__) */
