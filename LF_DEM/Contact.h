//
//  Contact.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Contact
 \brief Contact object, to be called from an Interaction object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__Contact__
#define __LF_DEM__Contact__
#include <iostream>
#include <vector>
#include <memory>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "ContactDashpot.h"
#include "ContactParams.h"
#include "PairwiseInteraction.h"
#include "PairVelocity.h"

namespace Interactions
{

struct contact_state {
	unsigned int p0;
	unsigned int p1;
	vec3d disp_tan;
	vec3d disp_rolling;
};

namespace Contact_ios {
	std::vector <struct contact_state> readStatesBStream(std::istream &input, unsigned int np);
	void writeStatesBStream(std::ostream &conf_export, const std::vector <struct contact_state> &cs);
}

class Contact {
private:
	/*********************************
	 *        Members                *
	 *********************************/
	PairwiseInteraction *interaction;
	FrictionModel friction_model;
	void (Contact::*frictionlaw)();

	double a_reduced;
	//===== forces and stresses ==================== //
	double kt_scaled;
	double kr_scaled;
	double kn_scaled;
	double mu_static;
	double mu_dynamic;
	double mu_rolling;
	/*********************************
	 *       Private Methods         *
	 *********************************/
	//===== forces and stresses computations =====//
	Sym2Tensor contact_stresslet_XF; //stress tensor of contact force
	double f_spring_normal_norm; // normal contact force
	double normal_load; // compressive load + cohesion. If it is positive, particies are in cohesive contact state.
	vec3d f_spring_normal; // normal contact force, spring only
	vec3d f_spring_tan; // tangential contact force, spring only
	vec3d f_spring_total; // spring only
	vec3d f_rolling;
	double ft_max; // friction_model = ft_max;
	double adhesion;
	double critical_load;

	vec3d rolling_velocity;
	int state;

	void incrementTangentialDisplacement(double dt, const struct PairVelocity &vel);
	void incrementRollingDisplacement(double dt, const struct PairVelocity &vel);

	/* state:
	 * 0 No contact
	 * 1 Friction is not activated (critical load model)
	 * 2 Static friction
	 * 3 Sliding
	 * 4 Infinite friction coeffient
	 * -2 Switching dynamic to static
	 */
	void frictionlaw_criticalload();
	void frictionlaw_criticalload_mu_inf();
	void frictionlaw_standard();
	void frictionlaw_ft_max();
	void frictionlaw_coulomb_max();
	void frictionlaw_infinity();
	void setTangentialForceNorm(double, double);
	void setRollingForceNorm(double, double);

	void calcContactStress(const struct PairVelocity &vel);

	vec3d prev_disp_tan; // useful for predictor-corrector method: disp_tan in the previous time step
	vec3d prev_disp_rolling;

	bool rolling_friction() const {return mu_rolling != 0;}

	void setInteractionData(const ContactParams &p);
	void setSpringConstants(const ContactParams &p);
	void setDashpotConstants(const ContactParams &p);

public:
	/*********************************
	 *       Public Methods          *
	 *********************************/
	//======= internal state =====================//
	Contact(PairwiseInteraction* interaction_, 
			const ContactParams &params,
			double norm_dashpot_coeff, 
			double tan_dashpot_coeff);

	std::unique_ptr<ContactDashpot> dashpot;
	
	vec3d disp_tan; // tangential displacement
	vec3d disp_rolling;
	
	void incrementDisplacements(double dt, const struct PairVelocity &vel);

	//===== forces/stresses  ========================== //
	void addUpForceTorque(std::vector<vec3d> &force_per_particle,
						  std::vector<vec3d> &torque_per_particle) const;
	void addUpForce(std::vector<vec3d> &force_per_particle) const;
	void calcContactSpringForce();
	vec3d getTotalForce(const struct PairVelocity &vel) const;
	vec3d getSpringForce() const;
	double getNormalForceValue(const struct PairVelocity &vel) const;
	double getNormalSpringForce() const;
	vec3d getTangentialForce(const struct PairVelocity &vel) const;
	vec3d getSlidingVelocity(const struct PairVelocity &vel) const;
	vec3d getRollingVelocity(const struct PairVelocity &vel) const;
	double get_normal_load() const;
	void addUpStress(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1,
					 const struct PairVelocity &vel);
	void addUpStressSpring(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1) const;
	void addUpStressDashpot(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1, 
						    const PairVelocity &pvel) const;
	Sym2Tensor getContactStressXF() const
	{
		return contact_stresslet_XF;
	}
	double calcEnergy() const;
	struct contact_state getState() const;
	void setState(const struct contact_state& cs);
	void saveState();
	void restoreState();

	bool is_frictional() const
	{
		return state >= 2;
	}
	int getFrictionState() const
	{
		return state;
	}
};

inline double calcContactRange(double a0, double a1)
{
	return a0 + a1;
}

} // namespace Interactions

#endif /* defined(__LF_DEM__Contact__) */
