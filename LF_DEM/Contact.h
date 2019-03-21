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
#include <iomanip>
#include <fstream>
#include <vector>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "ContactDashpot.h"

class System;
class Interaction;

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
	System *sys;
	Interaction *interaction;

	void (Contact::*frictionlaw)();
	bool active;
	unsigned int p0;
	unsigned int p1;
	double a0; // radii
	double a1;
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
	double ft_max; // friction_model = 5;
	vec3d rolling_velocity;
	int state;

	void incrementTangentialDisplacement();
	void incrementRollingDisplacement();

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
	inline void setTangentialForceNorm(double, double);
	inline void setRollingForceNorm(double, double);

public:
	/*********************************
	 *       Public Methods          *
	 *********************************/
	//======= internal state =====================//
	Contact():
	active(false),
	mu_static(0),
	mu_dynamic(0),
	mu_rolling(0),
	f_spring_normal_norm(0),
	normal_load(0),
	f_spring_normal(0),
	f_spring_tan(0),
	f_spring_total(0),
	f_rolling(0),
	ft_max(0)
	{};

	void init(System* sys_, Interaction* int_);
	ContactDashpot dashpot;
	void setInteractionData();
	void setSpringConstants();
	void setDashpotConstants();
	void activate();
	void deactivate();
	vec3d disp_tan; // tangential displacement
	vec3d disp_rolling;
	vec3d prev_disp_tan; // useful for predictor-corrector method: disp_tan in the previous time step
	vec3d prev_disp_rolling;
	void incrementDisplacements();
	double get_rcontact() const
	{
		return a0 + a1;
	}
	vec3d getSlidingVelocity() const;
	vec3d getRollingVelocity() const;

	//===== forces/stresses  ========================== //
	void addUpForceTorque(std::vector<vec3d> &force_per_particle,
						  std::vector<vec3d> &torque_per_particle) const;
	void addUpForce(std::vector<vec3d> &force_per_particle) const;
	void calcContactSpringForce();
	vec3d getTotalForce() const;
	vec3d getSpringForce() const;
	vec3d getNormalForce() const;
	double getNormalForceValue() const;
	double getNormalSpringForce() const;
	vec3d getTangentialForce() const;
	double get_normal_load() const;
	void calcContactStress();
	void addUpStress(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1);
	void addUpStressSpring(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1) const;
	Sym2Tensor getContactStressXF() const
	{
		return contact_stresslet_XF;
	}
	double calcEnergy() const;
	struct contact_state getState() const
	{
		struct contact_state cs;
		cs.p0 = p0;
		cs.p1 = p1;
		cs.disp_tan = disp_tan;
		cs.disp_rolling = disp_rolling;
		return cs;
	};
	void setState(const struct contact_state& cs)
	{
		p0 = cs.p0;
		p1 = cs.p1;
		disp_tan = cs.disp_tan;
		disp_rolling = cs.disp_rolling;
	}
	bool is_active() const
	{
		return active;
	}
	bool is_frictional() const
	{
		return state >= 2;
	}
	int getFrictionState() const
	{
		return state;
	}
};
#endif /* defined(__LF_DEM__Contact__) */
