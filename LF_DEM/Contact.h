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
#include "vec3d.h"
#include "StressTensor.h"
#include "ContactDashpot.h"

class System;
class Interaction;

struct contact_state {
	unsigned int p0;
	unsigned int p1;
	vec3d disp_tan;
	vec3d disp_rolling;
};

class Contact{
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
	StressTensor contact_stresslet_XF; //stress tensor of contact force
	double f_spring_normal_norm; // normal contact force
	double normal_load; // compressive load + cohesion. If it is positive, particies are in cohesive contact state.
	vec3d f_spring_normal; // normal contact force, spring only
	vec3d f_spring_tan; // tangential contact force, spring only
	vec3d f_spring_total; // spring only
	vec3d f_rolling;
	double ft_max; // friction_model = 5;
	vec3d rolling_velocity;
	int state;

	void calcRollingVelocities();
	void incrementTangentialDisplacement();
	void incrementRollingDisplacement();

	/* state:
	 * 0 No contact
	 * 1 Friction is not activated (critical load model)
	 * 2 Static friction
	 * 3 Sliding
	 * -2 Switching dynamic to static
	 */
	void frictionlaw_criticalload();
	void frictionlaw_criticalload_mu_inf();
	void frictionlaw_standard();
	void frictionlaw_ft_max();
	void frictionlaw_coulomb_max();

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
	ft_max(0),
	rolling_velocity(0),
	relative_surface_velocity_sqnorm(0)
	{};

	void init(System* sys_, Interaction* int_);
	ContactDashpot dashpot;
	void setInteractionData();
	void setSpringConstants();
	void activate();
	void deactivate();
	vec3d disp_tan; // tangential displacement
	vec3d disp_rolling;
	vec3d prev_disp_tan; // useful for predictor-corrector method: disp_tan in the previous time step
	vec3d prev_disp_rolling;
	double relative_surface_velocity_sqnorm;
	void incrementDisplacements();
	double get_rcontact()
	{
			return a0 + a1;
	}
	//===== forces/stresses  ========================== //
	void calcContactSpringForce();
	void addUpContactForceTorque();
	vec3d getTotalForce();
	vec3d getNormalForce();
	vec3d getTangentialForce();
	double get_normal_load();
	void calcContactStress();
	StressTensor getContactStressXF()
	{
		return contact_stresslet_XF;
	}
	double calcEnergy();
	struct contact_state getState() {
		struct contact_state cs;
		cs.p0 = p0;
		cs.p1 = p1;
		cs.disp_tan = disp_tan;
		cs.disp_rolling = disp_rolling;
		return cs;
	};
	void setState(const struct contact_state& cs) {
		p0 = cs.p0;
		p1 = cs.p1;
		disp_tan = cs.disp_tan;
		disp_rolling = cs.disp_rolling;
	}
	inline bool is_active()
	{
		return active;
	}
	inline bool is_frictional()
	{
		return state >= 2;
	}
	int getFrictionState()
	{
		return state;
	}
};
#endif /* defined(__LF_DEM__Contact__) */
