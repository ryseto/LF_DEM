//
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Contact.h"

namespace Interactions 
{

Contact::Contact(PairwiseInteraction* interaction_, 
				 const ContactParams &p, 
				 double norm_dashpot_coeff, 
				 double tan_dashpot_coeff) :
interaction(interaction_),
friction_model(p.friction_model)
{
	if (norm_dashpot_coeff != 0 || tan_dashpot_coeff != 0) {
		dashpot = std::unique_ptr<ContactDashpot>(new ContactDashpot (interaction_, norm_dashpot_coeff, tan_dashpot_coeff));
	}

	switch (friction_model) {
		case FrictionModel::frictionless:
			state = 1;
			break;
		case FrictionModel::Coulomb:
			frictionlaw = &Contact::frictionlaw_standard;
			state = 2;
			break;
		case FrictionModel::criticalload:
			frictionlaw = &Contact::frictionlaw_criticalload;
			critical_load = p.critical_load;
			state = 1; // critical load model
			break;
		case FrictionModel::criticalload_mu_inf:
			frictionlaw = &Contact::frictionlaw_criticalload_mu_inf;
			critical_load = p.critical_load;
			state = 1; // critical load model			
			break;
		case FrictionModel::infinity:
			frictionlaw = &Contact::frictionlaw_infinity;
			state = 2;
			break;
		case FrictionModel::ft_max:
			frictionlaw = &Contact::frictionlaw_ft_max;
			ft_max = p.ft_max;
			state = 2;
			break;
		case FrictionModel::Coulomb_max:
			frictionlaw = &Contact::frictionlaw_coulomb_max;
			ft_max = p.ft_max;
			state = 2;
			break;
	}

	setInteractionData(p);

	disp_tan.reset();
	disp_rolling.reset();
	prev_disp_tan.reset();
	prev_disp_rolling.reset();
	adhesion = p.adhesion;

	calcContactSpringForce();
}

void Contact::setSpringConstants(const ContactParams &p)
{
	double ro_12 = (interaction->a0+interaction->a1)/2;
	kn_scaled = ro_12*ro_12*p.kn; // F = kn_scaled * _reduced_gap;  <-- gap is scaled @@@@ Why use reduced_gap? Why not gap?
	kt_scaled = ro_12*p.kt; // F = kt_scaled * disp_tan <-- disp is not scaled
	kr_scaled = ro_12*p.kr; // F = kt_scaled * disp_tan <-- disp is not scaled
}

void Contact::setInteractionData(const ContactParams &p)
{
	/* [note]
	 * The reduced (or effective) radius is defined as
	 * 1/a_reduced = 1/a0 + 1/a1
	 * This definition comes from sphere vs half-plane geometory of contact mechanics.
	 * If sphere contacts with a half-plane (a1 = infinity), a_reduced = a0.
	 * For equal sized spheres, a_reduced = 0.5*a0 = 0.5*a1
	 */
	a_reduced = interaction->a0*interaction->a1/(interaction->a0+interaction->a1);
	setSpringConstants(p);
	mu_static = p.mu_static;
	mu_dynamic = p.mu_dynamic;
	mu_rolling = p.mu_rolling;
}

/*********************************
 *                                *
 *	   Contact Forces Methods     *
 *                                *
 *********************************/
vec3d Contact::getSlidingVelocity(const struct PairVelocity &vel) const
{
	vec3d rotational_deltav = -cross(interaction->a0*vel.O[0]+interaction->a1*vel.O[1], 
									 interaction->nvec);

	vec3d sliding_velocity = vel.U[1]-vel.U[0]+rotational_deltav;
 	sliding_velocity -= dot(sliding_velocity, interaction->nvec)*interaction->nvec;

	return sliding_velocity;
}

void Contact::incrementTangentialDisplacement(double dt, const struct PairVelocity &vel)
{
	/** Computes the tangential surface velocity difference between the two particles as
	delta_v = (surface_velo_p1 - surface_velo_p0).(Identity - nvec nvec),
	and increments the tangential spring stretch xi += delta_v*dt

	The surface velocity of p1 must be computed taking into account the Lees-Edwards periodic
	boundary conditions.

	 if p1 is upper, zshift = -1.
	 zshift = -1; //  p1 (z ~ lz), p0 (z ~ 0)

	******************************************************/
	disp_tan += getSlidingVelocity(vel)*dt;
}

vec3d Contact::getRollingVelocity(const struct PairVelocity &vel) const
{
	/**
	 Calculate rolling velocity
	 We follow Luding(2008).
	 The factor 2 is added.
	 cf. Book by Marshall and Li
	 equation 3.6.13 ??
	 */
	return 2*a_reduced*cross(vel.O[1]-vel.O[0], interaction->nvec);
}

void Contact::incrementRollingDisplacement(double dt, const struct PairVelocity &vel)
{
	disp_rolling += getRollingVelocity(vel)*dt; // always disp(t+1) = disp(t) + v*dt, no predictor-corrector headache :)
}

void Contact::incrementDisplacements(double dt, const struct PairVelocity &vel)
{
	/**
	 \brief Increment the tangential and rolling spring stretches from relative velocities, @b without checking the friction laws.

	 This should be called @b BEFORE updating the relative positions (ie normal and tangential vectors in the interaction).
	 This is because it needs the relative velocities at time t, which depend on a variable zshift at time t which deals with Lees-Edwards PBC.
	 This zshift is updated to their value at time t+1 whenever the relative positions are computed, so updating relative positions should be done after incrementing stretches.
	 */
	if (friction_model != FrictionModel::frictionless) {
		incrementTangentialDisplacement(dt, vel);
		if (rolling_friction()) {
			incrementRollingDisplacement(dt, vel);
		}
	}
}

void Contact::calcContactSpringForce()
{
	/**
		\brief Compute the contact forces and apply friction law, by rescaling forces and stretches if necessary.

	 f_something is the force acting on the particle 0.
	 disp_something is the relative displacement of the particle 1 from the particle 0.
	 Therefore, the sign of force is same as the one of the displacement.
	 */

	/* h < 0
	 * f_spring_normal_norm > 0 ..... repulsive force
	 * h > 0
	 * f_spring_normal_norm < 0 ..... attractive force
	 */
	f_spring_normal_norm = -kn_scaled*interaction->getReducedGap();
	f_spring_normal = -f_spring_normal_norm*interaction->nvec;
	if (friction_model != FrictionModel::frictionless) {
		disp_tan.vertical_projection(interaction->nvec);
		f_spring_tan = kt_scaled*disp_tan;
		if (rolling_friction()) {
			f_rolling = kr_scaled*disp_rolling;
		}
		(this->*frictionlaw)();
	}
	if (state == 1) {
		f_spring_total = f_spring_normal;
	} else {
		f_spring_total = f_spring_normal+f_spring_tan;
	}
}

vec3d Contact::getTotalForce(const PairVelocity &pvel) const
{
	/**
		\brief Compute the total contact forces (spring+dashpot).
		(This contains the dashpot, it is NOT only a static force.)
		@@@ Actually, the velocities are not correct after retrimming in the extensioanl flow simulation.
		@@@ It seems this affects only visualization data
		*/
	if (dashpot) {
		return f_spring_total + dashpot->getForceOnP0(pvel);
	} else {
		return f_spring_total;
	}
}

vec3d Contact::getSpringForce() const
{
	return f_spring_total;
}

void Contact::frictionlaw_standard()
{
	/**
	 \brief Friction law
	 */
	double supportable_tanforce = 0;
	double supportable_rollingforce = 0;
	double sq_f_tan;
	double sq_f_rolling;
	normal_load = f_spring_normal_norm;
	if (adhesion > 0) {
		normal_load += adhesion;
	}
	if (normal_load > 0) {
		if (state == 2) {
			// static friction in previous step
			supportable_tanforce = mu_static*normal_load;
		} else {
			// dynamic friction in previous step
			supportable_tanforce = mu_dynamic*normal_load;
		}
		if (rolling_friction()) {
			supportable_rollingforce = mu_rolling*normal_load;
		}
	}
	sq_f_tan = f_spring_tan.sq_norm();
	if (sq_f_tan < supportable_tanforce*supportable_tanforce) {
		state = 2; // static friction
	} else {
		state = 3; // dynamic friction
		supportable_tanforce = mu_dynamic*normal_load;
		// adjust the sliding spring for dynamic friction law
		setTangentialForceNorm(sqrt(sq_f_tan), supportable_tanforce);
	}
	if (rolling_friction()) {
		sq_f_rolling = f_rolling.sq_norm();
		if (sq_f_rolling > supportable_rollingforce*supportable_rollingforce) {
			setRollingForceNorm(sqrt(sq_f_rolling), supportable_rollingforce);
		}
	}
}

void Contact::frictionlaw_infinity()
{
	/**
	 \brief Friction law
	 */
	state = 2; // static friction
}

void Contact::setTangentialForceNorm(double current_force_norm,
									 double new_force_norm)
{
	disp_tan *= new_force_norm/current_force_norm;
	f_spring_tan = kt_scaled*disp_tan;
}

void Contact::setRollingForceNorm(double current_force_norm,
								  double new_force_norm)
{
	disp_rolling *= new_force_norm/current_force_norm;
	f_rolling = kr_scaled*disp_rolling;
}

void Contact::frictionlaw_criticalload()
{
	/* Since reduced_gap < 0, f_spring_normal_norm is always positive.
	 * f_spring_normal_norm = -kn_scaled*interaction->getReducedGap(); > 0
	 * F_normal = f_spring_normal_norm(positive) + dashpot_force_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	normal_load = f_spring_normal_norm-critical_load; // critical load model.
    if (normal_load < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_spring_tan.reset();
		if (rolling_friction()) {
			disp_rolling.reset();
			f_rolling.reset();
		}
	} else {
		double supportable_tanforce = mu_static*normal_load;
		double sq_f_tan = f_spring_tan.sq_norm();
		if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
			state = 3; // sliding
			setTangentialForceNorm(sqrt(sq_f_tan), supportable_tanforce);
		} else {
			state = 2; // static friction
		}
		if (rolling_friction()) {
			double supportable_rollingforce = mu_rolling*normal_load;
			double sq_f_rolling = f_rolling.sq_norm();
			if (sq_f_rolling > supportable_rollingforce*supportable_rollingforce) {
				setRollingForceNorm(sqrt(sq_f_rolling), supportable_rollingforce);
			}
		}
	}
}

void Contact::frictionlaw_criticalload_mu_inf()
{
	/* Since reduced_gap < 0, f_spring_normal_norm is always positive.
	 * f_spring_normal_norm = -kn_scaled*interaction->getReducedGap(); > 0
	 * F_normal = f_spring_normal_norm(positive) + dashpot_force_p0_normal
	 *
	 * supportable_tanforce = mu*(F_normal - critical_force)
	 *
	 */
	double supportable_tanforce = f_spring_normal_norm-critical_load; // critical load model.
	if (supportable_tanforce < 0) {
		state = 1; // frictionless contact
		disp_tan.reset();
		f_spring_tan.reset();
	} else {
		// Never rescaled.
		// [note]
		// The tangential spring constant may be rescaled to control maximum strain
		state = 2; // static friction
	}
}

void Contact::frictionlaw_ft_max()
{
	/**
	 \brief Friction law
	 */
	double sq_f_tan = f_spring_tan.sq_norm();
	if (sq_f_tan > ft_max*ft_max) {
		state = 3; // dynamic friction
		setTangentialForceNorm(sqrt(sq_f_tan), ft_max);
	} else {
		state = 2; // static friction
	}
}

void Contact::frictionlaw_coulomb_max()
{
	/**
	 \brief Friction law
	 *
	 * [NOTE]
	 * In the rate controlled simulation, the results look weird.
	 * There should be a bug.
	 *
	 * The calculated forces in calcContactInteraction() seem to be different
	 * between the rate-controlled and stress-controlled simulation.
	 * The current implementation is quite confusing, and should be fixed.
	 */
	double supportable_tanforce = 0;
	normal_load = f_spring_normal_norm;
	if (state == 2) {
		// static friction in previous step
		supportable_tanforce = mu_static*normal_load;
	} else if (state == 3) {
		// dynamic friction in previous step
		supportable_tanforce = mu_dynamic*normal_load;
	} else {
		std::cerr << "contact state" << state << std::endl;
		exit(1);
	}
	if (supportable_tanforce > ft_max) {
	 	supportable_tanforce = ft_max;
	}
	double sq_f_tan = f_spring_tan.sq_norm();
	if (sq_f_tan > supportable_tanforce*supportable_tanforce) {
		state = 3; // dynamic friction
		setTangentialForceNorm(sqrt(sq_f_tan), supportable_tanforce);
	} else {
		state = 2; // static friction
	}
}

void Contact::addUpForce(std::vector<vec3d> &force_per_particle) const
{
	/* Force
	 */
	force_per_particle[interaction->p0] += f_spring_total;
	force_per_particle[interaction->p1] -= f_spring_total;
}

void Contact::addUpForceTorque(std::vector<vec3d> &force_per_particle,
							   std::vector<vec3d> &torque_per_particle) const
{
	/* Force
	 */
	addUpForce(force_per_particle);
	/* Torque
	 */
	if (state >= 2) {
		vec3d t_ij = cross(interaction->nvec, f_spring_tan);
		torque_per_particle[interaction->p0] += interaction->a0*t_ij;
		torque_per_particle[interaction->p1] += interaction->a1*t_ij;
		if (rolling_friction()) {
			vec3d t_rolling = cross(interaction->nvec, f_rolling);
			torque_per_particle[interaction->p0] += interaction->a0*t_rolling;
			torque_per_particle[interaction->p1] -= interaction->a1*t_rolling;
		}
	}
}

void Contact::calcContactStress(const PairVelocity &pvel)
{
	/**
	 * The "xF" contact stress.
	 * The contact force F includes both the spring force and the dashpot force.
	 */
	contact_stresslet_XF = outer_sym(interaction->rvec, getTotalForce(pvel));
		// contact_stresslet_XF = outer_sym(interaction->rvec, f_spring_normal);
}

void Contact::addUpStress(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1,
						  const PairVelocity &pvel)
{
	calcContactStress(pvel);
	double r_ij = interaction->a0 + interaction->a1;
	stress_p0 += (interaction->a0/r_ij)*contact_stresslet_XF;
	stress_p1 += (interaction->a1/r_ij)*contact_stresslet_XF;
}

void Contact::addUpStressSpring(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1) const
{
	Sym2Tensor spring_stress;
	spring_stress = outer_sym(interaction->rvec, f_spring_total);
	double r_ij = interaction->a0 + interaction->a1;
	stress_p0 += (interaction->a0/r_ij)*spring_stress;
	stress_p1 += (interaction->a1/r_ij)*spring_stress;
}

void Contact::addUpStressDashpot(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1, 
								 const PairVelocity &pvel) const
{
	Sym2Tensor dashpot_stress;
	dashpot_stress = outer_sym(interaction->rvec, dashpot->getForceOnP0(pvel));
	double r_ij = interaction->a0 + interaction->a1;
	stress_p0 += (interaction->a0/r_ij)*dashpot_stress;
	stress_p1 += (interaction->a1/r_ij)*dashpot_stress;
}

double Contact::calcEnergy() const
{
	double overlap = -interaction->getReducedGap();
	double sq_tan_norm = disp_tan.sq_norm();
	double sq_disp_rolling = disp_rolling.sq_norm();
	/* normal */
	double energy = 0.5*kn_scaled*overlap*overlap;
	if (state >= 2) {
		/* sliding */
		energy += 0.5*kt_scaled*sq_tan_norm;
		if (rolling_friction()) {
			/* roling */
			energy += 0.5*kr_scaled*sq_disp_rolling;
		}
	}
	return energy;
}

double Contact::getNormalForceValue(const struct PairVelocity &vel) const
{
	return dot(interaction->nvec, getTotalForce(vel));
}

double Contact::getNormalSpringForce() const
{
	/* h < 0
	 * f_spring_normal_norm > 0 ..... repulsive force
	 * h > 0
	 * f_spring_normal_norm < 0 ..... attractive force
	 */
	return f_spring_normal_norm;
}

double Contact::get_normal_load() const
{
	return normal_load;
}

vec3d Contact::getTangentialForce(const struct PairVelocity &vel) const
{
	vec3d total_force = getTotalForce(vel);
	return total_force-dot(interaction->nvec, total_force)*interaction->nvec;
}


struct contact_state Contact::getState() const
{
	struct contact_state cs;
	cs.p0 = interaction->p0;
	cs.p1 = interaction->p1;
	cs.disp_tan = disp_tan;
	cs.disp_rolling = disp_rolling;
	return cs;
}

void Contact::setState(const struct contact_state& cs)
{
	disp_tan = cs.disp_tan;
	disp_rolling = cs.disp_rolling;
}

void Contact::saveState()
{
	prev_disp_tan = disp_tan;
	prev_disp_rolling = disp_rolling;
}

void Contact::restoreState()
{
	disp_tan = prev_disp_tan;
	disp_rolling = prev_disp_rolling;
}

} // namespace Interactions