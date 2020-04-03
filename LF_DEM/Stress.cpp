//
//  Stress.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <sstream>
#include <stdexcept>
#include <assert.h>
#include "System.h"
#include "StressComponent.h"
using namespace std;

void System::declareStressComponents()
{
	// It is essential that the stresses in stress_components sum up to the total stress.
	// exemple: if the total stress is S_TOTAL = S_A + S_B,
	// you should declare A and B, but not TOTAL.
	// If you need S_TOTAL for some other need, please use an other structure to store it.
	// The rate dependence is used only for the stress control algorithm.

	/*****************  GU stresses *********************/
	// From the velocity components
	if (Interactions::hasPairwiseResistanceStdInteraction(p)) {
		if (na_velo_components.empty()) {
			throw runtime_error(" System::declareStressComponents: No velocity components declared, you probably forgot it.");
		}
		for (const auto &vc: na_velo_components) {
			stress_components[vc.first] = StressComponent(StressType::velocity,
														  vc.second.vel.size(),
														  vc.second.rate_dependence,
														  vc.first);
		}
	}

	// Brownian
	if (p.brownian > 0) { // Brownian is different than other GU, needs predictor data too
		stress_components["brownian_predictor"] = StressComponent(StressType::brownian,
																  np, RateDependence::independent, "brownian");// rate dependent for now --> @ WORKING NOW @
	}

	/****************  ME stress ********************/
	if (Interactions::has_lubrication(p.lub)) {
		stress_components["M_E_hydro"] = StressComponent(StressType::velocitygrad, np, RateDependence::proportional, "hydro");
	}

	/****************  xF stresses *****************/
	if (control == Parameters::ControlVariable::rate) {
		stress_components["xF_contact"] = StressComponent(StressType::xf, np, RateDependence::dependent, "contact"); // rate dependent through xFdashpot
	} else { // stress controlled
		stress_components["xF_contact_rateprop"] = StressComponent(StressType::xf, np, RateDependence::proportional, "contact");
		stress_components["xF_contact_rateindep"] = StressComponent(StressType::xf, np, RateDependence::independent, "contact");
	}
	if (Interactions::has_repulsion(p.repulsion)) {
		stress_components["xF_repulsion"] = StressComponent(StressType::xf, np, RateDependence::independent, "repulsion");
	}

	// if (delayed_adhesion) {
	// 	stress_components["xF_delayed_adhesion"] = StressComponent(StressType::xf, np, RateDependence::independent, "delayed_adhesion");
	// }

	if (p.confinement.on) {
		stress_components["xF_confinement"] = StressComponent(StressType::xf, np, RateDependence::independent, "confinement");
	}

	if (dimer_manager) {
		if (control == Parameters::ControlVariable::rate) {
			stress_components["xF_dimer"] = StressComponent(StressType::xf, np, RateDependence::dependent, "dimer");
		} else if (control == Parameters::ControlVariable::stress) {
			stress_components["xF_dimer_rateprop"] = StressComponent(StressType::xf, np, RateDependence::proportional, "dimer");
			stress_components["xF_dimer_rateindep"] = StressComponent(StressType::xf, np, RateDependence::independent, "dimer");
		} else {
			throw runtime_error("Dimer stress component when control is not rate not stress?");
		}
	}

	if (control == Parameters::ControlVariable::stress) {
		for (const auto &sc: stress_components) {
			if (sc.second.rate_dependence == RateDependence::dependent) {
				ostringstream error_msg;
				error_msg << "Cannot run stress controlled simulation with stress component " << sc.first;
				error_msg << " as it is neither independent nor proportional to the shear rate.";
				throw runtime_error(error_msg.str());
			}
		}
	}

	for (const auto &sc: stress_components) {
		auto &group = sc.second.group;
		total_stress_groups[group] = Sym2Tensor();
	}
}

void System::gatherVelocitiesByRateDependencies(ParticleVelocity &rateprop_vel,
												ParticleVelocity &rateindep_vel) const
{
	/** Gather velocity components in rate proportional and rate independent parts.
			If there is a rate dependent (but not proportional), it is left out.

			[Note] The rate proportional part is a total velocity, and as such must be dealt with proper Lees-Edwards BC.
						 The rate independent part is always a non-affine velocity.
	*/
	assert(control == Parameters::ControlVariable::stress);// || control == viscnb);
	rateprop_vel.reset();
	rateindep_vel.reset();
	for (const auto &vc: na_velo_components) {
		if (vc.second.rate_dependence == RateDependence::proportional) {
			for (unsigned int i=0; i<rateprop_vel.vel.size(); i++) {
				rateprop_vel.vel[i] += vc.second.vel[i];
				rateprop_vel.ang_vel[i] += vc.second.ang_vel[i];
			}
		}
		if (vc.second.rate_dependence == RateDependence::independent) {
			for (unsigned int i=0; i<rateprop_vel.vel.size(); i++) {
				rateindep_vel.vel[i] += vc.second.vel[i];
				rateindep_vel.ang_vel[i] += vc.second.ang_vel[i];
			}
		}
	}
	if (vel_bg.rate_dependence == RateDependence::proportional) {
		for (unsigned int i=0; i<rateprop_vel.vel.size(); i++) {
			rateprop_vel.vel[i] += vel_bg.vel[i];
			rateprop_vel.ang_vel[i] += vel_bg.ang_vel[i];
		}
	}
	if (vel_bg.rate_dependence == RateDependence::independent) {
		for (unsigned int i=0; i<rateprop_vel.vel.size(); i++) {
			rateindep_vel.vel[i] += vel_bg.vel[i];
			rateindep_vel.ang_vel[i] += vel_bg.ang_vel[i];
		}
	}
}

void System::calcDimerXFPerParticleStressControlled()
{
	// spring part: easy
	// dashpot part: we have to split between rate proportional and rate independent parts.
	// a bit annoying :)

	ParticleVelocity rateprop_vel (np, VelocityType::total, RateDependence::proportional);
	ParticleVelocity rateindep_vel (np, VelocityType::nonaffine, RateDependence::independent);
	gatherVelocitiesByRateDependencies(rateprop_vel, rateindep_vel);

	auto &rateprop_XF = stress_components.at("xF_dimer_rateprop").particle_stress;
	auto &rateindep_XF = stress_components.at("xF_dimer_rateindep").particle_stress;

	dimer_manager->addUpStressSpring(rateindep_XF);
	dimer_manager->addUpStressDashpot(rateindep_XF, &rateindep_vel);
	dimer_manager->addUpStressDashpot(rateprop_XF, &rateprop_vel);
}

void System::calcContactXFPerParticleStressControlled()
{
	// spring part: easy
	// dashpot part: we have to split between rate proportional and rate independent parts.
	// a bit annoying :)

	ParticleVelocity rateprop_vel (np, VelocityType::total, RateDependence::proportional);
	ParticleVelocity rateindep_vel (np, VelocityType::nonaffine, RateDependence::independent);
	gatherVelocitiesByRateDependencies(rateprop_vel, rateindep_vel);

	auto &rateprop_XF = stress_components.at("xF_contact_rateprop").particle_stress;
	auto &rateindep_XF = stress_components.at("xF_contact_rateindep").particle_stress;

	interaction->addUpContactSpringStressXF(rateindep_XF);
	interaction->addUpContactDashpotStressXF(rateindep_XF, &rateindep_vel);
	interaction->addUpContactDashpotStressXF(rateprop_XF, &rateprop_vel);
}


void System::calcStressPerParticle()
{
	/**
	   \brief This method computes the stresses per particle, split by components (hydro, contact, ...).

	   From velocities \f$ V_{\mathrm{I}}\f$ associated with
	   interaction \f$\mathrm{I}\f$, this method gets the stresses \f$ - GV_{\mathrm{I}} \f$. (This corresponds to
	   \f$- GU_{\mathrm{I}} - H\Omega_{\mathrm{I}} \f$ in Jeffrey
	   notations \cite jeffrey_calculation_1992, and
	   \f$- R_{\mathrm{SU}} U_{\mathrm{I}} \f$ in Bossis and Brady
	   \cite brady_stokesian_1988 notations.)

	   For the hydrodynamic component, it also gets the \f$ M
	   E_{\infty}\f$ term (\f$R_{\mathrm{SE}} E_{\infty}\f$ is B&B
	   notations), so that all in all \f$ S_{\mathrm{H}} =
	   - GV_{\mathrm{H}} + M E_{\infty}\f$.

	   For the point forces, it also gets the \f$ -xF_{\mathrm{I}} \f$ term, so that \f$ S_{\mathrm{I}} =
	   - GV_{\mathrm{I}} - xF_{\mathrm{I}} \f$.

	   For the Brownian forces, it computes (in B&B notations) \f$ S_{\mathrm{B}} =
	   - kT \nabla\dot (R_{\mathrm{SU}}.R_{\mathrm{FU}}^{-1}) \f$ with the mid-step algorithm of Banchio and Brady
	   \cite banchio_accelerated_2003 with \f$ n=1 \f$.

	   In the Brownian mode, because of the mid-point scheme for the Brownian stress, you
	   should be careful when calling this method from outside of the
	   System::timeEvolutionPredictorCorrectorMethod method.
	*/
	for (auto &sc: stress_components) {
		sc.second.reset();
	}

	for (auto &sc: stress_components) {
		auto type = sc.second.type;
		const auto &component_name = sc.first;
		if (type == StressType::velocity) {
			interaction->addUpInteractionStressGU(sc.second.particle_stress,
									 			  &na_velo_components[component_name]);
		}
		if (type == StressType::velocitygrad) {  // there is only one
			interaction->addUpInteractionStressME(sc.second.particle_stress);
		}
	}
	if (control == Parameters::ControlVariable::rate) {
		auto &cstress_XF = stress_components.at("xF_contact").particle_stress;
		interaction->addUpContactStressXF(cstress_XF, &velocity);
	} else {
		calcContactXFPerParticleStressControlled();
	}

	if (Interactions::has_repulsion(p.repulsion)) {
		auto &rstress_XF = stress_components.at("xF_repulsion").particle_stress;
		interaction->addUpRepulsiveStressXF(rstress_XF);
	}

	// if (delayed_adhesion) {
	// 	auto &rstress_XF = stress_components.at("xF_delayed_adhesion").particle_stress;
	// 	addUpDelayedAdhesionStressXF(rstress_XF);
	// }

	if (dimer_manager) {
		if (control == Parameters::ControlVariable::rate) {
			auto &dstress_XF = stress_components.at("xF_dimer").particle_stress;
			dimer_manager->addUpStress(dstress_XF, &velocity);
		} else {
			calcDimerXFPerParticleStressControlled();
		}
	}
	if (p.confinement.on) {
		vec3d yvec = {0, 1, 0};
		auto &stress_XF = stress_components.at("xF_confinement").particle_stress;
		auto &force = force_components.at("confinement").force;

		for (unsigned i=0; i<stress_XF.size(); i++) {
			if (force[i].y > 0) { // boundary at y_min
				stress_XF[i] += outer_sym(-conf->radius[i]*yvec, force[i]);
			} else {                                       // boundary at y_max
				stress_XF[i] += outer_sym(conf->radius[i]*yvec, force[i]);
			}
		}
	}

	if (p.brownian > 0) {
		auto &bstress_predictor = stress_components.at("brownian_predictor").particle_stress;
		auto &bstress = stress_components.at("brownian").particle_stress;

		if (in_predictor) {
			for (unsigned int i=0; i<bstress.size(); i++) {
				bstress_predictor[i] = bstress[i];
			}
		} else {
			for (unsigned int i=0; i<bstress.size(); i++) {
				/*
				 * [ Banchio & Brady 2003 ] [ Ball & Melrose 1997 ]
				 */
				bstress[i] = 0.5*(bstress[i]-bstress_predictor[i]);
			}
		}
	}
}

void System::calcTotalStressPerParticle()
{
	for (unsigned int i=0; i<total_stress_pp.size(); i++) {
		total_stress_pp[i].reset();
	}
	for (const auto &sc: stress_components) {
		const auto &particle_stress = sc.second.particle_stress;
		for (unsigned int i=0; i<total_stress_pp.size(); i++) {
			total_stress_pp[i] += particle_stress[i];
		}
	}
}

void System::getStressCouette(int i,
							  double &stress_rr,
							  double &stress_tt,
							  double &stress_rt)
{
	// (0 xx, 1 xy, 2 xz, 3 yz, 4 yy, 5 zz)
	vec3d pos_normal = conf->position[i]-origin_of_rotation;
	double r_dist = pos_normal.norm();
	double ctheta = pos_normal.x/r_dist;
	double stheta = pos_normal.z/r_dist;
	double cc = ctheta*ctheta;
	double ss = stheta*stheta;
	double cs = ctheta*stheta;
	stress_rr =  cc*total_stress_pp[i].elm[0]+   2*cs*total_stress_pp[i].elm[2]+ss*total_stress_pp[i].elm[5];
	stress_tt =  ss*total_stress_pp[i].elm[0]-   2*cs*total_stress_pp[i].elm[2]+cc*total_stress_pp[i].elm[5];
	stress_rt = -cs*total_stress_pp[i].elm[0]+(cc-ss)*total_stress_pp[i].elm[2]+cs*total_stress_pp[i].elm[5];
}

void System::gatherStressesByRateDependencies(Sym2Tensor &rate_prop_stress,
											  Sym2Tensor &rate_indep_stress)
{
	rate_prop_stress.reset();
	rate_indep_stress.reset();
	for (const auto &sc: stress_components) {
		if (sc.second.rate_dependence == RateDependence::independent) {
			rate_indep_stress += sc.second.getTotalStress();
		}
		if (sc.second.rate_dependence == RateDependence::proportional) {
			rate_prop_stress += sc.second.getTotalStress();
		}
	}
	double system_volume = getSystemVolume();
	rate_prop_stress /= system_volume;
	rate_indep_stress /= system_volume;

	// suspending fluid viscosity
	rate_prop_stress += 2*imposed_flow->getSymGradU()/(6*M_PI);
}

void System::calcStress()
{
	for (auto &sc: total_stress_groups) {
		sc.second.reset();
	}
	for (const auto &sc: stress_components) {
		auto &group = sc.second.group;
		total_stress_groups[group] += sc.second.getTotalStress();
	}
	double system_volume = getSystemVolume();
	for (auto &sc: total_stress_groups) {
		sc.second /= system_volume;
	}

	total_stress_groups["hydro"] += 2*imposed_flow->getSymGradU()/(6*M_PI);

	total_stress.reset();
	for (const auto &sc: total_stress_groups) {
		total_stress += sc.second;
	}
	if (p.brownian > 0) {
		// take an averaged stress instead of instantaneous
		stress_avg.update(total_stress, get_time());
		total_stress = stress_avg.get();
	}

	if (wall_rheology) {
		if (z_top != -1) {
			shearstress_wall1 = force_tang_wall1/container.lx;
			shearstress_wall2 = force_tang_wall2/container.lx;
			normalstress_wall1 = force_normal_wall1/container.lx;
			normalstress_wall2 = force_normal_wall2/container.lx;
		} else {
			double wall_area_in = M_PI*2*(radius_in-radius_wall_particle);
			double wall_area_out = M_PI*2*(radius_out+radius_wall_particle);
			shearstress_wall1 = force_tang_wall1/wall_area_in;
			shearstress_wall2 = force_tang_wall2/wall_area_out;
			normalstress_wall1 = force_normal_wall1/wall_area_in;
			normalstress_wall2 = force_normal_wall2/wall_area_out;
		}
	}
}
