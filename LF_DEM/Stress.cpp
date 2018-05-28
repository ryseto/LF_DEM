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
#include "global.h"
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
	if (lubrication && na_velo_components.empty()) {
		throw runtime_error(" System::declareStressComponents: No velocity components declared, you probably forgot it.");
	}
	if (lubrication) {
		for (const auto &vc: na_velo_components) {
			stress_components[vc.first] = StressComponent(VELOCITY_STRESS,
														  vc.second.vel.size(),
														  vc.second.rate_dependence,
														  vc.first);
		}
	}

	// Brownian
	if (brownian) { // Brownian is different than other GU, needs predictor data too
		stress_components["brownian_predictor"] = StressComponent(BROWNIAN_STRESS,
		                                                          np, RATE_INDEPENDENT, "brownian");// rate dependent for now --> @ WORKING NOW @
	}

	/****************  ME stress ********************/
	if (lubrication) {
		stress_components["M_E_hydro"] = StressComponent(STRAIN_STRESS, np, RATE_PROPORTIONAL, "hydro");
	}

	/****************  xF stresses *****************/
	if (control == Parameters::ControlVariable::rate) {
		stress_components["xF_contact"] = StressComponent(XF_STRESS, np, RATE_DEPENDENT, "contact"); // rate dependent through xFdashpot
	} else { // stress controlled
		stress_components["xF_contact_rateprop"] = StressComponent(XF_STRESS, np, RATE_PROPORTIONAL, "contact");
		stress_components["xF_contact_rateindep"] = StressComponent(XF_STRESS, np, RATE_INDEPENDENT, "contact");
	}
	if (repulsiveforce) {
		stress_components["xF_repulsion"] = StressComponent(XF_STRESS, np, RATE_INDEPENDENT, "repulsion");
	}

	if (delayed_adhesion) {
		stress_components["xF_delayed_adhesion"] = StressComponent(XF_STRESS, np, RATE_INDEPENDENT, "delayed_adhesion");
	}

	if (control == Parameters::ControlVariable::stress) {
		for (const auto &sc: stress_components) {
			if (sc.second.rate_dependence == RATE_DEPENDENT) {
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

void System::addUpInteractionStressGU(std::vector<Sym2Tensor> &stress_comp,
                                      const std::vector<vec3d> &non_affine_vel,
                                      const std::vector<vec3d> &non_affine_ang_vel)
{
	if (!lubrication) {
		return;
	}

	for (const auto &inter: interaction) {
		if (inter.lubrication.is_active()) {
			unsigned int i, j;
			std::tie(i, j) = inter.get_par_num();
			inter.lubrication.addGUStresslet(non_affine_vel[i], non_affine_vel[j],
											 non_affine_ang_vel[i], non_affine_ang_vel[j],
											 stress_comp[i], stress_comp[j]);
		}
	}
}

void System::addUpInteractionStressME(std::vector<Sym2Tensor> &stress_comp)
{
	if (!lubrication) {
		return;
	}
	for (const auto &inter: interaction) {
		if (inter.lubrication.is_active()) {
			unsigned int i, j;
			std::tie(i, j) = inter.get_par_num();
			inter.lubrication.addMEStresslet(E_infinity,
											 stress_comp[i],
											 stress_comp[j]); // R_SE:Einf-R_SU*v
		}
	}
}

void System::gatherVelocitiesByRateDependencies(vector<vec3d> &rateprop_vel,
                                                vector<vec3d> &rateprop_ang_vel,
                                                vector<vec3d> &rateindep_vel,
                                                vector<vec3d> &rateindep_ang_vel) const
{
	/** Gather velocity components in rate proportional and rate independent parts.
			If there is a rate dependent (but not proportional), it is left out.

			[Note] The rate proportional part is a total velocity, and as such must be dealt with proper Lees-Edwards BC.
						 The rate independent part is always a non-affine velocity.
	*/
	assert(control == Parameters::ControlVariable::stress);// || control == viscnb);
	for (unsigned int i=0; i<rateprop_vel.size(); i++) {
		rateprop_vel[i].reset();
		rateprop_ang_vel[i].reset();
		rateindep_vel[i].reset();
		rateindep_ang_vel[i].reset();
	}
	for (const auto &vc: na_velo_components) {
		if (vc.second.rate_dependence == RATE_PROPORTIONAL) {
			for (unsigned int i=0; i<rateprop_vel.size(); i++) {
				rateprop_vel[i] += vc.second.vel[i];
				rateprop_ang_vel[i] += vc.second.ang_vel[i];
			}
		}
		if (vc.second.rate_dependence == RATE_INDEPENDENT) {
			for (unsigned int i=0; i<rateprop_vel.size(); i++) {
				rateindep_vel[i] += vc.second.vel[i];
				rateindep_ang_vel[i] += vc.second.ang_vel[i];
			}
		}
	}
	for (unsigned int i=0; i<rateprop_vel.size(); i++) {
		rateprop_vel[i] += u_inf[i];
		rateprop_ang_vel[i] += omega_inf;
	}
}

void System::calcContactXFPerParticleStressControlled()
{
	// spring part: easy
	// dashpot part: we have to split between rate proportional and rate independent parts.
	// a bit annoying :)

	vector<vec3d> rateprop_vel (np);
	vector<vec3d> rateprop_ang_vel (np);
	vector<vec3d> rateindep_vel (np);
	vector<vec3d> rateindep_ang_vel (np);
	gatherVelocitiesByRateDependencies(rateprop_vel, rateprop_ang_vel,
	                                   rateindep_vel, rateindep_ang_vel);

	auto &rateprop_XF = stress_components.at("xF_contact_rateprop").particle_stress;
	auto &rateindep_XF = stress_components.at("xF_contact_rateindep").particle_stress;

	for (const auto &inter: interaction) {
		unsigned int i, j;
		std::tie(i, j) = inter.get_par_num();
		if (inter.contact.is_active()) {
			inter.contact.addUpStressSpring(rateindep_XF[i], rateindep_XF[j]); // - rF_cont
		}
		if (inter.contact.dashpot.is_active()) {
			// rate_prop_vel is a full velocity (not non affine)
			Sym2Tensor rateprop_stress = outer_sym(inter.rvec,
			                                       inter.contact.dashpot.getForceOnP0(rateprop_vel[i],
			                                                                          rateprop_vel[j],
			                                                                          rateprop_ang_vel[i],
			                                                                          rateprop_ang_vel[j]));
			double r_ij = radius[i]+radius[j];
			rateprop_XF[i] += (radius[i]/r_ij)*rateprop_stress;
			rateprop_XF[j] += (radius[j]/r_ij)*rateprop_stress;

			Sym2Tensor rateindep_stress = outer_sym(inter.rvec,
													inter.contact.dashpot.getForceOnP0_nonaffine(rateindep_vel[i],
																								 rateindep_vel[j],
			                                                                                     rateindep_ang_vel[i],
			                                                                                     rateindep_ang_vel[j]));
			rateindep_XF[i] += (radius[i]/r_ij)*rateindep_stress;
			rateindep_XF[j] += (radius[j]/r_ij)*rateindep_stress;
		}
	}
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
    if (lubrication) {
        for (auto &inter: interaction) {
            if (!inter.lubrication.tangential) {
                inter.lubrication.calcXFunctionsStress();
            } else {
                inter.lubrication.calcXYFunctionsStress();
            }
        }
    }

	for (auto &sc: stress_components) {
		sc.second.reset();
	}

    for (auto &sc: stress_components) {
        auto type = sc.second.type;
        const auto &component_name = sc.first;
        if (type == VELOCITY_STRESS) {
            addUpInteractionStressGU(sc.second.particle_stress,
                                     na_velo_components[component_name].vel,
                                     na_velo_components[component_name].ang_vel);
        }
        if (type == STRAIN_STRESS) {
            addUpInteractionStressME(sc.second.particle_stress);
        }
    }
    if (control == Parameters::ControlVariable::rate) {
        auto &cstress_XF = stress_components.at("xF_contact").particle_stress;
		for (auto &inter: interaction) {
			if (inter.contact.is_active()) {
				unsigned int i, j;
				std::tie(i, j) = inter.get_par_num();
				inter.contact.addUpStress(cstress_XF[i], cstress_XF[j]); // - rF_cont
			}
		}
	} else {
		calcContactXFPerParticleStressControlled();
	}

	if (repulsiveforce) {
		auto &rstress_XF = stress_components.at("xF_repulsion").particle_stress;
		for (auto &inter: interaction) {
			unsigned int i, j;
			std::tie(i, j) = inter.get_par_num();
			inter.repulsion.addUpStressXF(rstress_XF[i], rstress_XF[j]); // - rF_rep
		}
	}

	if (delayed_adhesion) {
		auto &rstress_XF = stress_components.at("xF_delayed_adhesion").particle_stress;
		for (auto &inter: interaction) {
			unsigned int i, j;
			std::tie(i, j) = inter.get_par_num();
			inter.delayed_adhesion->addUpStressXF(rstress_XF[i], rstress_XF[j], inter.rvec); // - rF_rep
		}
	}

	if (brownian) {
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
	vec3d pos_normal = position[i]-origin_of_rotation;
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
		if (sc.second.rate_dependence == RATE_INDEPENDENT) {
			rate_indep_stress += sc.second.getTotalStress();
		}
		if (sc.second.rate_dependence == RATE_PROPORTIONAL) {
			rate_prop_stress += sc.second.getTotalStress();
		}
	}
	rate_prop_stress /= system_volume;
	rate_indep_stress /= system_volume;

	if (!zero_shear) {
		// suspending fluid viscosity
		rate_prop_stress += 2*E_infinity/(6*M_PI);
	}
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
	for (auto &sc: total_stress_groups) {
		sc.second /= system_volume;
	}

	if (!zero_shear) {
		total_stress_groups["hydro"] += 2*E_infinity/(6*M_PI);
	}

	total_stress.reset();
	for (const auto &sc: total_stress_groups) {
		total_stress += sc.second;
	}
	if (brownian && brownian_dominated) {
		// take an averaged stress instead of instantaneous
		stress_avg.update(total_stress, get_time());
		total_stress = stress_avg.get();
	}

    if (wall_rheology) {
        if (z_top != -1) {
            shearstress_wall1 = force_tang_wall1/lx;
            shearstress_wall2 = force_tang_wall2/lx;
            normalstress_wall1 = force_normal_wall1/lx;
            normalstress_wall2 = force_normal_wall2/lx;
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
