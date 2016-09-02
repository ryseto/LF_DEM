//
//  Stress.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <sstream>
#include <stdexcept>
#include "System.h"
#include "global.h"
#include "StressComponent.h"


using namespace std;


void System::declareStressComponents() {
	// It is essential that the stresses in stress_components sum up to the total stress.
	// exemple: if the total stress is S_TOTAL = S_A + S_B,
	// you should declare A and B, but not TOTAL.
	// If you need S_TOTAL for some other need, please use an other structure to store it.
	// The rate dependence is used only for the stress control algorithm.

	/*****************  GU stresses *********************/
	// From the velocity components
	if (lubrication && velocity_components.empty()) {
		throw runtime_error(" System::declareStressComponents: No velocity components declared, you probably forgot it.");
	}
	if (lubrication) {
		for (const auto &vc: velocity_components) {
			if (vc.first != "brownian") {
				stress_components[vc.first] = StressComponent(VELOCITY_STRESS,
				                                                vc.second.vel.size(),
			                                                  vc.second.rate_dependence,
			                                                  vc.first);
			}
		}
	}

	// Brownian
	if (brownian) { // Brownian is different than other GU, needs predictor data too
		stress_components["brownian"] = StressComponent(BROWNIAN_STRESS,
		                                                   np, RATE_DEPENDENT, "brownian");// rate dependent for now
		stress_components["brownian_predictor"] = StressComponent(BROWNIAN_STRESS,
		                                                             np, RATE_DEPENDENT, "brownian");// rate dependent for now
	}

	/****************  ME stress ********************/
	if (lubrication) {
		stress_components["M_E_hydro"] = StressComponent(STRAIN_STRESS, np, RATE_PROPORTIONAL, "hydro");
	}

	/****************  xF stresses *****************/
	if (rate_controlled) {
		stress_components["xF_contact"] = StressComponent(XF_STRESS, np, RATE_DEPENDENT, "contact"); // rate dependent through xFdashpot
	} else { // stress controlled
		stress_components["xF_contact_rateprop"] = StressComponent(XF_STRESS, np, RATE_PROPORTIONAL, "contact");
		stress_components["xF_contact_rateindep"] = StressComponent(XF_STRESS, np, RATE_INDEPENDENT, "contact");
	}
	if (repulsiveforce) {
		stress_components["xF_repulsion"] = StressComponent(XF_STRESS, np, RATE_INDEPENDENT, "repulsion");
	}

	if (stress_controlled) {
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
		total_stress_groups[group] = StressTensor();
	}
}


void System::addUpInteractionStressGU(std::vector<StressTensor> &stress_comp,
                                      const std::vector<vec3d> &non_affine_vel,
                                      const std::vector<vec3d> &non_affine_ang_vel)
{
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			unsigned int i, j;
			std::tie(i, j) = interaction[k].get_par_num();
			if (lubrication) {
				if (interaction[k].lubrication.is_active()) {
					interaction[k].lubrication.addGUStresslet(non_affine_vel[i], non_affine_vel[j],
				                                            non_affine_ang_vel[i], non_affine_ang_vel[j],
				                                            stress_comp[i], stress_comp[j]);
				}
			}
		}
	}
}

void System::addUpInteractionStressME(std::vector<StressTensor> &stress_comp)
{
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			unsigned int i, j;
			std::tie(i, j) = interaction[k].get_par_num();
			if (lubrication) {
				if (interaction[k].lubrication.is_active()) {
					interaction[k].lubrication.addMEStresslet(costheta_shear,
					                                          sintheta_shear,
					                                          shear_rate,
					                                          stress_comp[i], stress_comp[j]); // R_SE:Einf-R_SU*v
				}
			}
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
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				if (!interaction[k].lubrication.tangential) {
					interaction[k].lubrication.calcXFunctionsStress();
				} else {
					interaction[k].lubrication.calcXYFunctionsStress();
				}
			}
		}
	}

	for(auto &sc: stress_components) {
		sc.second.reset();
	}

	for(auto &sc: stress_components) {
		auto type = sc.second.type;
		const auto &component_name = sc.first;
		if (type == VELOCITY_STRESS) {
			addUpInteractionStressGU(sc.second.particle_stress,
			                         velocity_components[component_name].vel,
															 velocity_components[component_name].ang_vel);
		}
		if (type == STRAIN_STRESS) {
			addUpInteractionStressME(sc.second.particle_stress);
		}
	}

	auto &cstress_XF = stress_components["xF_contact"].particle_stress;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			unsigned int i, j;
			std::tie(i, j) = interaction[k].get_par_num();
			if (interaction[k].contact.is_active()) {
				interaction[k].contact.addUpStress(cstress_XF[i], cstress_XF[j]); // - rF_cont
			}
		}
	}

	if (repulsiveforce) {
		auto &rstress_XF = stress_components["xF_repulsion"].particle_stress;
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				unsigned int i, j;
				std::tie(i, j) = interaction[k].get_par_num();
				// if (interaction[k].repulsion.is_active()) {
				interaction[k].repulsion.addUpStressXF(rstress_XF[i], rstress_XF[j]); // - rF_rep
				// }
			}
		}
	}

	if (brownian) {
		auto &bstress_predictor = stress_components["GU_brownian_predictor"].particle_stress;
		auto &bstress = stress_components["GU_brownian"].particle_stress;

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
		for (int i=0; i<np_mobile; i++) { // @@@@ np_mobile or total_stress_pp.size()??
			total_stress_pp[i] + particle_stress[i];
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

void System::gatherRateDependences()
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

	// suspending fluid viscosity
	if (p.cross_shear) {
		total_stress_groups["hydro"].elm[2] += costheta_shear*shear_rate/6./M_PI;
		total_stress_groups["hydro"].elm[3] += sintheta_shear*shear_rate/6./M_PI;
	}	else {
		total_stress_groups["hydro"].elm[2] += shear_rate/6./M_PI;
	}

	total_stress.reset();
	for (const auto &sc: total_stress_groups) {
		total_stress += sc.second;
	}
	if (brownian && lowPeclet) {
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
