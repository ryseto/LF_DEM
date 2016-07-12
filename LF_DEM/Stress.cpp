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

using namespace std;

void System::stressReset()
{
	/**
	   \brief Sets stresses arrays to zero.

	   To be called by System::calcStressPerParticle()
	*/
	for (int i=0; i<np; i++) {
		lubstress[i].reset();
		contactstressGU[i].reset();
	}
	if (mobile_fixed) {
		for (int i=0; i<np; i++) {
			hydrofromfixedstressGU[i].reset();
		}
	}
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			repulsivestressGU[i].reset();
		}
	}
	if (brownian) {
		for (int i=0; i<np; i++) {
			brownianstressGU[i].reset();
		}
	}
}

void System::calcStressPerParticle()
{
	/**
	   This method computes the stresses per particle, split by components (hydro, contact, ...).

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
	stressReset();
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (p.lubrication_model > 0) {
				if (p.lubrication_model == 1 || p.lubrication_model == 3) {
					interaction[k].lubrication.calcXFunctionsStress();
				} else if (p.lubrication_model == 2) {
					interaction[k].lubrication.calcXYFunctionsStress();
				} else {
					ostringstream error_str;
					error_str << "lubrication_model = " << p.lubrication_model << endl << "lubrication_model > 3 is not implemented" << endl;
					throw runtime_error(error_str.str());
				}
				interaction[k].lubrication.addStressME(); // R_SE:Einf-R_SU*v
				if (interaction[k].lubrication.is_active()) {
					interaction[k].lubrication.addStressesGU(); // R_SE:Einf-R_SU*v
				}
			}
			interaction[k].contact.calcContactStress(); // - rF_cont
			if (repulsiveforce) {
				interaction[k].repulsion.calcStressXF(); // - rF_rep
			}
		}
	}
	if (brownian) {
		if (in_predictor) {
			for (int i=0; i<np; i++) {
				brownianstressGU_predictor[i] = brownianstressGU[i];
			}
		} else {
			for (int i=0; i<np; i++) {
				/*
				 * [ Banchio & Brady 2003 ] [ Ball & Melrose 1997 ]
				 */
				brownianstressGU[i] = 0.5*(brownianstressGU[i]-brownianstressGU_predictor[i]);
			}
		}
	}
}

void System::calcTotalStressPerParticle()
{
	for (int i=0; i<np_mobile; i++) {
		total_stress_pp[i].reset();
		total_stress_pp[i] += lubstress[i];
		total_stress_pp[i] += (contactstressXF[i]+contactstressGU[i]);
		if (mobile_fixed) {
			total_stress_pp[i] += hydrofromfixedstressGU[i];
		}
		if (repulsiveforce) {
			total_stress_pp[i] += (repulsivestressXF[i]+repulsivestressGU[i]);
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

void System::calcStress()
{
	// Lubrication stress
	total_hydro_stress.reset();
	for (int i=0; i<np; i++) {
		total_hydro_stress += lubstress[i];
	}
	total_hydro_stress /= system_volume;
	total_hydro_stress += shear_rate/6./M_PI; // suspending fluid viscosity
	// Stress from contact force
	// GU contribution
	total_contact_stressGU.reset();
	for (int i=0; i<np; i++) {
		total_contact_stressGU += contactstressGU[i];
	}
	total_contact_stressGU /= system_volume;
	// XF contribution
	total_contact_stressXF.reset();
	bool cstress_per_particle = false;
	if (p.out_particle_stress.find('c') != string::npos || couette_stress) {
		cstress_per_particle = true;
		for (int i=0; i<np; i++) {
			contactstressXF[i].reset();
		}
	}
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].contact.is_active()) {
			total_contact_stressXF += interaction[k].contact.getContactStressXF();
			if (cstress_per_particle) {
				StressTensor sc = interaction[k].contact.getContactStressXF();
				unsigned int i, j;
				std::tie(i, j) = interaction[k].get_par_num();
				double r_ij = interaction[k].ro;
				contactstressXF[i] += (radius[i]/r_ij)*sc;
				contactstressXF[j] += (radius[j]/r_ij)*sc;
			}
		}
	}
	total_contact_stressXF /= system_volume;
	// Stress from repulsive force
	if (repulsiveforce) {
		// XF contribution
		total_repulsive_stressXF.reset();
		bool rstress_per_particle = false;
 		if (p.out_particle_stress.find('r') != string::npos || couette_stress){
			rstress_per_particle = true;
			for (int i=0; i<np; i++) {
				repulsivestressXF[i] = 0;
			}
		}
		for (int k=0; k<nb_interaction; k++) {
			total_repulsive_stressXF += interaction[k].repulsion.getStressXF();
			if (rstress_per_particle) {
				/* NOTE:
					As the repulsive force is not a contact force, there is an ambiguity defining the stress per particle. Here we make the choice of attributing 1/2 of the interaction stress to each particle.
				*/
				StressTensor sc = 0.5*interaction[k].repulsion.getStressXF();
				unsigned int i, j;
				std::tie(i, j) = interaction[k].get_par_num();
				repulsivestressXF[i] += sc;
				repulsivestressXF[j] += sc;
			}
		}
		total_repulsive_stressXF /= system_volume;
		// GU contribution
		total_repulsive_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_repulsive_stressGU += repulsivestressGU[i];
		}
		total_repulsive_stressGU /= system_volume;
	}
	// Stress from Brownian force
	if (brownian) {
		// GU contribution
		total_brownian_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_brownian_stressGU += brownianstressGU[i];
		}
		total_brownian_stressGU /= system_volume;
	}
	if (mobile_fixed) {
		total_hydrofromfixed_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_hydrofromfixed_stressGU += hydrofromfixedstressGU[i];
		}
		total_hydrofromfixed_stressGU /= system_volume;
	}
	total_stress = total_hydro_stress;
	total_stress += total_contact_stressXF;
	total_stress += total_contact_stressGU; // added (Aug 15 2013)
	if (repulsiveforce) {
		total_stress += total_repulsive_stressXF;
		total_stress += total_repulsive_stressGU;
	}
	if (brownian) {
		total_stress += total_brownian_stressGU;
		if (lowPeclet) {
			// take an averaged stress instead of instantaneous
			stress_avg.update(total_stress, get_time());
			total_stress = stress_avg.get();
		}
	}
	if (mobile_fixed) {
		total_stress += total_hydrofromfixed_stressGU;
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
