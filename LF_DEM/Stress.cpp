//
//  Stress.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
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


void
System::calcStressPerParticle()
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
			if (lubrication_model == 1) {
				interaction[k].lubrication.calcXFunctionsStress();
			} else if (lubrication_model == 2) {
				interaction[k].lubrication.calcXYFunctionsStress();
			} else {
				cerr << "lubrication_model = " << lubrication_model << endl;
				cerr << "lubrication_model = 3 is not implemented" << endl;
				exit(1);
			}
			interaction[k].lubrication.addHydroStress(); // R_SE:Einf-R_SU*v
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

void
System::calcStress()
{
	//////////////////////////////////////////////////////////////
	total_hydro_stress.reset();
	for (int i=0; i<np; i++) {
		total_hydro_stress += lubstress[i];
	}
	total_hydro_stress /= system_volume;
	//////////////////////////////////////////////////////////////
	total_contact_stressGU.reset();
	for (int i=0; i<np; i++) {
		total_contact_stressGU += contactstressGU[i];
	}
	total_contact_stressGU /= system_volume;
	//////////////////////////////////////////////////////////////
	total_contact_stressXF_normal.reset();
	total_contact_stressXF_tan.reset();
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			total_contact_stressXF_normal += interaction[k].contact.getContactStressXF_normal();
			total_contact_stressXF_tan += interaction[k].contact.getContactStressXF_tan();
		}
	}
	total_contact_stressXF_normal /= system_volume;
	total_contact_stressXF_tan /= system_volume;
	total_contact_stressXF = total_contact_stressXF_normal + total_contact_stressXF_tan;
	
	//////////////////////////////////////////////////////////////
	if (repulsiveforce) {
		total_repulsive_stressXF.reset();
		for (int k=0; k<nb_interaction; k++) {
			total_repulsive_stressXF += interaction[k].repulsion.getStressXF();
		}
		total_repulsive_stressXF /= system_volume;
		//////////////////////////////////////////////////////////
		total_repulsive_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_repulsive_stressGU += repulsivestressGU[i];
		}
		total_repulsive_stressGU /= system_volume;
	}
	//////////////////////////////////////////////////////////////
	if (brownian) {
		total_brownian_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_brownian_stressGU += brownianstressGU[i];
		}
		total_brownian_stressGU /= system_volume;
	}

	/* NOTE:
	 *
	 * The total stress DID not include the contact GU terms,
	 * because we consider that the relative motion is not expected hard spheres
	 * and artificial in the soft-sphere contact model.
	 * [Aug 15, 2013]
	 * In the contact model, force is divided into two parts (spring and dash-pot).
	 * In physics, the total force is important.
	 * Therefore, both should be included for the stress calculation.
	 *
	 */

	total_stress = total_hydro_stress;
	total_stress += total_contact_stressXF;
	total_stress += total_contact_stressGU; // added (Aug 15 2013)
	if (repulsiveforce) {
		total_repulsive_stress = total_repulsive_stressXF+total_repulsive_stressGU;
		total_stress += total_repulsive_stress;
	}
	if (brownian) {
		total_stress += total_brownian_stressGU;
		if (lowPeclet) { // take an averaged stress instead of instantaneous
			stress_avg->update(total_stress, time);
			//		cout << time << " " << total_stress.getStressXZ() << " ";
			total_stress = stress_avg->get();
			//		cout << total_stress.getStressXZ() << endl;
		}
	}
	einstein_stress = einstein_viscosity*shear_rate; // should we include that in the hydro_stress definition?
}

