//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__System__
#define __LF_DEM__System__
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <queue>
#include <list>
#include <string>
#include "Stresslet.h"
#include "Interaction.h"
#include "vec3d.h"
#include "BrownianForce.h"
#include "BoxSet.h"
#include "StokesSolver.h"
#include "cholmod.h"
#include "MersenneTwister.h"
using namespace std;

class Simulation;
class Interaction;
class BrownianForce;
class BoxSet;

class System{
private:
	int _np;
	int np3;
	int maxnum_interactionpair;
	BoxSet boxset;
	double _lx;
	double _ly;
	double _lz;
	double _lx_half; // =lx/2
	double _ly_half; // =ly/2
	double _lz_half; // =lz/2
	double system_volume;
	double radius_max;
	double sq_lub_max;
	double shear_strain;
	int linalg_size;
	int linalg_size_per_particle;
	int dof;
	int max_lub_int;
	void timeEvolutionBrownian();
	void timeEvolutionEulersMethod();
	void timeEvolutionPredictorCorrectorMethod();
	void deltaTimeEvolution();
	void deltaTimeEvolutionCorrector();
	void deltaTimeEvolutionPredictor();

	void setContactForceToParticle();
	void setColloidalForceToParticle();
	void buildLubricationTerms(bool rhs=true);
	void buildLubricationRHS();
	void buildContactTerms();
	void buildColloidalForceTerms();

	void addStokesDrag();
	void updateResistanceMatrix();
	void print_res();
	void calcStressesHydroContactBrownian();
	double *lub_cont_forces_init;
	void calcStressesHydroContact();
	double evaluateMaxOverlap();
	double evaluateMaxDispTan();
	void evaluateMaxContactVelocity();
	double evaluateMaxVelocity();
	double evaluateMaxAngVelocity();
	
protected:
public:
	System();
	~System();
	double *v_hydro;
	double *v_cont;
	double *v_colloidal;
	double *v_lub_cont;
	double *v_lub_cont_mid;
	double *v_Brownian_init;
	double *v_Brownian_mid;
	bool in_predictor;
	bool in_corrector;
	int ts; // time steps
	int dimension;
	vec3d *position;
	Interaction *interaction;
	double *radius;
	double *radius_cubic;
	double *angle; // for 2D visualization
	vec3d *velocity;
	vec3d *velocity_predictor;
	vec3d *ang_velocity;
	vec3d *ang_velocity_predictor;
	vec3d *contact_force;
	vec3d *contact_torque;
	vec3d *colloidal_force;
	stresslet* lubstress; // G U + M E
	stresslet* bgfstress; // by particle
	stresslet* contactstressGU; // by particle
	stresslet* colloidalstressGU; // by particle
	stresslet* brownianstress; // by particle
	int brownianstress_calc_nb;
	stresslet total_hydro_stress;
	stresslet total_contact_stressXF_normal;
	stresslet total_contact_stressXF_tan;
	stresslet total_contact_stressGU;
	stresslet total_colloidal_stressXF;
	stresslet total_colloidal_stressGU;
	stresslet total_brownian_stress;
	double kn;
	double kt;
	double lub_max;
	double lub_coeff_contact;
	double mu_static; // static friction coefficient.
	bool friction;
	bool colloidalforce;
	bool brownian;
	int integration_method; // 0: Euler's method 1: PredictorCorrectorMethod
	double diag_stokes_drag;
	double bgf_factor;

	
	int num_interaction;
	double d_strain;
	/*
	 * Leading term of lubrication force is 1/gap_nondim, 
	 * with gap_nondim the gap
	 * gap_nondim = 2r/(a0+a1) - 2.
	 * we set a cutoff for the lubrication interaction,
	 * such that the lub term is proportional to:
	 * 
	 * 1/(gap_nondim+lub_reduce_parameter) 
	 * when gap_nondim > 0.
	 */
	double lub_reduce_parameter;
	double contact_relaxzation_time;
	BrownianForce *fb;
	double shear_disp;
	/* For non-Brownian suspension:
	 * dimensionless_shear_rate = 6*pi*mu*a^2*shear_rate/F_col(0)
	 * For Brownian suspension, it should be Peclet number
	 */
	double dimensionless_shear_rate;
	/* Colloidal force to stabilize suspension
	 * (This gives simple shear-rate depenedence.)
	 */
	double colloidalforce_amplitude; // colloidal force dimensionless
	double colloidalforce_length; // colloidal force length (dimensionless)
	double kb_T;
	double volume_fraction;
	double vel_difference;
	double dt;
	double max_velocity;
	double max_ang_velocity;
	double min_gap_nondim;
	double max_overlap; // = ro-r
	double max_disp_tan;
	
	double overlap_target;
	double disp_tan_target;
	queue<int> deactivated_interaction;

	double max_contact_velo_tan;
	double max_contact_velo_normal;
	double ave_overlap;
	int contact_nb;
	int cnt_monitored_data;
	double average_Fc_normal_norm;
	double max_Fc_normal_norm;
	bool draw_rotation_2d;
	string simu_name;
	ofstream fout_int_data;
//	
	double total_energy;
	
	
	void setSystemVolume();
	void setupSystemForGenerateInit();
	void setupSystem(const vector <vec3d> &initial_positions,
					 const vector <double> &radii);
	void allocateRessources();
	void timeEvolution(double strain_interval);
	void timeEvolutionRelax(int time_step);
	void displacement(int i, const vec3d &dr);

	void checkNewInteraction();
	void checkInteractionEnd();
	void updateInteractions();
	double sq_distance(int i, int j);
	double distance(int i, int j);
	double lubricationForceFactor(int i, int j);
	int periodize(vec3d &);
	void periodize_diff(vec3d &);
	void periodize_diff(vec3d &, int &);
	void updateVelocityLubrication();
	void updateVelocityRestingFluid();
	void forceReset();
	void torqueReset();
	void stressReset();
	void stressBrownianReset();
	void calcStress();
	void analyzeState();
	void computeBrownianStress();
	StokesSolver stokes_solver;
	void lubricationStress(int i, int j);
	void initializeBoxing();
	void calcLubricationForce(); // for visualization of force chains
	set <Interaction*> *interaction_list;
	set <int> *interaction_partners;
	void openFileInteractionData();
	void adjustContactModelParameters(int nb_average);
//	void adjustTimeStep();
	void calcTotalPotentialEnergy();

	void setupShearFlow(bool activate){
		if (activate) {
			vel_difference = _lz;
		} else {
			vel_difference = 0;
		}
	}
	/*************************************************************/
	inline void lx(double length){
		_lx = length;
		_lx_half = 0.5*_lx;
	}

	inline void ly(double length){
		_ly = length;
		_ly_half = 0.5*_ly;
	}

	inline void lz(double length){
		_lz = length;
		_lz_half = 0.5*_lz;
	}

	inline void setRadiusMax(double _radius_max){
		radius_max = _radius_max;
	}

	inline double valSystemVolume(){
		return system_volume;
	}
	
	double getParticleContactNumber(){
		return (double)2*contact_nb/_np;
	}
	
	inline double lx(){
		return _lx;
	}

	inline double ly(){
		return _ly;
	}

	inline double lz(){
		return _lz;
	}

	inline double lx_half(){
		return _lx_half;
	}

	inline double ly_half(){
		return _ly_half;
	}

	inline double lz_half(){
		return _lz_half;
	}

	inline void np(int val){
		_np = val;
		np3 = 3*_np;
	}

	inline int np(){
		return _np;
	}

	inline double strain(){
		return shear_strain;
	}
	
};
#endif /* defined(__LF_DEM__System__) */
