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
#include "StressTensor.h"
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
	int np;
	int np3;
	int maxnum_interactionpair;
	BoxSet boxset;
	double lx;
	double ly;
	double lz;
	double lx_half; // =lx/2
	double ly_half; // =ly/2
	double lz_half; // =lz/2
	double system_volume;
	double radius_max;
	double sq_lub_max;
	double shear_strain;
	double kn;
	double kt;
	double lub_max;
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
	StressTensor* lubstress; // G U + M E
	StressTensor* bgfstress; // by particle
	StressTensor* contactstressGU; // by particle
	StressTensor* colloidalstressGU; // by particle
	StressTensor* brownianstress; // by particle
	int brownianstress_calc_nb;
	StressTensor total_hydro_stress;
	StressTensor total_contact_stressXF_normal;
	StressTensor total_contact_stressXF_tan;
	StressTensor total_contact_stressGU;
	StressTensor total_colloidal_stressXF;
	StressTensor total_colloidal_stressGU;
	StressTensor total_brownian_stress;
	
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
	double total_energy;
	
	void setSystemVolume();
	void setConfiguration(const vector <vec3d> &initial_positions,
						  const vector <double> &radii,
						  double lx_, double ly_, double lz_,
						  double volume_fraction_);
	void setupSystemForGenerateInit();
	void setupSystem();
	void allocatePositionRadius();
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
	void calcTotalPotentialEnergy();
	void setupShearFlow(bool activate){
		if (activate) {
			vel_difference = lz;
		} else {
			vel_difference = 0;
		}
	}
	/*************************************************************/
	inline void setBoxSize(double lx_, double ly_, double lz_){
		lx = lx_;
		lx_half = 0.5*lx;
		ly = ly_;
		ly_half = 0.5*ly;
		lz = lz_;
		lz_half = 0.5*lz;
		cerr << "box: " << lx << ' ' << ly << ' ' << lz << endl;
	}

	inline void Radius_max(double _radius_max){
		radius_max = _radius_max;
	}

	inline double System_volume(){
		return system_volume;
	}
	
	double getParticleContactNumber(){
		return (double)2*contact_nb/np;
	}
	
	inline double Lx(){
		return lx;
	}

	inline double Ly(){
		return ly;
	}

	inline double Lz(){
		return lz;
	}

	inline double Lx_half(){
		return lx_half;
	}

	inline double Ly_half(){
		return ly_half;
	}

	inline double Lz_half(){
		return lz_half;
	}

	inline void Np(int val){
		np = val;
		np3 = 3*np;
	}

	inline int Np(){
		return np;
	}

	inline double Shear_strain(){
		return shear_strain;
	}
	
	inline void Kn(double val){
		kn = val;
	}

	inline double Kn(){
		return kn;
	}

	inline void Kt(double val){
		kt = val;
	}
	
	inline double Kt(){
		return kt;
	}
	
	inline void Lub_max(double val){
		lub_max = val;
	}
	
	inline double Lub_max(){
		return 	lub_max;
	}
	
};
#endif /* defined(__LF_DEM__System__) */
