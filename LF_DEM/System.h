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
#include "BoxSet.h"
#include "StokesSolver.h"
#include "cholmod.h"
#include "MersenneTwister.h"
using namespace std;

class Simulation;
class Interaction;
class BoxSet;

class System{
private:
	int np;
	int maxnb_interactionpair;
	int maxnb_interactionpair_per_particle;
	int nb_of_active_interactions;
	int ts; // time steps
	double time;
	double dt; // <=== It should be called d_strain.
	double disp_max;
	double lx;
	double ly;
	double lz;
	double lx_half; // =lx/2
	double ly_half; // =ly/2
	double lz_half; // =lz/2
	/* l_periodic_threshold = lz - a1*lub_max (a1 is larger particle)
	 * threshold distance which can intaract each other as periodic image.
	 */
	double system_volume;
	double volume_fraction;
	double sq_lub_max;
	double shear_strain;
	double lub_max;
	double sd_coeff;
	double mu_static; // static friction coefficient.
	double kb_T; // dimensionless kb_T = 1/Pe
	int linalg_size;
	int linalg_size_per_particle;
	int dof;
	int max_lub_int;
	double repulsiveforce_length; // repulsive force length (dimensionless)
	int integration_method; // 0: Euler's method 1: PredictorCorrectorMethod
	/* data */
	void (System::*timeEvolutionDt)(bool);
	void timeEvolutionEulersMethod(bool calc_stress);
	void timeEvolutionPredictorCorrectorMethod(bool calc_stress);
	void timeStepMove();
	void timeStepMoveCorrector();
	void timeStepMovePredictor();
	void timeStepBoxing();
	void setContactForceToParticle();
	void setRepulsiveForceToParticle();
	void buildHydroTerms(bool, bool);
	void (System::*buildLubricationTerms)(bool, bool);
	void buildLubricationTerms_squeeze(bool mat, bool rhs); // lubrication_model = 1
	void buildLubricationTerms_squeeze_tangential(bool mat, bool rhs); // lubrication_model = 2
	//void buildBrownianTerms();
	void generateBrownianForces();
	void buildContactTerms(bool);
	void buildRepulsiveForceTerms(bool);
	void print_res();
	double evaluateMinGap();
	double evaluateMaxDispTan();
	double evaluateMaxFcNormal();
	double evaluateMaxFcTangential();
	void evaluateMaxContactVelocity();
	double evaluateMaxVelocity();
	double evaluateMaxAngVelocity();
	void countNumberOfContact();
	MTRand *r_gen;
	double *radius_cubed;
	bool strain_controlled;
	bool stress_controlled;

protected:
public:
	System();
	~System();
	bool brownian;	
	bool in_predictor;
	bool twodimension;
	/* zero_shear:
	 * To be used for relaxation to generate initial configuration.
	 */
	bool zero_shear;
	/* When Pe is small, the simulation parameters need to be adjusted differently.
	 * For high Pe, spring constants are scaled with shear rate being proportional to Pe.
	 * In dimensionless simulation, kn and kt appear the same.
	 * For low Pe, they need to be scaled with Pe.
	 *
	 */
	double Pe_switch;
	double shear_strain_end;
	double critical_normal_force;
	double scale_factor_SmallPe;
	vec3d *position;
	Interaction *interaction;
	BoxSet boxset;
	double *radius;
	double *angle; // for 2D visualization
	double *resistance_matrix_dblock;
	vec3d *velocity;
	vec3d *velocity_predictor;
	vec3d *na_velocity;
	vec3d *ang_velocity;
	vec3d *ang_velocity_predictor;
	vec3d *na_ang_velocity;
	vec3d *vel_repulsive;
	vec3d *ang_vel_repulsive;
	vec3d *vel_contact;
	vec3d *ang_vel_contact;
	vec3d *vel_hydro;
	vec3d *ang_vel_hydro;
	vec3d *vel_brownian;
	vec3d *ang_vel_brownian;
	vec3d *contact_force;
	vec3d *contact_torque;
	vec3d *repulsive_force;
	double *brownian_force;
	StressTensor* lubstress; // G U + M E
	StressTensor* contactstressGU; // by particle
	StressTensor* repulsivestressGU; // by particle
	StressTensor* brownianstressGU; // by particle
	StressTensor* brownianstressGU_predictor; // by particle
	StressTensor contactstressXF_normal;
	StressTensor contactstressXF_tan;
	StressTensor total_hydro_stress;
	StressTensor total_contact_stressXF_normal;
	StressTensor total_contact_stressXF_tan;
	StressTensor total_contact_stressGU;
	StressTensor total_repulsive_stressXF;
	StressTensor total_repulsive_stressGU;
	StressTensor total_brownian_stressGU;
	double dt_max;
	double kn;
	double kt;
	double kr;
	double dt_lowPeclet;
	double kn_lowPeclet;
	double kt_lowPeclet;
	double kr_lowPeclet;
	int friction_model;
	bool friction;
	bool rolling_friction;
	bool repulsiveforce;
	double lub_coeff_contact;
	double cohesive_force;
	// resistance coeffient for normal mode
	double log_lub_coeff_contact_tan_dashpot;
	double log_lub_coeff_contact_tan_lubrication;
	double log_lub_coeff_contact_tan_total;
	set <Interaction*> *interaction_list;
	set <int> *interaction_partners;
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication
	 * 3 1/xi lubrication for h>0 and tangential dashpot.
	 */
	int lubrication_model;
	int nb_interaction;
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
	double contact_relaxation_time;
	double contact_relaxation_time_tan;
	double shear_disp;
	/* For non-Brownian suspension:
	 * dimensionless_shear_rate = 6*pi*mu*a^2*shear_rate/F_repulsive(0)
	 * For Brownian suspension, it should be Peclet number
	 */
	double dimensionless_shear_rate;
	double repulsiveforce_amplitude; // dimensionless
	/* Velocity difference between top and bottom
	 * in Lees-Edwards boundary condition
	 * vel_difference = shear_rate * lz
	 */
	double vel_difference;
	double max_velocity;
	double max_relative_velocity;
	double max_sliding_velocity;
	double max_ang_velocity;
	double min_gap_nondim;
	double max_disp_tan;
	double overlap_target;
	double disp_tan_target;
	double max_kn;
	queue<int> deactivated_interaction;
	double max_contact_velo_tan;
	double max_contact_velo_normal;
	double ave_contact_velo_tan;
	double ave_contact_velo_normal;
	double ave_sliding_velocity;
	int contact_nb; // gap < 0
	int fric_contact_nb; // fn > f* in the critical load model
	double average_fc_normal;
	double max_fc_normal;
	double max_fc_tan;
	double strain_interval_output_data;
	double strain_interval_output;
	string simu_name;
	bool kn_kt_adjustment;
	double target_stress_input;
	double target_stress;
	void setSystemVolume(double depth = 0);
	void setConfiguration(const vector <vec3d> &initial_positions,
						  const vector <double> &radii,
						  double lx_, double ly_, double lz_);
	void setInteractions_GenerateInitConfig();
	void setupSystem(string control);
	void setupBrownian();
	void allocatePositionRadius();
	void allocateRessources();
	void timeEvolution(double strain_interval);
	void displacement(int i, const vec3d &dr);
	void checkNewInteraction();
	void checkInteractionEnd();
	void updateInteractions();
	double lubricationForceFactor(int i, int j);
	int periodize(vec3d &);
	void periodize_diff(vec3d &, int &);
	void computeVelocities(bool divided_velocities);
	void forceReset();
	void torqueReset();
	void stressReset();
	void stressBrownianReset();
	void calcStress();
	void calcStressPerParticle();
	void analyzeState();
	StokesSolver stokes_solver;
	void lubricationStress(int i, int j);
	void initializeBoxing();
	void calcLubricationForce(); // for visualization of force chains
	int adjustContactModelParameters();
	void setupShearFlow(bool activate){
		if (activate) {
			/* In dimensionless simulations for non Browninan simulation,
			 * shear rate is always 1.
			 */
			vel_difference = lz;
		} else {
			vel_difference = 0;
		}
	}
	/*************************************************************/
	void setBoxSize(double lx_, double ly_, double lz_){
		lx = lx_;
		lx_half = 0.5*lx;
		ly = ly_;
		ly_half = 0.5*ly;
		lz = lz_;
		lz_half = 0.5*lz;
		cerr << "box: " << lx << ' ' << ly << ' ' << lz << endl;
	}
	//	inline double System_volume(){
	//		return system_volume;
	//	}
	double getParticleContactNumber(){
		return (double)2*contact_nb/np;
	}
	void set_integration_method(int val){integration_method = val;}
	void set_lubrication_model(int val){lubrication_model = val;}
	double get_lx(){return lx;}
	double get_ly(){return ly;}
	double get_time(){return time;}
	inline double get_lz(){return lz;}
	inline double Lx_half(){return lx_half;}
	inline double Ly_half(){return ly_half;}
	inline double Lz_half(){return lz_half;}
	inline void set_np(int val){np = val;}
	inline int get_np(){return np;}
	inline double get_shear_strain(){return shear_strain;}
	inline double get_kn(){return kn;}
	inline void set_kt(double val){kt = val;}
	inline double get_kt(){return kt;}
	inline void set_kr(double val){kr = val;}
	inline double get_kr(){return kr;}
	inline void set_lub_max(double val){lub_max = val;}
	inline double get_lub_max(){return lub_max;}
	inline void set_dt(double val){dt = val;}
	void set_disp_max(double val){disp_max = val;}
	inline double get_dt(){return dt;}
	/* inline void set_repulsiveforce_amplitude(double val){ */
	/* 	repulsiveforce_amplitude = val;} */
	inline double get_repulsiveforce_amplitude(){return repulsiveforce_amplitude;}
	void set_repulsiveforce_length(double val){repulsiveforce_length = val;}
	void set_sd_coeff(double val){sd_coeff = val;}
	inline double get_repulsiveforce_length(){return repulsiveforce_length;}
	void set_mu_static(double val){mu_static = val;}
	inline double get_mu_static(){return mu_static;}
	//inline double get_lub_coeff_contact(){return lub_coeff_contact;}
	//	inline double get_log_lub_coeff_dynamicfriction(){
	/* In a sliding state, the resistance coeffient is the sum of
	 * lubrication and dashpot.
	 * In our standard model, we do not set the dashpot.
	 * Only tangential lubrications are considered.
	 */
	//		return log_lub_coeff_contact_tan_total;
	//	}
	inline double get_nb_of_active_interactions(){return nb_of_active_interactions;}

};
#endif /* defined(__LF_DEM__System__) */
