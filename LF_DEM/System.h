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
	int maxnb_interactionpair;
	BoxSet boxset;
	int ts; // time steps
	double dt;
	double dt_max;
	double disp_max;
	double lx;
	double ly;
	double lz;
	double lx_half; // =lx/2
	double ly_half; // =ly/2
	double lz_half; // =lz/2
	double system_volume;
	double sq_lub_max;
	double shear_strain;
	double kn;
	double kt;
	double lub_max;
	double lub_coeff_contact; // resistance coeffient for normal mode
	double tang_coeff_contact; // resistance coeffient for tangential mode
	double mu_static; // static friction coefficient.
	double bgf_factor;
	double kb_T;

	int linalg_size;
	int linalg_size_per_particle;
	int dof;
	int max_lub_int;
	double colloidalforce_amplitude; // colloidal force dimensionless
	double colloidalforce_length; // colloidal force length (dimensionless)
	int integration_method; // 0: Euler's method 1: PredictorCorrectorMethod
	bool twodimension;
	/* data */
	int intr_max_fc_normal;
	int intr_max_fc_tan;
	vector<double> max_fc_normal_history; // for kn-kt ajusting algorithm
	vector<double> max_fc_tan_history; // for kn-kt ajusting algorithm
	vector<double> velocity_history;
	/* cnt_parameter_changed is used for kn-kt ajusting algorithm.
	 * Spring constants are changed in the simulation.
	 * This may cause large contact forces.
	 * These values are not considered for the next deteremination.
	 */
	bool after_parameter_changed; //
	void timeEvolutionBrownian();
	void timeEvolutionEulersMethod();
	void timeEvolutionPredictorCorrectorMethod();
	void deltaTimeEvolution();
	void deltaTimeEvolutionRelax();
	void deltaTimeEvolutionCorrector();
	void deltaTimeEvolutionPredictor();

	void setContactForceToParticle();
	void setColloidalForceToParticle();
	void buildLubricationTerms(bool rhs=true);
	
//	void buildLubricationRHS();
	void buildContactTerms();
	void buildColloidalForceTerms();

	void addStokesDrag();
	void updateResistanceMatrix();
	void print_res();
	//void calcStressesHydroContactBrownian();
	double *lub_cont_forces_init;
	void calcStressesHydroContact();


	double evaluateMaxOverlap();
	double evaluateMaxDispTan();
	void evaluateMaxContactVelocity();
	double evaluateMaxVelocity();
	double evaluateMaxAngVelocity();
protected:
public:
	~System();
	double *v_total; // Velocities and angular velocities of all particles.
	double *v_hydro; // Velocities balancing with the 
	double *v_cont;
	double *v_colloidal;

	//	double *v_lub_cont_mid;
	//	double *v_Brownian_init;
	//	double *v_Brownian_mid;
	bool in_predictor;
	bool in_corrector;
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
	bool friction;
	bool colloidalforce;
	bool brownian;
	
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication (only squeeze mode)
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
	double vel_difference;
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
	double average_fc_normal;
	double max_fc_normal;
	double max_fc_tan;
	string simu_name;
	ofstream fout_int_data;
	double total_energy;
	
	void setSystemVolume(double depth = 0);
	void setConfiguration(const vector <vec3d> &initial_positions,
						  const vector <double> &radii,
						  double lx_, double ly_, double lz_);
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
	void adjustContactModelParameters();
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

	inline double System_volume(){
		return system_volume;
	}
	
	double getParticleContactNumber(){
		return (double)2*contact_nb/np;
	}
	void set_integration_method(int val){integration_method = val;}
	void set_lubrication_model(int val){lubrication_model = val;}
	void set_bgf_factor(int val){bgf_factor = val;}
	void set_kb_T(double val){kb_T = val;}
	double get_kb_T(){return kb_T;}
	double get_lx(){return lx;}
	double get_ly(){return ly;}
	inline double get_lz(){return lz;}
	inline double Lx_half(){return lx_half;}
	inline double Ly_half(){return ly_half;}
	inline double Lz_half(){return lz_half;}
	inline void set_np(int val){np = val, np3 = 3*val;}
	inline int get_np(){return np;}
	inline double get_shear_strain(){return shear_strain;}
	inline void set_kn(double val){kn = val;}
	inline double get_kn(){return kn;}
	inline void set_kt(double val){kt = val;}
	inline double get_kt(){return kt;}
	inline void set_lub_max(double val){lub_max = val;}
	inline double get_lub_max(){return lub_max;}
	inline void set_dt(double val){dt = val;}
	void set_dt_max(double val){dt_max = val;}
	void set_disp_max(double val){disp_max = val;}
	inline double get_dt(){return dt;}
	inline void set_colloidalforce_amplitude(double val){colloidalforce_amplitude = val;}
	inline double get_colloidalforce_amplitude(){return colloidalforce_amplitude;}
	inline void set_colloidalforce_length(double val){colloidalforce_length = val;}
	inline double get_colloidalforce_length(){return colloidalforce_length;}
	inline void set_mu_static(double val){mu_static = val;}
	inline double get_mu_static(){return mu_static;}
	inline double get_lub_coeff_contact(){return lub_coeff_contact;}
	inline double get_tang_coeff_contact(){return tang_coeff_contact;}

};
#endif /* defined(__LF_DEM__System__) */
