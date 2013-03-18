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
#include <string>
#include "common.h"
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
	queue<int> deactivated_interaction;
	double _lx;
	double _ly;
	double _lz;
	double _lx2; // =lx/2
	double _ly2; // =ly/2
	double _lz2; // =lz/2
	double system_volume;
	double radius_max;
	double sq_lub_max;
	double _rho; // N/V
	double shear_strain;
	double _time;
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
	void displacement(int i, const vec3d &dr);
	void setContactForceToParticle();
	void buildLubricationTerms(bool rhs=true);
	void buildLubricationRHS();
	void buildContactTerms();
	void addStokesDrag();
	void updateResistanceMatrix();
	void print_res();
	void calcStressesHydroContactBrownian();
	double *lub_cont_forces_init;
	void calcStressesHydroContact();
	inline void increment_time(double _dt){
		_time += _dt;
	}
protected:
public:
	double *v_hydro;
	double *v_cont;
	double *v_lub_cont;
	double *v_lub_cont_mid;
	double *v_Brownian_init;
	double *v_Brownian_mid;
	bool in_predictor;
    /* For DEMsystem
     */
	System(){};
	~System();
	int ts; // time steps
	int dimension;
	vec3d *position;
	double *radius;
	double *radius_cubic;
	double *angle; // for 2D visualization
	vec3d *velocity;
	vec3d *velocity_predictor;
	vec3d *ang_velocity;
	vec3d *ang_velocity_predictor;
	vec3d *contact_force;
	vec3d *contact_torque;
	stresslet* lubstress; // G U + M E
	stresslet* bgfstress;
	stresslet* contactstressXF;
	stresslet* contactstressGU;
	stresslet* brownianstress;
	int brownianstress_calc_nb;
	double total_hydro_stress[5];
	double total_contact_stressXF[5];
	double total_contact_stressGU[5];
	double total_brownian_stress[5];
	double total_hydro_stress2[5];
	double total_contact_stress2[5];
	double kn;
	double kt;
	double lub_max;
	double mu_static; // static friction coefficient.
	bool lubrication;
	bool friction;
	bool brownian;
	int integration_method; // 0: Euler's method 1: PredictorCorrectorMethod
	double diag_stokes_drag;
	double bgf_factor;
	bool shearrate_scale_Fc_normal;
	bool poly;
	Interaction *interaction;
	int num_interaction;
	/*
	 * Leading term of lubrication force is 1/ksi, with ksi the gap
	 * ksi = 2r/(a0+a1) - 2.
	 * we set a cutoff for the lubrication interaction,
	 * such that the lub term is proportional to:
	 * 
	 * 1/ksi when ksi > gap_cutoff*(a0+a1)/2. = ksi_cutoff
	 * 1/ksi_cutoff when h <= ksi_cutoff
	 */
	double lub_reduce_parameter;
	BrownianForce *fb;
	double shear_disp;
	double shear_rate;
	double kb_T;
	double volume_fraction;
	double vel_difference;
	double dt;
	double dt_mid;
	double dt_ratio;
	double minvalue_gap_nondim;
	double ave_overlap;
	double average_contact_time;
	int contact_nb;
	double average_nearing_time;
	int nearing_nb;
	bool draw_rotation_2d;
	string simu_name;
	void setSystemVolume();
	void setupSystem(const vector<vec3d> &initial_positions,
					 const vector <double> &radii);
	void allocateRessources();
	void timeEvolution(int time_step);
	void checkNewInteraction();
	void checkInteractionEnd();
	void updateInteractions(bool _in_predictor=false);
	double sq_distance(int i, int j);
	double distance(int i, int j);
	double lubricationForceFactor(int i, int j);
	int periodize(vec3d &);
	void periodize_diff(vec3d &);
	void periodize_diff(vec3d &, int &);
	void updateVelocity();
	void updateVelocityLubrication();
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
	ofstream fout_trajectory;
	
	/*************************************************************/
	inline void lx(double length){
		_lx=length;
		_lx2=0.5*_lx;
	}
	inline void ly(double length){
		_ly=length;
		_ly2=0.5*_ly;
	}
	inline void lz(double length){
		_lz=length;
		_lz2=0.5*_lz;
	}
	inline void setRadiusMax(double _radius_max){
		radius_max = _radius_max;
	}
	inline double valSystemVolume(){
		return system_volume;
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
	inline double lx2(){
		return _lx2;
	}
	inline double ly2(){
		return _ly2;
	}
	inline double lz2(){
		return _lz2;
	}
	inline void np(int val){
		_np = val;
		np3=3*_np;
	}
	inline int np(){
		return _np;
	}
	inline double rho(){ // N/V
		return _rho;
	}
	inline double time(){
		return _time;
	}
	inline double strain(){
		return shear_strain;
	}
};
#endif /* defined(__LF_DEM__System__) */
