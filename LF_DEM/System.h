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
//#include <Accelerate/Accelerate.h>
#include "Interaction.h"

#include "vec3d.h"
#include "BrownianForce.h"
#include "BoxSet.h"
#include "StokesSolver.h"

using namespace std;

class Simulation;
class Interaction;
class BrownianForce;
class BoxSet;



class System{
private:
	int np3;
	int maxnum_interactionpair;
	queue<int> deactivated_interaction;

	double _lx;
	double _ly;
	double _lz;
	double _lx2; // =lx/2
	double _ly2; // =ly/2
	double _lz2; // =lz/2
	double system_volume;
	double radius_max;

	void buildLubricationTerms(bool);
	void buildBrownianTerms();
	void buildContactTerms();
	void addStokesDrag();
	void factorizeResistanceMatrix();

	void updateResistanceMatrix();

	int linalg_size;
	int linalg_size_per_particle;
	int dof;
	int max_lub_int;

	// rhs vector building
	double *v_lub_cont;
	double *v_Brownian_init;
	double *v_Brownian_mid;

	BoxSet* boxset;
	void print_res();

	
protected:
public:
    /* For DEMsystem
     */
	System(){};
	~System();
	int np; // number of particles
	int ts; // time steps
	int dimension;
	vec3d *position;
	double *radius;
	double *angle; // for 2D visualization
	vec3d *velocity;
	vec3d *relative_velocity;
	vec3d *relative_velocity_lub_cont;
	vec3d *relative_velocity_brownian;
	vec3d *ang_velocity;
	vec3d *total_force;
	vec3d *lubrication_force;
	vec3d *contact_force;
	vec3d *brownian_force;
	vec3d *total_velocity;
	vec3d *lubrication_velocity;
	vec3d *contact_velocity;
	vec3d *brownian_velocity;
	vec3d *torque; // right now only contact torque
	vec3d *lub_force; // Only for outputing data
	vector <stresslet> lubstress; // G U + M E
	vector <stresslet> lubstress2; // r * F_lub
	vector <stresslet> contactstress;
	vector <stresslet> brownianstress;
	int brownianstress_calc_nb;
	double total_stress_bgf;
	double total_lub_stress[5];
	double total_contact_stress[5];
	double total_brownian_stress[5];
	double kn;
	double kt;
	double eta;
	double lub_max;
	double sq_lub_max;
	double mu_static; // static friction coefficient.
	double mu_dynamic;// dynamic friction coefficient.
	double dynamic_friction_critical_velocity;
	double sq_critical_velocity;
	bool lubrication;
	bool friction;
	bool brownian;
	double diag_stokes_drag;
	double bgf_factor;
	bool poly;
	Interaction *interaction;
	int num_interaction;
	double shear_strain;
	
	/*
	 * Leading term of lubrication force is 1/ksi, with ksi the gap
	 * ksi = 2r/(a0+a1) - 2.
	 * we set a cutoff for the lubrication interaction,
	 * such that the lub term is proportional to:
	 * 
	 * 1/ksi when ksi > gap_cutoff*(a0+a1)/2. = ksi_cutoff
	 * 1/ksi_cutoff when h <= ksi_cutoff
	 */
	double gap_cutoff;
	BrownianForce *fb;
	double shear_disp;
	double shear_rate;
	double kb_T;
	double volume_fraction;
	double vel_difference;
	double dt;
	double dt_mid;
	double dt_ratio;
	double gap_min;
	double ave_overlap;
	bool draw_rotation_2d;
	vector <int> lubparticle;
	vector <double> lubparticle_vec[3];
	string simu_name;
	/* The definition of contact is
	 * when the gap is smaller than "dist_near"
	 */
	int cnt_nearing_number[10];
	double max_nearing_time;
	double ave_nearing_time;
	double dist_near;
	bool near;
	int num_nearing;
	vector<double> nearing_time_record;
	
	
	/*
	 * nearing_number[i] means
	 * the number of nearing of particle i.
	 */
	vector<int> nearing_number;
		
	/*************************************************************/
	void lx(double length){
		_lx=length;
		_lx2=0.5*_lx;
	}
	void ly(double length){
		_ly=length;
		_ly2=0.5*_ly;
	}
	void lz(double length){
		_lz=length;
		_lz2=0.5*_lz;
	}
	void setRadiusMax(double _radius_max){
		radius_max = _radius_max;
	}
	double valSystemVolume(){
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

	void set_np(int _np){
		np=_np;
		np3=3*np;
	}
	//	void prepareSimulation();
	void setSystemVolume();
	void setupSystem(const vector<vec3d> &initial_positions,
					 const vector <double> &radii);
	void allocateRessources();
	void timeEvolution(int time_step);
	void checkNewInteraction();
	void checkInteractionEnd();
	void updateInteractions();
	void calcContactForces();
	double sq_distance(int i, int j);
	double distance(int i, int j);
	double lubricationForceFactor(int i, int j);
	void displacement(int i, const double &dx_, const double &dy, const double &dz);
	void periodize(vec3d &);
	void periodize_diff(vec3d &);
	void periodize_diff(vec3d &, int &);
	void updateVelocity();
	void updateVelocityLubrication();
	void updateVelocityLubricationBrownian();
	void deltaTimeEvolution();
	void forceReset();
	void torqueReset();
	void stressReset();
	void stressBrownianReset();
	void calcStress();

	void analyzeState();
	void computeBrownianStress();


	int numpart(){
		return np;
	}

	StokesSolver *stokes_solver;

	void lubricationStress(int i, int j);
	void initializeBoxing();
	void calcLubricationForce(); // for visualization of force chains
	// interactions
	bool out_pairtrajectory;
	set <Interaction*> *interaction_list;
	set <int> *interaction_partners;
	ofstream fout_trajectory;


};
#endif /* defined(__LF_DEM__System__) */
