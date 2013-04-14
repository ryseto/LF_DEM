//
//  Simulation.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__Simulation__
#define __LF_DEM__Simulation__
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include "System.h"

class Simulation{
private:
	System sys;
	string filename_addition;
	int num_of_particle;
	bool dt_adjustment;
	bool kn_kt_adjustment;
	double shear_strain_end;
	bool import_positions;
	vector<vec3d> initial_positions;
	vector<double> radii;
	string filename_import_positions;
	string filename_parameters;

	double unit_of_force;
	double unit_of_velocity;
	double unit_of_length;
	double viscosity_solvent; // It can be real value.
	double radius_of_particle; // It can be real value.
	/* Colloidal force parameter
	 * A exp(-h/lambda)
	 * cf_amp_dl0 = A/F_0 at shear rate = 1.0
	 */
	double cf_amp_dl0;
	/*
	 * Resultant data
	 */
	double Viscosity;
	double N1;
	double N2;
	double Viscosity_h;
	double N1_h;
	double N2_h;
	double Viscosity_c_XF;
	double N1_c_XF;
	double N2_c_XF;
	double Viscosity_c_GU;
	double N1_c_GU;
	double N2_c_GU;
	double Viscosity_col_XF;
	double Viscosity_col_GU;
	double Viscosity_b;
	double N1_b;
	double N2_b;
  	double N1_2;
	double N2_2;
	/*
	 * For output data.
	 */
	ofstream fout_rheo;
	ofstream fout_particle;
	ofstream fout_interaction;
	bool out_data_particle;
	bool out_data_interaction;
	bool origin_zero_flow;
	void timeEvolution();
	void evaluateData();
	/*
	 * For inputs
	 */
	void setDefaultParameters();
	void readParameterFile();
	void openOutputFiles();
	void prepareSimulationName();
	void autoSetParameters(const string &keyword,
						   const string &value);
	void importInitialPositionFile();
	/*
	 * For outputs
	 */
	void setUnits();
	void output_yap();
	void output_vel();
	void outputDataHeader(ofstream &fout);
	void initContactPair();
	void outputRheologyData();
	void outputConfigurationData();
	vec3d shiftUpCoordinate(double x, double y, double z);
public:
    /* For DEMsystem
     */
	Simulation();
	~Simulation();
	void SimulationMain(int argc, const char * argv[]);
	void RelaxationZeroShear(vector<vec3d> &positions,
							  vector<double> &radii,
							  double lx, double ly, double lz);
};
#endif /* defined(__LF_DEM__Simulation__) */
