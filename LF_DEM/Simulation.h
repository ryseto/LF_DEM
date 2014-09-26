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
#include <string>
#include "System.h"

class Simulation{
private:
	System sys;
	vector<double> radius;
	string filename_import_positions;
	string filename_parameters;
	string filename_sequence;
	ostringstream string_control_parameters;
	double strain_interval_knkt_adjustment;
	double volume_fraction;
	string import_line[2];
	string control_var;
	bool user_sequence;
	double shear_rate_expectation;
	double time_interval_output_data;
	double time_interval_output_config;		
	/*
	 * Resultant data
	 */
	double viscosity;
	double normalstress_diff_1;
	double normalstress_diff_2;
	double particle_pressure;
	StressTensor total_stress;
	StressTensor total_contact_stressXF;
	StressTensor total_repulsive_stress;
	double viscosity_hydro; // Only lubrication...
	double normalstress_diff_1_hydro;
	double normalstress_diff_2_hydro;
	double viscosity_cont_XF;
	double normalstress_diff_1_cont_XF;
	double normalstress_diff_2_cont_XF;
	double particle_pressure_cont;
	double viscosity_friction; // Fc_tan contribution.
	double normalstress_diff_1_friction;
	double normalstress_diff_2_friction;
	double viscosity_cont_GU;
	double normalstress_diff_1_cont_GU;
	double normalstress_diff_2_cont_GU;
	double viscosity_repulsive_XF;
	double normalstress_diff_1_repulsive_XF;
	double normalstress_diff_2_repulsive_XF;
	double particle_pressure_repulsive;
	double viscosity_repulsive_GU;
	double normalstress_diff_1_repulsive_GU;
	double normalstress_diff_2_repulsive_GU;
	double viscosity_brownian;
	double normalstress_diff_1_brownian;
	double normalstress_diff_2_brownian;
	/*
	 * For output data.
	 */
	ofstream fout_rheo;
	ofstream fout_particle;
	ofstream fout_interaction;
	ofstream fout_st;
	bool out_data_particle;
	bool out_data_interaction;
	bool origin_zero_flow;
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
	void contactForceParameter(string filename);
	void importPreSimulationData(string filename);
	/*
	 * For outputs
	 */
	void evaluateData();
	void outputDataHeader(ofstream &fout);
	void outputRheologyData();
	void outputStressTensorData();
	void outputConfigurationData();
	void outputFinalConfiguration();
	vec3d shiftUpCoordinate(double x, double y, double z);
	void setupSimulationSteadyShear(vector<string> &input_files, double peclet_num,
						 double scaled_repulsion, double scaled_cohesion, double scaled_critical_load, string control_variable);
public:
	/* For DEMsystem*/
	Simulation();
	~Simulation();
	void simulationSteadyShear(vector<string> &input_files, double peclet_num, double scaled_repulsion,
							   double scaled_cohesion,
							   double scaled_critical_load, string control_variable);
	void simulationUserDefinedSequence(string seq_type, vector<string> &input_files, string control_variable);
};
#endif /* defined(__LF_DEM__Simulation__) */

