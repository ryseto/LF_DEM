/**
 \class Simulation
 \brief Class launching the simulation by setting up the System class and performing predefined shear evolution
 \author Ryohei Seto
 \author Romain Mari
 */
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
#include <ctime>
#include <map>
#include <algorithm>
#include "global.h"
#include "System.h"
#include "ParameterSet.h"
#include "InputValue.h"
#include "OutputData.h"
#include "Events.h"

class Simulation
{
private:
	System sys;
	ParameterSet p_initial;
	std::map <std::string, std::string> input_force_units;   // pairs: (force_type, unit)
	std::map <std::string, double> input_force_values;   // pairs: (force_type, value)
	std::map <std::string, double> dimensionless_numbers; // pairs: (force_type_1/force_type_2, force_value_1/force_value_2)
	std::map <std::string, std::string> unit_longname;
	std::list <InputValue> input_values;
	double volume_or_area_fraction;
	std::string header_imported_configulation[2];
	std::string control_var;
	double shear_rate_expectation;
	double time_interval_output_data;
	double time_interval_output_config;
	double strain_interval_output_data;
	double strain_interval_output_config;
	double strain_end;
	double time_end;
	/*
	 * Resultant data
	 */
	std::string internal_unit_scales;
	std::string output_unit_scales;
	double target_stress_input;
	double input_rate;
	std::string input_rate_unit;
	int time_strain_0;
	int time_strain_1;
	int time_strain_end;
	int timestep_1;
	int timestep_end;
	/*
	 * For output data.
	 */
	std::ofstream fout_particle;
	std::ofstream fout_interaction;
	std::ofstream fout_time;
	std::ofstream fout_input;
	OutputData outdata;
	OutputData outdata_st;
	OutputData outdata_pst;
	/*
	 * For inputs
	 */
	void setDefaultParameters();
	void readParameterFile(const std::string &);
	void openOutputFiles(bool, const std::string &, const std::string &, const std::string &);
	void prepareSimulationName(bool, const std::string &, const std::string &, const std::string &);
	void echoInputFiles(std::string in_args, std::vector<std::string> &input_files);
	void autoSetParameters(const std::string &keyword, const std::string &value);
	void contactForceParameter(std::string filename);
	void contactForceParameterBrownian(std::string filename);
	void importPreSimulationData(std::string filename);
	void importConfiguration(const std::string &);
	void importConfigurationBinary(std::ifstream &file_import, const std::string &);
	void importContactsBinary(std::ifstream &file_import);
	void exportForceAmplitudes();
	void setLowPeclet();
	void convertForceValues(std::string new_long_unit);
	void convertInputValues(std::string new_long_unit);
	void resolveUnitSystem(std::string long_unit);
	void setUnitScaleRateControlled();
	void setUnitScaleMagnetic();
	void setupNonDimensionalization(double dimensionlessnumber, std::string input_scale);
	void convertInputForcesRateControlled(double dimensionlessnumber, std::string input_scale);
	void convertInputForcesStressControlled(double dimensionlessnumber, std::string input_scale);
	void convertInputForcesMagnetic(double dimensionlessnumber, std::string rate_unit);
	void catchSuffixedValue(std::string type, std::string keyword, std::string value_str, double *value_ptr);
	/*
	 * For outputs
	 */
	void evaluateData();
	void createDataHeader(std::stringstream &data_header);
	void outputDataHeader(std::ofstream &fout);
	void getSnapshotHeader(std::stringstream &snapshot_header);
	void outputData();
	void outputDataMagnetic();
	void outputConfigurationData();
	void outputFinalConfiguration(const std::string &);
	void outputConfigurationBinary();
	void outputConfigurationBinary(std::string);
	double getRate();
	vec3d shiftUpCoordinate(double x, double y, double z);
	void outputComputationTime();
	bool kill;

	/*********** Events  ************/
	void handleEventsShearJamming();
	void handleEventsFragility();

 public:
	/* For DEMsystem*/
	Simulation();
	~Simulation();
	void simulationSteadyShear(std::string in_args, std::vector<std::string> &input_files, bool binary_conf,
							   double dimensionless_number, std::string input_scale, std::string control_variable);
	// void simulationfinedSequence(std::string seq_type, std::string in_args, std::vector<std::string> &input_files, bool binary_conf, std::string control_variable);

	void simulationInverseYield(std::string in_args,
								std::vector<std::string> &input_files,
								bool binary_conf,
								double dimensionless_number,
								std::string input_scale,
								std::string control_variable);

	void simulationMagnetic(std::string in_args, std::vector<std::string> &input_files,
							bool binary_conf, double dimensionless_number,
							std::string input_scale, std::string control_variable);
	void setupSimulation(std::string in_args,
												 std::vector<std::string> &input_files,
												 bool binary_conf,
												 double dimensionlessnumber,
												 std::string input_scale);
	void setControlVariable(const std::string & var){
		control_var = var;
	};
	ParameterSet p;
	bool keepRunning();
	void timeEvolution(double next_output_data);
	void generateOutput(double &next_output_data, double &next_output_config, int &binconf_counter);
	/*********** Events  ************/
	std::list <Event> events;
	void setupEvents();
	void handleEvents();
	System &getSys(){ return sys;}
};
#endif /* defined(__LF_DEM__Simulation__) */
