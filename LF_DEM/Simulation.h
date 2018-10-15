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
//  Copyright (c) 2012-2016 Ryohei Seto and Romain Mari. All rights reserved.
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
#include <set>
#include <algorithm>
#include "global.h"
#include "System.h"
#include "ParameterSet.h"
#include "ParameterSetFactory.h"
#include "DimensionalQty.h"
#include "OutputData.h"
#include "Events.h"
#include "Timer.h"

class Simulation
{
private:
	System sys;
	Parameters::ParameterSet p_initial;
	std::string header_imported_configulation[2];
	Parameters::ControlVariable control_var;
	/*
	 * Resultant data
	 */
	Dimensional::Unit output_unit;
	Dimensional::UnitSystem system_of_units;

	double target_stress_input;
	double input_rate;
	double dimensionless_rate;
	double viscosity;
	double normal_stress_diff1;
	double normal_stress_diff2;
	bool restart_from_chkp;
	double jamming_strain; // Strain stroke to get jamming (simulation_mode = 2)
	time_t time_strain_0;
	time_t time_strain_1;
	time_t time_strain_end;
	int timestep_1;
	int timestep_end;
	Sym2Tensor stress_basis_0;
	Sym2Tensor stress_basis_3;
	Sym2Tensor Einf_base; // This original Einf is kept and can be used when flow is stopped.
	vec3d Omegainf_base; // Same as Einf_base.
	/*
	 * For output data.
	 */
	std::ofstream fout_time;
	std::ofstream fout_input;
	OutputData outdata;
	OutputData outdata_st;
	OutputData outdata_pst;
	OutputData outdata_par;
	OutputData outdata_int;
	std::ofstream fout_boxing;
	/*
	 * For inputs
	 */
	void setupOptionalSimulation(std::string indent);
	std::vector<Sym2Tensor> getParticleStressGroup(std::string group);
	
public:
	/* For DEMsystem*/
	Simulation(State::BasicCheckpoint chkp = State::zero_time_basicchkp);
	void simulationSteadyShear(std::string in_args,
							   std::vector<std::string>& input_files,
							   bool binary_conf,
							   Parameters::ControlVariable control_variable,
							   Dimensional::DimensionalQty<double> control_value,
							   std::string flow_type,
							   std::string simu_identifier);


	void setupSimulation(std::string in_args,
						 std::vector<std::string>& input_files,
						 bool binary_conf,
						 Dimensional::DimensionalQty<double> control_value,
						 std::string flow_type,
						 std::string simu_identifier);
	void setupFlow(Dimensional::DimensionalQty<double> control_value);
	void setConfigToSystem(bool binary_conf, const std::string &filename);
	TimeKeeper initTimeKeeper();
	Parameters::ParameterSet p;
	bool keepRunning();
	// void timeEvolution(double& next_output_data);
	void generateOutput(const std::set<std::string> &output_events, int& binconf_counter);
	/*********** Events  ************/
	std::list <Event> events;
	void setupEvents();
	void handleEvents();
	System &getSys()
	{
		return sys;
	}

	void assertParameterCompatibility();
	void openOutputFiles();
	std::string prepareSimulationName(bool binary_conf,
									  const std::string& filename_import_positions,
									  const std::string& filename_parameters,
									  const std::string& simu_identifier,
									  Dimensional::DimensionalQty<double> control_value);
	void echoInputFiles(std::string in_args,
						std::vector<std::string>& input_files);
	void contactForceParameter(std::string filename);
	void contactForceParameterBrownian(std::string filename);
	void setupNonDimensionalization(Dimensional::DimensionalQty<double> control_value, 
									Parameters::ParameterSetFactory &PFact);
	void stopShearing(TimeKeeper &tk); //simulation mode 22
	void stressReversal(std::set<std::string> &output_events); //simulation mode 2
	/*
	 * For outputs
	 */
	void createDataHeader(std::stringstream& data_header);
	void getSnapshotHeader(std::stringstream& snapshot_header);
	void outputData();
	void outputConfigurationData();
	void outputFinalConfiguration(const std::string&);
	void outputIntFileTxt();
	void outputParFileTxt();
	void outputPstFileTxt();
	void outputConfigurationBinary(std::string);
	void checkpoint();
	vec3d shiftUpCoordinate(double x, double y, double z);
	void relativePositionView(std::vector<vec3d> &pos, std::vector<vec3d> &vel);
	void outputComputationTime();
	bool kill;
	bool force_to_run;
	bool diminish_output;
	/*********** Events  ************/
	void handleEventsShearJamming();
	void handleEventsFragility();
	std::string gitVersion();
	std::string simu_name;
	void timeEvolutionUntilNextOutput(const TimeKeeper &tk);
	void printProgress();
};
#endif /* defined(__LF_DEM__Simulation__) */
