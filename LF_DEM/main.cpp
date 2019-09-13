//
//  main.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <iostream>
#include <string>
#include <stdexcept>
#include <csignal>
#include <getopt.h>
#include "global.h"
#include "Simulation.h"
#include "GenerateInitConfig.h"
#ifndef GIT_VERSION
#include "VersionInfo.h"
#endif
#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace std;

#ifdef SIGINT_CATCH
volatile sig_atomic_t sig_caught = 0;
#endif

int mainConventional(int argc, char **argv);
int mainLammpsLike(int argc, char **argv);

std::string prepareSimulationNameFromChkp(const std::string& filename_chkp)
{
	/**
	 \brief Determine simulation name when starting from checkpoint
	 */
	string::size_type chkp_name_start = filename_chkp.rfind("chk_") + 4;
	string::size_type chkp_name_end = filename_chkp.rfind(".dat");

	return filename_chkp.substr(chkp_name_start, chkp_name_end-chkp_name_start);
}

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

int main(int argc, char **argv)
{
	cout << endl << "LF_DEM version " << GIT_VERSION << endl << endl;
	string usage = "(1) Simulation\n $ LF_DEM [-r Rate] [-s Stress] [-R Rate_Sequence] [-S body-force/repulsive force (for sedimentation)]\
	[-e] [-m ?] [-k kn_kt_File] [-v Simulation_Identifier] [-i Provisional_Data] [-n]\
	Configuration_File Parameter_File \
	\n ('-r inf' for infinite Pe and '-r zero' for zero shear simulations) \
	\n\n OR \n\n(2) Generate initial configuration\n $ LF_DEM [-a Random_Seed] [-p Volume_Fraction] [-b cluster radius] -g[c/w/s/b/f]\n";

	int generate_init = 0;
	string type_init_config = "normal";

	int random_seed = 1;
	double volume_frac_gen = 0;
	double cluster_phi = 0;
	bool binary_conf = false;
	bool force_to_run = false;
	string simulation_type = "shear rheology";
	string chkp_filename = "";
	string simu_name;
	
	Dimensional::DimensionalQty<double> control_value;
	Parameters::ControlVariable control_variable = Parameters::ControlVariable::rate;
	string simu_identifier = "";
	vector<string> str_vec;
	map <string, string> input_files;

	const struct option longopts[] = {
		{"rate-controlled",   required_argument, 0, 'r'},
		{"rate-infty",        required_argument, 0, '8'},
		{"stress-controlled", required_argument, 0, 's'},
		{"channel-flow",      required_argument, 0, 'C'},
		{"sedimentation",     required_argument, 0, 'S'},
		{"generate",          optional_argument, 0, 'g'},
		{"generate_basic",    optional_argument, 0, 'G'},
		{"random-seed",       required_argument, 0, 'a'},
		{"cluster-radius",    required_argument, 0, 'b'},
		{"volume-fraction",   required_argument, 0, 'p'},
		{"kn-kt-file",        required_argument, 0, 'k'},
		{"binary",            no_argument,       0, 'n'},
		{"name",              no_argument,       0, 'N'},
		{"identifier",        required_argument, 0, 'v'},
		{"force-to-run",      no_argument,       0, 'f'},
		{"diminish-output",   no_argument,       0, 'd'},
		{"checkpoint-file",   required_argument, 0, 'c'},
		{"dimer-file",   	  required_argument, 0, 'D'},
		{"help",              no_argument,       0, 'h'},
		{0, 0, 0, 0},
	};
	int index;
	int c;
	while ((c = getopt_long(argc, argv, "hn80fds:t:r:C:S:g::G:p:a:b:k:i:v:c:N:D:", longopts, &index)) != -1) {
		switch (c) {
			case 's':
				simulation_type = "shear rheology";
				control_variable = Parameters::ControlVariable::stress;
				control_value = Dimensional::str2DimensionalQty(Dimensional::Dimension::Stress, optarg, "shear stress");
				break;
			case 'r':
				simulation_type = "shear rheology";
				control_variable = Parameters::ControlVariable::rate;
				if (strcmp(optarg, "inf") == 0) {
					control_value = {Dimensional::Dimension::Force, 1, Dimensional::Unit::hydro};
					cout << "Rate control, infinite shear rate (hydro + hard contacts only)" << endl;
				} else if (strcmp(optarg, "zero") == 0) {
					control_value = {Dimensional::Dimension::Force, 0, Dimensional::Unit::kn};
					cout << "Rate control, zero shear rate (hydro + hard contacts only)" << endl;
				} else {
					control_value = Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, optarg, "shear rate");
				}
				break;
			case '8':
				cout << "'-r inf' for infinite shear rate (hydro + hard contacts only)" << endl;
				return 1;
			case '0':
				cout << "'-r zero' for zero shear rate (hydro + hard contacts only)" << endl;
				break;
			case 'C':
				simulation_type = "channel flow";
				control_variable = Parameters::ControlVariable::pressure_drop;
				control_value = Dimensional::str2DimensionalQty(Dimensional::Dimension::Stress, optarg, "pressure drop");
				break;
			case 'S':
				simulation_type = "sedimentation";
				control_variable = Parameters::ControlVariable::force;
				if (strcmp(optarg, "inf") == 0) {
					control_value = {Dimensional::Dimension::Force, 1, Dimensional::Unit::bodyforce};
					cerr << "Body force/repulsive force = infinity" << ' ' << control_value.value << endl;
				} else {
					cerr << "control variable = body force/repulsive force = " << ' ' << control_value.value << endl;
					control_value = Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, optarg, "force");
				}
				break;
			case 'k':
				input_files["knkt"] = optarg;
				break;
			case 'c':
				chkp_filename = optarg;
				break;
			case 'i':
				input_files["stress_rate"] = optarg;
				break;
			case 'g':
				generate_init = 1; // normal
				if (optarg) {
					if (optarg[0] == 'c') {
						generate_init = 2; // circular wide gap
					} else if (optarg[0] == 'w') {
						generate_init = 3; // simple shear with wall
					} else if (optarg[0] == 's') {
						generate_init = 4; // winding
					} else if (optarg[0] == 'b') {
						generate_init = 5; // bottom
					} else if (optarg[0] == 'f') {
						generate_init = 6; // filter mesh
					}
				}
				break;
			case 'G':
				generate_init = 7;
				str_vec = split(optarg, ':');
				break;
			case 'p':
				volume_frac_gen = atof(optarg);
				break;
			case 'a':
				random_seed = atoi(optarg);
				break;
			case 'b':
				/* To generate initial configuration. Radius of initial cluster.
				 * positive value: one cluster
				 * negative value: two clusters
				 */
				cluster_phi = atof(optarg);
				if (cluster_phi > 0) {
					cout << "To form one cluster\n";
				} else if (cluster_phi < 0) {
					cout << "To form two clusters\n";
				}
				break;
			case 'n':
				binary_conf = true;
				break;
			case 'v':
				simu_identifier = optarg;
				break;
			case 'f':
				force_to_run = true;
				break;
			case 'N':
				simu_name = optarg;
				break;
			case 'D':
				input_files["dimers"] = optarg;
				break;
			case 'h':
				cerr << usage << endl;
				exit(1);
			case '?':
				/* getopt already printed an error message. */
				break;
			default:
				abort();
		}
	}
	ostringstream in_args;
	for (int i=0; i<argc; i++) {
		in_args << argv[i] << " ";
	}
	if (generate_init >= 1) {
		GenerateInitConfig generate_init_config;
		if (generate_init != 7) {
			generate_init_config.generate(random_seed, volume_frac_gen, cluster_phi,
										  generate_init);
		} else {
			generate_init_config.generateBasic(random_seed, stof(str_vec[0]), stoi(str_vec[1]), stoi(str_vec[2]) == 2);
		}
	} else {
#ifdef SIGINT_CATCH
		std::signal(SIGINT, sigint_handler);
#endif
		if (optind == argc-2) {
			input_files["config"] = argv[optind++];
			input_files["params"] = argv[optind++];
		} else {
			cout << usage << endl;
			exit(1);
		}

		State::BasicCheckpoint state = State::zero_time_basicchkp;
		
		if (!chkp_filename.empty()) {
			state = State::readBasicCheckpoint(chkp_filename);
		}
		Simulation simulation(state);
		if (!chkp_filename.empty()) {
			simulation.simu_name = prepareSimulationNameFromChkp(chkp_filename);
		} else {
			simulation.simu_name = simu_name;
		}

		simulation.force_to_run = force_to_run;

		try {
			if (simulation_type == "shear rheology") {
				simulation.simulationSteadyShear(in_args.str(), input_files, binary_conf,
												 control_variable, control_value,
												 simu_identifier);
			} else {
				simulation.simulationFlowField(simulation_type,
											   in_args.str(), input_files, binary_conf,
											   control_variable, control_value,
											   simu_identifier);
			}
		} catch (runtime_error& e) {
			cout << e.what() << endl;
			return 1;
		}
	}
	cout << " Job done ok" << endl;
	return 0;
}
