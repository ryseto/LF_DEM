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

void mainConventional(int argc, char **argv);
void mainLammpsLike(int argc, char **argv);

std::string prepareSimulationNameFromChkp(const std::string& filename_chkp)
{
	/**
	 \brief Determine simulation name when starting from checkpoint
	 */
	string::size_type chkp_name_start = filename_chkp.rfind("chk_") + 4;
	string::size_type chkp_name_end = filename_chkp.rfind(".dat");

	return filename_chkp.substr(chkp_name_start, chkp_name_end-chkp_name_start);
}

int main(int argc, char **argv)
{
	if (argc <= 3) {
		mainLammpsLike(argc, argv);
	} else {
		mainConventional(argc, argv);
	}
	return 0;
}

void mainLammpsLike(int argc, char **argv)
{
	bool binary_conf = false;
	bool force_to_run = false;
	string chkp_filename = "";
	const struct option longopts[] = {
		{"binary",          no_argument, 0, 'n'},
		{"force-to-run",    no_argument, 0, 'f'},
		{"checkpoint-file", no_argument, 0, 'c'}
	};
	int index;
	int c;
	while ((c = getopt_long(argc, argv, "nfc", longopts, &index)) != -1) {
		switch (c) {
			case 'n':
				binary_conf = true;
				break;
			case 'c':
				chkp_filename = optarg;
				break;
			case 'f':
				force_to_run = true;
				break;
			default:
				abort();
		}
	}
	State::BasicCheckpoint state = State::zero_time_basicchkp;
	if (!chkp_filename.empty()) {
		state = State::readBasicCheckpoint(chkp_filename);
	}
	Simulation simulation(state);
	simulation.force_to_run = force_to_run;
	simulation.simulationMain(argv[optind], binary_conf);
	return;
}

void mainConventional(int argc, char **argv)
{
	cout << endl << "LF_DEM version " << GIT_VERSION << endl << endl;
	string usage = "(1) Simulation\n $ LF_DEM [-r Rate] [-s Stress] [-R Rate_Sequence] [-S Stress_Sequence]\
	[-e] [-m ?] [-k kn_kt_File] [-v Simulation_Identifier] [-i Provisional_Data] [-n]\
	Configuration_File Parameter_File \
	\n\n OR \n\n(2) Generate initial configuration\n $ LF_DEM [-a Random_Seed] [-p Volume_Fraction] -g[c/w/s]\n";

	int generate_init = 0;
	string type_init_config = "normal";

	int random_seed = 1;
	double volume_frac_gen = 0;
	bool binary_conf = false;
	bool force_to_run = false;
	string config_filename = "not_given";
	string param_filename = "not_given";
	string knkt_filename = "not_given";
	string stress_rate_filename = "not_given";
	string chkp_filename = "";
	string simu_name;
	
	Dimensional::DimensionalQty<double> control_value;
	Parameters::ControlVariable control_variable = Parameters::ControlVariable::rate;
	string simu_identifier = "";
	const struct option longopts[] = {
		{"rate-controlled",   required_argument, 0, 'r'},
		{"rate-infty",        required_argument, 0, '8'},
		{"stress-controlled", required_argument, 0, 's'},
		{"generate",          optional_argument, 0, 'g'},
		{"random-seed",       required_argument, 0, 'a'},
		{"volume-fraction",   required_argument, 0, 'p'},
		{"kn-kt-file",        required_argument, 0, 'k'},
		{"binary",            no_argument,       0, 'n'},
		{"name",              no_argument,       0, 'N'},
		{"identifier",        required_argument, 0, 'v'},
		{"force-to-run",      no_argument,       0, 'f'},
		{"diminish-output",   no_argument,       0, 'd'},
		{"checkpoint-file",   no_argument,       0, 'c'},
		{"help",              no_argument,       0, 'h'},
		{0, 0, 0, 0},
	};
	
	int index;
	int c;
	while ((c = getopt_long(argc, argv, "hn80fds:t:r:g::p:a:k:i:v:c:N:", longopts, &index)) != -1) {
		switch (c) {
			case 's':
				control_variable = Parameters::ControlVariable::stress;
				control_value = Dimensional::str2DimensionalQty(Dimensional::Dimension::Stress, optarg, "shear stress");
				break;
			case 'r':
				control_variable = Parameters::ControlVariable::rate;
				control_value = Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, optarg, "shear rate");
				break;
			case '8':
				control_variable = Parameters::ControlVariable::rate;
				control_value = {Dimensional::Dimension::Force, 1, Dimensional::Unit::hydro};
				cout << "Rate control, infinite shear rate (hydro + hard contacts only)" << endl;
				break;
			case '0':
				control_variable = Parameters::ControlVariable::rate;
				control_value = {Dimensional::Dimension::Force, 0, Dimensional::Unit::kn};
				cout << "Rate control, zero shear rate (hydro + hard contacts only)" << endl;
				break;
			case 'k':
				knkt_filename = optarg;
				break;
			case 'c':
				chkp_filename = optarg;
				break;
			case 'i':
				stress_rate_filename = optarg;
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
					}
				}
				break;
			case 'p':
				volume_frac_gen = atof(optarg);
				break;
			case 'a':
				random_seed = atoi(optarg);
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
		generate_init_config.generate(random_seed, volume_frac_gen, generate_init);
	} else {
#ifdef SIGINT_CATCH
		std::signal(SIGINT, sigint_handler);
#endif
		if (optind == argc-2) {
			config_filename = argv[optind++];
			param_filename = argv[optind++];
		} else {
			cerr << usage << endl;
			exit(1);
		}
		vector <string> input_files(5, "not_given");
		input_files[0] = config_filename;
		input_files[1] = param_filename;
		input_files[2] = knkt_filename;
		input_files[3] = stress_rate_filename;
		
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
			simulation.simulationSteadyShear(in_args.str(), input_files, binary_conf,
											 control_variable, control_value,
											 simu_identifier);
		} catch (runtime_error& e) {
			cerr << e.what() << endl;
			return 1;
		}
	}
	cerr << " Job done ok" << endl;
	return;
}
