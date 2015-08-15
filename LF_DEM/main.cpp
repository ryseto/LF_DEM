//
//  main.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <iostream>
#include <string>
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

int main(int argc, char **argv)
{
	cout << endl << "LF_DEM version " << GIT_VERSION << endl << endl;
	string usage = "(1) Simulation\n $ LF_DEM [-r Rate ] [-s Stress ] \
	[-R Rate_Sequence ] [-S Stress_Sequence ] [-m ?] [-k kn_kt_File] [-i Provisional_Data] [-n] \
	Configuration_File Parameter_File \n\n OR \n\n(2) Generate initial configuration\n $ LF_DEM -g Random_Seed\n";
	
	double dimensionless_number = 0;
	string numeral, suffix;
		
	bool generate_init = false;
	int random_seed = 1;

	bool binary_conf = false;
	
	string config_filename = "not_given";
	string param_filename = "not_given";
	string knkt_filename = "not_given";
	string stress_rate_filename = "not_given";
	string seq_filename = "not_given";
	string seq_type;
	string rheology_control = "rate";
	const struct option longopts[] = {
		{"rate-controlled", required_argument, 0, 'r'},
		{"rate-infty", required_argument, 0, '8'},
		{"rate-seq-file", required_argument, 0, 'R'},
		{"stress-controlled", required_argument, 0, 's'},
		{"stress-seq-file", required_argument, 0, 'S'},
		{"magnetic", required_argument, 0, 'm'},
		{"magnetic", no_argument, 0, 'M'},
		{"generate", required_argument, 0, 'g'},
		{"kn-kt-file", required_argument, 0, 'k'},
		{"binary", no_argument, 0, 'n'},
		{"help", no_argument, 0, 'h'},
		{0,0,0,0},
	};

	int index;
	int c;
	while ((c = getopt_long(argc, argv, "hnM8m:g:s:S:t:r:R:k:i:h:", longopts, &index)) != -1) {
		switch (c) {
			case 's':
				rheology_control = "stress";
				if (getSuffix(optarg, numeral, suffix)) {
					dimensionless_number = atof(numeral.c_str());
					cout << "Stress control: " << dimensionless_number << endl;
				} else {
					errorNoSuffix("shear stress");
				}
				break;
			case 'S':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				rheology_control = "stress";
				seq_filename = optarg;
				seq_type = "s";
				cout << "Stress sequence, file " << seq_filename << endl;
				break;
			case 't':
				rheology_control = "stress";
				if (getSuffix(optarg, numeral, suffix)) {
					dimensionless_number = atof(numeral.c_str());
					cout << "Stress control: " << dimensionless_number << endl;
				} else {
					errorNoSuffix("shear stress");
				}
				seq_type = "iy";
				break;
			case 'r':
				rheology_control = "rate";
				if (getSuffix(optarg, numeral, suffix)) {
					dimensionless_number = atof(numeral.c_str());
					cout << "Rate control: " << dimensionless_number << endl;
				} else {
					errorNoSuffix("shear rate");
				}
				break;
			case 'R':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				rheology_control = "rate";
				seq_filename = optarg;
				seq_type = "r";
				cout << "Rate sequence, file " << seq_filename << endl;
				break;
			case '8':
				rheology_control = "rate";
				dimensionless_number = 1;
				suffix = "h";
				cout << "Rate control, infinite shear rate (hydro + hard contacts only)" << endl;
				break;
			case 'm':
				/*
				 * magnetic moment m
				 * Typical magnetic force: (3 mu m^2)/(4 pi (2a)^4)
				 * Typical Brownian force: kT/a
				 * Dimensionless_number (Pe number) can be defined as the ratio between these two forces:
				 * Typical magnetic force/Typical Brownian force 
				 * dimensionless_number = Pe_M = (3 mu m^2) / (64 pi kT a^3)
				 */
				rheology_control = "magnetic";
				if (getSuffix(optarg, numeral, suffix)) {
					dimensionless_number = atof(numeral.c_str());
					cout << "Magnetic field control: " << dimensionless_number << endl;
				} else {
					errorNoSuffix("magnetic field");
				}
				cerr << "Magnetic simulation" << endl;
				break;
			case 'k':
				knkt_filename = optarg;
				break;
			case 'i':
				stress_rate_filename = optarg;
				break;
 			case 'g':
				generate_init = true;
				random_seed = atoi(optarg);
				break;
			case 'M':
				rheology_control = "magnetic";
				cerr << "Magnetic simulation" << endl;
				break;
 			case 'n':
				binary_conf = true;
				break;
			case 'h':
				cerr << usage << endl;
				exit(1);
			case '?':
				/* getopt already printed an error message. */
				break;
			default:
				abort ();
		}
	}
	ostringstream in_args;
	for (int i=0; i<argc; i++) {
		in_args << argv[i] << " ";
	}
	if (generate_init) {
		GenerateInitConfig generate_init_config;
		bool magnetic_config = false;
		if (rheology_control == "magnetic") {
			magnetic_config = true;
		}
		generate_init_config.generate(random_seed, magnetic_config);
	} else {
		if (optind == argc-2) {
			config_filename = argv[optind++];
			param_filename = argv[optind++];
		} else {
			cerr << usage << endl;
			exit(1);
		}
		vector <string> input_files;
		input_files.resize(5);
		input_files[0] = config_filename;
		input_files[1] = param_filename;
		input_files[2] = knkt_filename;
		input_files[3] = stress_rate_filename;
		input_files[4] = seq_filename;
		Simulation simulation;
		if (rheology_control == "magnetic") {
			simulation.simulationMagnetic(in_args.str(), input_files, binary_conf,
										  dimensionless_number, suffix, rheology_control);
		} else if (seq_type == "iy") {
			simulation.simulationInverseYield(in_args.str(), input_files, binary_conf,
											  dimensionless_number, suffix, rheology_control);
			
		} else if (seq_filename == "not_given") {
			simulation.simulationSteadyShear(in_args.str(), input_files, binary_conf,
											 dimensionless_number, suffix, rheology_control);
		} else {
		  cerr << " User def sequence temporarily disabled " << endl;
		  //		  simulation.simulationUserDefinedSequence(seq_type, in_args.str(), input_files, binary_conf, rheology_control);
		}
	}
	return 0;
}
