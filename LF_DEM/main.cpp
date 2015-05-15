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
#include "Simulation.h"
#include "GenerateInitConfig.h"
#ifndef GIT_VERSION
#include "VersionInfo.h"
#endif
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//void incompatibility_exiting(string a, string b)
//{
//	cerr << a << " and " << b << " not compatible " << endl;
//	exit(1);
//}


int main(int argc, char **argv)
{
	cerr << endl << "LF_DEM version " << GIT_VERSION << endl << endl;
	string usage = "(1) Simulation\n $ LF_DEM [-r Rate ] [-s Stress ] \
	[-R Rate_Sequence ] [-S Stress_Sequence ] [-k kn_kt_File] [-i Provisional_Data] [-n] \
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
		{"repulsion", required_argument, 0, 'r'},
		{"rep-seq-file", required_argument, 0, 'R'},
		{"stress-controlled", required_argument, 0, 's'},
		{"stress-seq-file", required_argument, 0, 'S'},
		{"generate", no_argument, 0, 'g'},
		{"kn-kt-file", required_argument, 0, 'k'},
		{"binary", no_argument, 0, 'n'},
		{"magnetic", no_argument, 0, 'm'},
		{"help", no_argument, 0, 'h'},
		{0,0,0,0},
	};

	int index;
	int c;
	while ((c = getopt_long(argc, argv, "hng:s:S:r:R:k:i:h:", longopts, &index)) != -1) {
		switch (c) {
			case 's':
				rheology_control = "stress";
				dimensionless_number = atof(optarg);
				cerr << "Stress control: " << dimensionless_number << endl;
				break;
			case 'S':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				rheology_control = "stress";
				seq_filename = optarg;
				seq_type = "s";
				cerr << "Stress sequence, file " << seq_filename << endl;
				break;
			case 'r':
				rheology_control = "rate";
				getSuffix(optarg, numeral, suffix);
				dimensionless_number = stof(numeral);
				cerr << "Rate control: " << dimensionless_number << endl;
				break;
			case 'R':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				rheology_control = "rate";
				seq_filename = optarg;
				seq_type = "r";
				cerr << "Rate sequence, file " << seq_filename << endl;
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
		generate_init_config.generate(random_seed);
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

		if (seq_filename == "not_given") {
			simulation.simulationSteadyShear(in_args.str(), input_files, binary_conf,
											 dimensionless_number, suffix, rheology_control);
		} 
		else {
		  cerr << " User def sequence temporarily disabled " << endl;
		  //		  simulation.simulationUserDefinedSequence(seq_type, in_args.str(), input_files, binary_conf, rheology_control);
		}
	}
	return 0;
}
