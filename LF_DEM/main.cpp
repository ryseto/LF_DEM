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
#define no_argument 0
#define required_argument 1
#define optional_argument 2

void incompatibility_exiting(string a, string b){
	cerr << a << " and " << b << " not compatible " << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	string usage = "(1) Simulation\n $ LF_DEM [-p Peclet_Num ] [-c Scaled_Critical_Load ] \
	[-r Scaled_Repulsion ] [-s Stress ] [-a Scaled_Cohesion ] \
	[-S Stress_Sequence ] [-k kn_kt_File] [-i Provisional_Data] \
	Configuration_File Parameter_File \n\n OR \n\n(2) Generate initial configuration\n $ LF_DEM -g Random_Seed\n";
	
	bool peclet = false;
	double peclet_num = 0;

	bool repulsion = false;
	double ratio_repulsion = 0;

	bool critical_load = false;
	double ratio_critical_load = 0;

	bool cohesion = false;
	double ratio_cohesion = 0;
	
	bool generate_init = false;
	int random_seed = 1;
	string config_filename = "not_given";
	string param_filename = "not_given";
	string knkt_filename = "not_given";
	string stress_rate_filename = "not_given";
	string seq_filename = "not_given";
	string seq_type;
	string rheology_control = "rate";
	const struct option longopts[] = {
		{"peclet", required_argument, 0, 'p'},
		{"pe-seq-file", required_argument, 0, 'P'},
		{"critical-load", required_argument, 0, 'c'},
		{"cload-seq-file", required_argument, 0, 'C'},
		{"repulsion", required_argument, 0, 'r'},
		{"rep-seq-file", required_argument, 0, 'R'},
		{"cohesion", required_argument, 0, 'a'},
		{"coh-seq-file", required_argument, 0, 'A'},
		{"generate", no_argument, 0, 'g'},
		{"kn-kt-file", required_argument, 0, 'k'},
		{"stress-controlled", required_argument, 0, 's'},
		{"stress-seq-file", required_argument, 0, 'S'},
		{"help", no_argument, 0, 'h'},
		{0,0,0,0},
	};
	int index;
	int c;
	while ((c = getopt_long(argc, argv, "hg:s:S:p:P:r:R:c:C:a:A:b:B:k:i:", longopts, &index)) != -1) {
		switch (c) {
			case 'p':
				peclet = true;
				peclet_num = atof(optarg);
				cerr << "Brownian, Peclet number " << peclet_num << endl;
				break;
			case 'P':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				peclet = true;
				seq_filename = optarg;
				seq_type = "p";
				cerr << "Brownian, Peclet sequence, file " << seq_filename << endl;
				break;
			case 's':
				repulsion = true;
				rheology_control = "stress";
				ratio_repulsion = atof(optarg);
				cerr << "Repulsion, stress " << ratio_repulsion << endl;
				break;
			case 'S':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				repulsion = true;
				rheology_control = "stress";
				seq_filename = optarg;
				seq_type = "s";
				cerr << "Stress sequence, file " << seq_filename << endl;
				break;
			case 'r':
				//if (!stress_controlled) {
				rheology_control = "rate";
				ratio_repulsion = atof(optarg);
				repulsion = true;
				cerr << "Repulsion, shear rate " << ratio_repulsion << endl;
				//} else {
				//	cerr << "option -r ignored for stress controlled simulations " << endl;
				//}
				break;
			case 'R':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				rheology_control = "rate";
				repulsion = true;
				seq_filename = optarg;
				seq_type = "r";
				cerr << "Repulsion, sequence, file " << seq_filename << endl;
				break;
			case 'a':
				cohesion = true;
				rheology_control = "rate";
				ratio_cohesion = atof(optarg);
				cerr << "Cohesion, shear rate " << ratio_cohesion << endl;
				break;
			case 'A':
				if (seq_filename != "not_given") { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				cohesion = true;
				rheology_control = "rate";
				seq_filename = optarg;
				seq_type = "a";
				cerr << "Cohesion, sequence, file " << seq_filename << endl;
				break;
			case 'b':
				cohesion = true;
				rheology_control = "stress";
				ratio_cohesion = atof(optarg);
				cerr << "Cohesion, shear rate " << ratio_cohesion << endl;
				break;

			case 'c':
				critical_load = true;
				ratio_critical_load = atof(optarg);
				cerr << "Critical load, shear rate " << ratio_critical_load << endl;
				break;
			case 'C':
				if(seq_filename != "not_given"){ cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				critical_load = true;
				seq_filename = optarg;
				seq_type = "c";
				cerr << "Critical load, sequence, file " << seq_filename << endl;
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
	// Incompatibilities
	if (peclet && rheology_control == "stress") {
		incompatibility_exiting("peclet", "stress_controlled");
//	} else if (cohesion && rheology_control == "stress") {
//		incompatibility_exiting("cohesion", "stress_controlled");
	} else if (critical_load && rheology_control == "stress") {
		incompatibility_exiting("critical_load", "stress_controlled");
	} else if (critical_load && repulsion) {
		incompatibility_exiting("critical_load", "repulsion");
	} else if (peclet && cohesion) {
		incompatibility_exiting("peclet", "cohesion");
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
			simulation.simulationSteadyShear(input_files, peclet_num,
											 ratio_repulsion, ratio_cohesion, ratio_critical_load,
											 rheology_control);
		} else {
			simulation.simulationUserDefinedSequence(seq_type, input_files, rheology_control);
		}
	}
	return 0;
}
