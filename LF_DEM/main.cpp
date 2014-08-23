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


int main(int argc, char **argv)
{
	string usage = "(1) Simulation\n $ LF_DEM [-p Peclet_Num ] [-c Scaled_Critical_Load ] [-r Scaled_Repulsion ] [-a Scaled_Cohesion] [-k kn_kt_File] Configuration_File Parameter_File \n\n OR \n\n (2) Generate initial configuration\n $ LF_DEM -g\n";
	
	double peclet_num = 0;
	double scaled_repulsion = 0;
	double scaled_cohesion = 0;
	double scaled_critical_load = 0;
	
	bool generate_init = false;
	bool knkt = false;
	bool sequence = false;
	
	string config_filename;
	string param_filename;
	string knkt_filename;
	bool strain_controlled = true;
	bool stress_controlled = !strain_controlled;
	
	string seq_filename;
	string seq_type;
	
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
		{"stress-controlled", no_argument, 0, 's'},
		{"help", no_argument, 0, 'h'},
		{0,0,0,0},
	};
	int index;
	int c;
	while ((c = getopt_long(argc, argv, "ghsp:P:r:R:c:C:a:A:k:", longopts, &index)) != -1) {
		switch (c) {
			case 'p':
				peclet_num = atof(optarg);
				cerr << "Peclet number " << peclet_num << endl;
				break;
			case 'P':
				if (sequence) { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				sequence = true;
				seq_filename = optarg;
				seq_type = "p";
				cerr << "Peclet sequence, file " << seq_filename << endl;
				break;
			case 'r':
				scaled_repulsion = atof(optarg);
				cerr << "scaled repulsion " << scaled_repulsion << endl;
				break;
			case 'R':
				if (sequence) { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				sequence = true;
				seq_filename = optarg;
				seq_type = "r";
				cerr << "scaled repulsion sequence, file " << seq_filename << endl;
				break;
			case 'a':
				scaled_cohesion = atof(optarg);
				cerr << "scaled cohesion " << scaled_cohesion << endl;
				break;
			case 'A':
				if (sequence) { cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				sequence = true;
				seq_filename = optarg;
				seq_type = "a";
				cerr << "scaled cohesion sequence, file " << seq_filename << endl;
				break;
			case 'c':
				scaled_critical_load = atof(optarg);
				cerr << "scaled critical load " << scaled_critical_load << endl;
				break;
			case 'C':
				if(sequence){ cerr << " Only one parameter sequence allowed " << endl; exit(1);};
				sequence = true;
				seq_filename = optarg;
				seq_type = "c";
				cerr << "scaled critical load sequence, file " << seq_filename << endl;
				break;
			case 'k':
				knkt = true;
				knkt_filename = optarg;
				break;
			case 'g':
				generate_init = true;
				break;
			case 's':
				stress_controlled = true;
				strain_controlled = !stress_controlled;
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
	
	if (generate_init) {
		GenerateInitConfig generate_init_config;
		generate_init_config.generate();
	} else {
		if (optind == argc-2) {
			config_filename = argv[optind++];
			param_filename = argv[optind++];
		} else {
			cerr << usage << endl;
			exit(1);
		}
		vector <string> input_files;
		input_files.push_back(config_filename);
		input_files.push_back(param_filename);
		
		if (knkt) {
			input_files.push_back(knkt_filename);
		}
		if (sequence) {
			input_files.push_back(seq_filename);
		}
		
		if (scaled_repulsion > 0 && scaled_critical_load > 0) {
			cerr << " Repulsion AND Critical Load cannot be used at the same time" << endl;
			exit(1);
		}
		
		Simulation simulation;
		if (!sequence) {
			if (strain_controlled) {
				simulation.simulationSteadyShear(input_files, peclet_num,
												 scaled_repulsion, scaled_cohesion, scaled_critical_load, "strain");
			} else if (stress_controlled) {
				simulation.simulationSteadyShear(input_files, peclet_num,
												 scaled_repulsion, scaled_cohesion, scaled_critical_load, "stress");
			}
		} else {
			if (strain_controlled) {
				simulation.simulationUserDefinedSequence(seq_type, input_files, "strain");
			} else if (stress_controlled) {
				simulation.simulationUserDefinedSequence(seq_type, input_files, "stress");
			}
		}
		
	}
}
