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

int main(int argc, char **argv)
{
	int c;
	string usage = "(1) Simulation\n $ LF_DEM [-p Peclet ] [-c Scaled_Critical_Load ] [-r Scaled_Repulsion ] [-k kn_kt_File] Configuration_File Parameter_File \n\n OR (2) Generate initial configuration\n $ LF_DEM -g\n\n Note: $ LF_DEM -c -1 corresponds to infinite shear rate";
	
	double Peclet = 0;
	double scaled_repulsion = 0;
	double scaled_critical_load = 0;
	
	bool generate_init = false;
	bool knkt = false;
	
	string config_filename;
	string param_filename;
	string knkt_filename;
	
	while ((c = getopt(argc, argv, "ghp:r:c:k:")) != -1) {
		switch (c) {
			case 'p':
				Peclet = atof(optarg);
				cerr << "Peclet " << Peclet << endl;
				break;
			case 'r':
				scaled_repulsion = atof(optarg);
				cerr << "scaled repulsion " << scaled_repulsion << endl;
				break;
			case 'c':
				scaled_critical_load = atof(optarg);
				cerr << "scaled critical load " << scaled_critical_load << endl;
				break;
			case 'k':
				knkt = true;
				knkt_filename = optarg;
				break;
			case 'g':
				generate_init = true;
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
		string *input_files;
		int fnb;
		if (knkt) {
			fnb = 3;
			input_files = new string [fnb];
			input_files[0] = config_filename;
			input_files[1] = param_filename;
			input_files[2] = knkt_filename;
		} else {
			fnb = 2;
			input_files = new string [fnb];
			input_files[0] = config_filename;
			input_files[1] = param_filename;
		}
		if (scaled_repulsion>0 && scaled_critical_load > 0) {
			cerr << " Repulsion AND Critical Load cannot be used at the same time" << endl;
			exit(1);
		}
		Simulation simulation;
		simulation.simulationConstantShearRate(fnb, input_files, Peclet,
											   scaled_repulsion, scaled_critical_load);
	}
}
