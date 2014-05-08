//
//  main.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <iostream>
#include "Simulation.h"
#include "GenerateInitConfig.h"

int main(int argc, const char * argv[])
{
	if (argc == 1){
		cerr << "usage:" << endl;
		cerr << "(1) Simulation" << endl;
		cerr << "  $ LF_DEM ARG1 ARG2 ARG3" << endl;
		cerr << " ARG1: file name of initial configuration" << endl;
		cerr << " ARG2: parameter file" << endl;
		cerr << " ARG3: dimensionless shear rate" << endl;
		cerr << "       (`-1' is considered as infinity)" << endl;
		cerr << "(2) Generate initial configuration" << endl;
		cerr << "  $ LF_DEM g" << endl;
		return 0;
	}
	if (argv[1][0] == 'g'){
		GenerateInitConfig generate_init_config;
		generate_init_config.generate();
	} else {
		Simulation simulation;
		simulation.simulationConstantShearRate(argc, argv);
	}
    return 0;
}
