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
		cerr << "LF_DEM POSITIONFILE LUBCUTOFF [f/nf]" << endl;
		exit(1);
	}
	
	if (argv[1][0] == 'g'){
		generateInitialConfiguration(argc, argv);
	} else {
		Simulation simulation;
		simulation.SimulationMain(argc, argv);
	}
    return 0;
}
