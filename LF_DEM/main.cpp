//
//  main.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include <iostream>
#include "Simulation.h"
using namespace std;
int main(int argc, const char * argv[])
{
    
	Simulation simulation;
	simulation.SimulationMain(argc, argv);
    return 0;
}
