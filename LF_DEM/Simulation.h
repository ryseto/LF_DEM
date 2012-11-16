//
//  Simulation.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Simulation__
#define __LF_DEM__Simulation__

#include <iostream>
#include <fstream>
#include <queue>
#include "Interaction.h"
#include "System.h"

class Simulation{
private:
	System sys;
	int dimension;
	int num_particle;
	int ts_max;
	double eta;
	double dt;
	/*
	 * Interparticle interactions
	 */
	Interaction *interaction;
	int **interacting_pair; // Table
	int max_num_interaction; // Allowed length of interaction array.
	int num_interaction; // Length of used interaction array.
	queue<int> deactivated_interaction;
	/*
	 *  Simulation parameters
	 */
	double cutoff_distance;
	/*
	 * For output data.
	 */
	ofstream fout_tmp;

	/*********************************************/
	void output_tmp();
	void initInteractingPair();
	void checkBreak();
	void checkContact();
	void timeEvolution();
	void checkPeriodicBoundary();

public:
    /* For DEMsystem
     */
	Simulation();
	~Simulation();
	void run();

};


#endif /* defined(__LF_DEM__Simulation__) */











