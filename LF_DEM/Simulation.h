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
	int num_particle;
	int ts_max;
	double shear_strain;
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
	bool draw_rotation_2d;
	ofstream fout_yap;
	ofstream fout_vel;


	/*********************************************/
	void SetParameters(int argc, const char * argv[]);
	void output_yap();
	void output_vel();
	void initInteractingPair();
	void checkBreak();
	void checkContact();
	void timeEvolution();
	vec3d shiftUpCoordinate(double x, double y, double z);
public:
    /* For DEMsystem
     */
	Simulation();
	~Simulation();
	void SimulationMain(int argc, const char * argv[]);

};


#endif /* defined(__LF_DEM__Simulation__) */











