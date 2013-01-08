//
//  Simulation.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__Simulation__
#define __LF_DEM__Simulation__
#include <iostream>
#include <fstream>
#include <queue>
#include "System.h"

class Simulation{
private:
	System sys;

	double shear_strain;
	bool import_positions;
	vector< vec3d> initial_positions;
	string filename_import_positions;
	string filename_parameters;
	/*
	 *  Simulation parameters
	 */
	int ts_max;
	int interval_snapshot;
	/*
	 * For output data.
	 */
	ofstream fout_yap;
	ofstream fout_vpy;
	ofstream fout_rheo;
	double yap_force_factor;
	bool origin_zero_flow;
	void SetDefaultParameters(int argc, const char * argv[]);
	void ReadParameterFile(int argc, const char * argv[]);
	void SetParametersPostProcess();


	void AutoSetParameters(const string &keyword,
						   const string &value);
	void importInitialPositionFile();
	void output_yap();
	void output_vpython(double);
	void output_vel();
	void initContactPair();
	void outputRheologyData();
	void timeEvolution();
	vec3d shiftUpCoordinate(double x, double y, double z);
	void drawLine2(char type , vec3d pos1, vec3d pos2, ofstream &fout);
	void drawLine(char type , vec3d pos, vec3d vec, ofstream &fout);
	void drawLine(double x0, double y0, double z0,
				  double x1, double y1, double z1,
				  ofstream &fout);
	/*
	 * Genrate initial configuration
	 * LF_DEM g phi lx ly lz
	 */
public:
    /* For DEMsystem
     */
	Simulation();
	~Simulation();
	void SimulationMain(int argc, const char * argv[]);
};
#endif /* defined(__LF_DEM__Simulation__) */
