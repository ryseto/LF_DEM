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
#include <sstream>
#include "System.h"

class Simulation{
private:
	System sys;
	string filename_addition;
	double shear_strain_end;
	bool import_positions;
	vector< vec3d> initial_positions;
	vector <double> radii;
	string filename_import_positions;
	string filename_parameters;
	int np_a;
	int np_b;
	double radius_a;
	double radius_b;
	/*
	 * Resultant data
	 */
	double Viscosity;
	double N1;
	double N2;
	double Viscosity_h;
	double N1_h;
	double N2_h;
	double Viscosity_c;
	double N1_c;
	double N2_c;
	double Viscosity_b;
	double N1_b;
	double N2_b;
	double Viscosity_2;
	double Viscosity_2_h;
	double Viscosity_2_c;
	double N1_2;
	double N2_2;
	
	/*
	 *  Simulation parameters
	 */
	double strain_interval_out;
	/*
	 * For output data.
	 */
	ofstream fout_vpy;
	ofstream fout_rheo;
	ofstream fout_particle;
	ofstream fout_interaction;
	bool out_data_particle;
	bool out_data_interaction;
	bool origin_zero_flow;
	/*
	 *
	 */
	void timeEvolution();
	void evaluateData();
	/*
	 * For inputs
	 */
	void setDefaultParameters();
	void readParameterFile();
	void openOutputFiles();
	void prepareSimulationName();
	void autoSetParameters(const string &keyword,
						   const string &value);
	void importInitialPositionFile();
	/*
	 * For outputs
	 */
	void output_yap();
	void output_vel();
	void outputDataHeader(ofstream &fout);
	void initContactPair();
	void outputRheologyData();
	void outputConfigurationData();
	
	vec3d shiftUpCoordinate(double x, double y, double z);
	void drawLine2(char type , const vec3d &pos1, const vec3d &pos2, ofstream &fout);
	void drawLine(char type , const vec3d &pos, const vec3d &vec, ofstream &fout);
	void drawLine(double x0, double y0, double z0,
				  double x1, double y1, double z1,
				  ofstream &fout);
public:
    /* For DEMsystem
     */
	Simulation();
	~Simulation();
	void SimulationMain(int argc, const char * argv[]);
};
#endif /* defined(__LF_DEM__Simulation__) */
