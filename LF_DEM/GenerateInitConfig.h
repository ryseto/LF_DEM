//
//  GenerateInitConfig.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/6/12.
//  Copyright (c) 2012-2014 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class GenerateInitConfig
 \brief Class to generate initial configurations (option -g). Not used by shear simulations.
 \author Ryohei Seto
 \author Romain Mari
 */


#ifndef __LF_DEM__GenerateInitConfig__
#define __LF_DEM__GenerateInitConfig__
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "vec3d.h"
#include "System.h"
#include "ParameterSet.h"
#include "Events.h"

class GenerateInitConfig{
private:
	System sys;
	ParameterSet p;
	std::list<Event> events;

	char disperse_type;
	double volume_fraction;
	double volume_fraction1;
	double volume_fraction2;
	double vf_ratio; // = volume_fraction1/volume_fraction;
	double lx_lz;
	double ly_lz;
	double lx;
	double ly;
	double lz;
	double lx_half;
	double ly_half;
	double lz_half;
    double z_top;
    double z_bot;
	double a1;
	double a2;
	vec3d* grad;
	vec3d* prev_grad;
	double gradientDescent();
	double computeGradient();
	void moveAlongGradient(vec3d*, int);
	void storeGradient();
	double step_size;
	int rand_seed;
	/*
	 * This can generate monodisperse or binary configurations consisting of two types of 
	 * particles. 
	 * It is indicated by "bidisperse" option.
	 * Ones have magnetic susceptibility 1
	 * The others have magnetic susceptibility -1.
	 */
	bool magnetic_config;
	bool circulargap_config;
    bool parallel_wall_config;
	double radius_in;
	double radius_out;
	int np_wall1;
	int np_wall2;
	int np_fix;
	int np_movable;
#ifndef USE_DSFMT
	MTRand rand_gen;
#endif
#ifdef USE_DSFMT
	dsfmt_t rand_gen;
#endif
	double zeroTMonteCarloSweep();
	int overlapNumber(int);
	double particleEnergy(int);
	void updateInteractions(int);
	std::vector<vec3d> position;
	std::vector<double> radius;
	int np;
	int np1;
	int np2;
	vec3d dr;
	inline vec3d randUniformSphere(double r);
	inline vec3d randUniformCircle(double r);
	double sqContactDistance(int i, int j, double contact_distance);
	void putRandom();
	void setParameters();
	void outputPositionData();

public:
	GenerateInitConfig():
	sys(System(p, events)),
	magnetic_config(false),
	circulargap_config(false),
    parallel_wall_config(false) {};
	
	int generate(int rand_seed_, int config_type);
	/* config_type = 1 -- noraml
	 * config_type = 2 -- circular wide gap
     * config_type = 3 -- simple shear with wall
     * config_type = 10 -- magnetic
	 */
};
#endif /* defined(__LF_DEM__GenerateInitConfig__) */
