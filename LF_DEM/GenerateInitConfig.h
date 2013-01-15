//
//  GenerateInitConfig.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/6/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__GenerateInitConfig__
#define __LF_DEM__GenerateInitConfig__
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "vec3d.h"
using namespace std;

class GenerateInitConfig{
private:
	int random_seed;
	int dimension;
	double volume_fraction;
	double lx_lz, ly_lz;
	double lx, ly, lz;
	double lx2, ly2, lz2;

	double number_ratio;

	double a1;
	double a2;
	vector<vec3d> position;
	vector<double> radius;
	int np, np1, np2;
	vec3d dr;

	inline vec3d randUniformSphere(double r);
	double sqContactDistance(int i, int j, double contact_distance);
	bool is_contact(int i, int j);
	bool checkOverlap();
	void putRandom();
	void solveOverlap();
	void setParameters(int argc, const char * argv[]);
	void outputPositionData();
public:
	GenerateInitConfig(){
		dr.set(0,0,0);
	};
	int generate(int argc, const char * argv[]);

};

#endif /* defined(__LF_DEM__GenerateInitConfig__) */
