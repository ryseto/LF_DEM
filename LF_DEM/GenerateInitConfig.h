//
//  GenerateInitConfig.h
//  LF_DEM
//
//  Created by Ryohei Seto on 12/6/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__GenerateInitConfig__
#define __LF_DEM__GenerateInitConfig__
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "vec3d.h"
using namespace std;
int random_seed;
int dimension;
double volume_fraction;
double lx, ly, lz;
double lx2, ly2, lz2;
vector<vec3d> position;
int np;
vec3d dr;

inline vec3d randUniformSphere(double r){
	double z = 2*drand48() - 1.0;
	double phi = 2*M_PI*drand48();
	double sin_theta = sqrt(1.0-z*z);
	return vec3d( r*sin_theta*cos(phi),  r*sin_theta*sin(phi), r*z);
}

double sqContactDistance(int i, int j){
	dr.x = position[i].x - position[j].x;
	dr.z = position[i].z - position[j].z;
	if (dr.z > 3 ){
		dr.z -= lz;
	} else if (dr.z < -3){
		dr.z += lz;
	}
	if (abs(dr.z) < 2){
		while(dr.x > 3){
			dr.x -= lx;
		}
		while(dr.x < - lx2){
			dr.x += lx;
		}
		if (abs(dr.x) < 2){
			if (dimension == 3){
				dr.y = position[i].y - position[j].y;
				if (dr.y > 3 ){
					dr.y -= ly;
				} else if (dr.y < -3){
					dr.y += ly;
				}
				if (abs(dr.y) < 2){
					return dr.sq_norm();
				}
			} else {
				return dr.sq_norm();
			}
		}
	}
	return 100;
}

bool checkOverlap(){
	static int i_previous = 0;
	static int j_previous = 1;
	if ( sqContactDistance(i_previous, j_previous) < 4){
		return true;
	}
	for (int i = 0; i < np ; i++){
		for (int j = i+1; j < np ; j++){
			if ( sqContactDistance(i, j) < 4){
				i_previous = i ;
				j_previous = j ;
				return true;
			}
		}
	}
	return false;
}

void putRandom(){
	srand48(random_seed);
	for (int i=0; i < np; i++){
		position[i].x = lx*drand48();
		position[i].z = lz*drand48();
		if (dimension == 2){
			position[i].y = ly2;
		} else {
			position[i].y = ly*drand48();
		}
	}
}

void solveOverlap(){
	int cc = 0;
	vector<int> previous_overlap;
	previous_overlap.resize(np);
	double dd = 0.01;
	vec3d delta_translation;
	while (true){
		int i = lrand48() % np;
		if (dimension == 2){
			double rand_angle = 2*M_PI*drand48();
			delta_translation.set(dd*cos(rand_angle), 0, dd*sin(rand_angle));
		} else {
			delta_translation = randUniformSphere(dd);
		}
		position[i] += delta_translation;
		position[i].periodicBoundaryBox(lx, ly, lz);

		
		
		int overlap = -1;
		if (sqContactDistance(i, previous_overlap[i]) < 4){
			overlap =  previous_overlap[i];
		} else {
			for (int j = 0; j < np ; j++){
				if (j != i){
					if ( sqContactDistance(i, j) < 4){
						previous_overlap[i] = j;
						overlap = j;
						break;
					}
				}
			}
		}
		if (overlap == -1){
			if (cc > 10000){
				if (cc % 1000 == 0 ){
					if ( checkOverlap() == false){
						break;
					}
				}
			}
			cc ++;
		} else {
			position[i] += dd*dr;
			position[i].periodicBoundaryBox(lx, ly, lz);
			position[overlap] -= dd*dr;
			position[overlap].periodicBoundaryBox(lx, ly, lz);
		}
	}
}

void outputPositionData(){
	ofstream fout;
	ostringstream ss_posdatafilename;
	ss_posdatafilename << "D" << dimension;
	if (dimension == 2){
		ss_posdatafilename << "L" << lx << "_" << lz;
	} else {
		ss_posdatafilename << "L" << lx << "_" << ly << "_" << lz;
	}
	ss_posdatafilename << "vf" << volume_fraction << ".dat";
	cerr << ss_posdatafilename.str() << endl;
	fout.open(ss_posdatafilename.str().c_str());
	for (int i = 0; i < np ; i++){
		fout << position[i].x << ' ';
		fout << position[i].y << ' ';
		fout << position[i].z << endl;
	}
	fout.close();
}

void setParameters(int argc, const char * argv[]){
	if (argc == 6){
		dimension = 2;
		volume_fraction = atof(argv[2]);
		lx = atof(argv[3]);
		ly = 0;
		lz = atof(argv[4]);
		np = (int)(volume_fraction*lx*lz / M_PI);
		random_seed = atoi(argv[5]);
	} else if (argc == 7){
		dimension = 3;
		volume_fraction = atof(argv[2]);
		lx = atof(argv[3]);
		ly = atof(argv[4]);
		lz = atof(argv[5]);
		double box_volume  = lx*ly*lz;
		double sphere_volume = (4.0/3)*M_PI;
		np = (int)(volume_fraction*box_volume /sphere_volume);
		random_seed = atoi(argv[6]);
	} else {
		cerr << "2D: LF_DEM g vf lx lz randomseed" << endl;
		cerr << "3D: LF_DEM g vf lx ly lz randomseed" << endl;
		exit(1);
	}
	lx2 = lx/2;
	ly2 = ly/2;
	lz2 = lz/2;
	cerr << "np = " << np << endl;
	cerr << "vf = " << volume_fraction << endl;
	cerr << "box =" << lx << ' ' << ly << ' ' << lz << endl;
}

int generateInitialConfiguration(int argc, const char * argv[]){
	setParameters(argc, argv);
	position.resize(np);
	putRandom();
	solveOverlap();
	outputPositionData();
	return 0;
}


#endif /* defined(__LF_DEM__GenerateInitConfig__) */
