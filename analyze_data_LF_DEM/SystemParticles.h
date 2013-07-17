//
//  SystemParticles.h
//  LF_DEM
//
//  Created by Ryohei Seto on 6/2/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__SystemParticles__
#define __LF_DEM__SystemParticles__
#include <iostream>
#include <vector>
#include "vec3d.h"

using namespace std;

class SystemParticles{
private:
	double d;
	int ix;
	int iz;
	double** g;
//	int cnt;
	double a1;
	double a2;
	int ix_max;
	int iz_max;

	double xrange;
	double zrange;
	double v_cell;
	double v1;
	double v2;
	double impact1;
	double impact2;
public:
	SystemParticles(){
//		d = 0.02;
		d = 0.05;
//		cnt = 0;
		xrange = 6;
		zrange = 6;
		a1 = 1;
		a2 = 1.4;
		v1 = M_PI*(4.0/3)*a1*a1*a1;
		v2 = M_PI*(4.0/3)*a2*a2*a2;
		v_cell = d*d*(2*d);
	}
	
	int np;
	int cnt_data;
	double volume_fraction;
	double lx;
	double ly;
	double lz;
	double lx_half;
	double ly_half;
	double lz_half;
	
	double strain;
	double shear_disp;
	double dimensionless_shear_rate;
	vector <double> radius;
	vector <vec3d> pos;
	vector <double> contri;
	void eval_pair_correlation();
	void allocate_g(){
		ix_max = 2*xrange/d+1;
		iz_max = zrange/d+1;
		g = new double* [ix_max];
		for (int ix=0; ix < ix_max; ix++){
			g[ix] = new double [iz_max];
		}
		for (int iz=0; iz < iz_max; iz++){
			for (int ix=0; ix < ix_max; ix++){
				g[ix][iz] = 0;
			}
		}

	}
	
	void calc_pair_correlation(){
		for (int iz=0; iz < iz_max; iz++){
			for (int ix=0; ix < ix_max; ix++){
				if (g[ix][iz] == 0){
					cout << 0 << ' ';
				} else {
					cout << log(g[ix][iz]/(cnt_data*volume_fraction)) << ' ';
				}
			}
			cout << endl;
		}
	}
	
};
#endif /* defined(__LF_DEM__SystemParticles__) */





