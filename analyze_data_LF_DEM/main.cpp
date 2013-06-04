//
//  main.cpp
//  analyze_data_LF_DEM
//
//  Created by Ryohei Seto on 6/2/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#include <iostream>
#include <vector>
#include <fstream>
#include "SystemParticles.h"

using namespace std;
int main(int argc, const char * argv[])
{
	
	SystemParticles sys;
	char buf[1024];
	ifstream fin;
	fin.open(argv[1]);
	fin.getline(buf, 1024);
	fin >> buf >> buf >> sys.np;
	fin >> buf >> buf >> sys.volume_fraction;
	fin >> buf >> buf >> sys.lx;
	fin >> buf >> buf >> sys.ly;
	fin >> buf >> buf >> sys.lz;
	sys.lx_half = sys.lx/2;
	sys.ly_half = sys.ly/2;
	sys.lz_half = sys.lz/2;
	int i;
	sys.pos.resize(sys.np);
	sys.radius.resize(sys.np);
	sys.allocate_g();
	sys.cnt_snapshot = 0;
	
	while (fin >> buf >> sys.strain >> sys.shear_disp >> sys.dimensionless_shear_rate){
		cerr << sys.strain << endl;
		for (int j=0; j < sys.np; j ++){
			fin >> i >> sys.radius[i] >> sys.pos[i].x >> sys.pos[i].y >> sys.pos[i].z;
			fin.getline(buf, 1024);
		}
		if (sys.strain > 3){
			sys.cnt_snapshot ++;
			sys.eval_pair_correlation();
		}
	}
	
	sys.calc_pair_correlation();
	
	
	
	
    return 0;
}

