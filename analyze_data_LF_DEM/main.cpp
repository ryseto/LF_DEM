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
	ifstream fin_int;
	fin.open(argv[1]);
	fin.getline(buf, 1024);
	fin >> buf >> buf >> sys.np;
	fin >> buf >> buf >> sys.volume_fraction;
	fin >> buf >> buf >> sys.lx;
	fin >> buf >> buf >> sys.ly;
	fin >> buf >> buf >> sys.lz;
	fin_int.open(argv[2]);
	fin_int.getline(buf, 1024);
	fin_int.getline(buf, 1024);
	fin_int.getline(buf, 1024);
	fin_int.getline(buf, 1024);
	fin_int.getline(buf, 1024);
	
	sys.lx_half = sys.lx/2;
	sys.ly_half = sys.ly/2;
	sys.lz_half = sys.lz/2;
	int i;
	sys.pos.resize(sys.np);
	sys.radius.resize(sys.np);
	sys.allocate_g();
	sys.cnt_data = 0;
	sys.contri.resize(sys.np);
	while (fin >> buf >> sys.strain >> sys.shear_disp >> sys.dimensionless_shear_rate){
		cerr << sys.strain << endl;
		for (int j=0; j < sys.np; j ++){
			fin >> i >> sys.radius[i] >> sys.pos[i].x >> sys.pos[i].y >> sys.pos[i].z;
			fin.getline(buf, 1024);
			sys.contri[j] = 0;
		}
		int num_bond;
//		fin_int.getline(buf, 1024);
		fin_int >> buf >> buf >> num_bond;
//		cerr << buf << ' ' << num_bond << endl;
		int i0, i1;
		double f_cont_n;
		double f_cont_t;

		for (int j=0; j < num_bond; j ++){
			fin_int >> i0 >> i1 >> buf >> buf >> buf >> buf;
			fin_int >> buf >> buf >> f_cont_n >> f_cont_t;
			fin.getline(buf, 1024);
			double val = sqrt(f_cont_n*f_cont_n + f_cont_t+f_cont_t);
			sys.contri[i0] += val ;
			sys.contri[i1] += val ;
		}

//		fout_interaction << i << ' ' << j << ' '; // 1, 2
//		fout_interaction << sys.interaction[k].is_contact() << ' '; // 3
//		fout_interaction << nr_vec.x << ' '; // 4
//		fout_interaction << nr_vec.y << ' '; // 5
//		fout_interaction << nr_vec.z << ' '; // 6
//		fout_interaction << sys.interaction[k].Gap_nondim() << ' '; // 7
//		fout_interaction << sys.interaction[k].getLubForce() << ' '; // 8
//		fout_interaction << sys.interaction[k].getFcNormal() << ' '; // 9
//		fout_interaction << sys.interaction[k].getFcTan() << ' '; // 10
//		fout_interaction << sys.interaction[k].getColloidalForce() << ' '; // 11
//		fout_interaction << 6*M_PI*stress_contact.getStressXZ() << ' '; // 12
//		fout_interaction << 6*M_PI*stress_contact.getNormalStress1() << ' '; // 13
//		fout_interaction << 6*M_PI*stress_contact.getNormalStress2() << endl; // 14
//
		if (sys.strain > 3){
			sys.eval_pair_correlation();
		}
	}
	
	sys.calc_pair_correlation();
	
	
	
	
    return 0;
}

