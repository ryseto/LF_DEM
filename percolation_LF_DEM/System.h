//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 7/24/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__System__
#define __LF_DEM__System__

#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include <list>
#include <string>
#include "vec3d.h"

class System{
private:
	double volume_fraction;
	double lx;
	double ly;
	double lz;
	double lx_half;
	double ly_half;
	double lz_half;

	double shear_disp;
	double dimensionless_shear_rate;
	ifstream fin_p;
	ifstream fin_i;
	ofstream fout_yap;
	ofstream fout_data;
	int **contact_table;
	vector < vector<int> > contact_list;
	
protected:
public:
	int np;
	int max_cluster_size;
	double strain;
	~System(){};
	vector<vec3d> position;
	vector<vec3d> velocity;
	vector<vec3d> ang_velocity;
	vector<double> radius;
	
	void setImportFile_par(string filename){
		fin_p.open(filename.c_str());
		char buf[128];
		fin_p.getline(buf, 1024);
		fin_p >> buf >> buf >> np;
		fin_p >> buf >> buf >> volume_fraction;
		fin_p >> buf >> buf >> lx;
		fin_p >> buf >> buf >> ly;
		fin_p >> buf >> buf >> lz;
		lx_half = lx/2;
		ly_half = ly/2;
		lz_half = lz/2;
		cerr << "np = " << np << endl;
		cerr << lx << ' ' << ly << ' ' << lz << endl;
		position.resize(np);
		radius.resize(np);
		velocity.resize(np);
		ang_velocity.resize(np);
		
		
		contact_table = new int* [np];
		for (int i=0; i<np; i++){
			contact_table[i] = new int [np];
		}
		contact_list.resize(np);
	}
	
	void setImportFile_int(string filename){
		fin_i.open(filename.c_str());
		char buf[128];
		fin_i.getline(buf, 1024);
		fin_i.getline(buf, 1024);
		fin_i.getline(buf, 1024);
		fin_i.getline(buf, 1024);
		fin_i.getline(buf, 1024);
		fin_i.getline(buf, 1024);
		cerr << "======" << endl;
	}
	
	void openDataFile(string filename){
		fout_data.open(filename.c_str());
	}
	
	void importConfiguration(){
		char buf[1024];
		fin_p >> buf >> strain >> shear_disp >> dimensionless_shear_rate;
		for (int j=0; j < np; j ++){
			int i;
			fin_p >> i >> radius[i];
			fin_p >> position[i].x >> position[i].y >> position[i].z;
//			fin_p >> velocity[i].x >> velocity[i].y >> velocity[i].z;
//			fin_p >> ang_velocity[i].x >> ang_velocity[i].y >> ang_velocity[i].z;
			fin_p.getline(buf, 1024);
		}
		double strain_;
		int num_interaction;
		fin_i >> buf >> strain_ >> num_interaction;
		cerr << buf << ' ' << strain_ << ' ' << num_interaction << endl;
		cerr << strain << ' ' << strain_ << ' ' << "num_i" << num_interaction<< endl;
		for (int i=0; i<np; i++){
			contact_list[i].clear();
		}

		for (int j=0; j<num_interaction; j++){
			int i0, i1;
			bool contact = false;
			fin_i >> i0 >> i1 >> contact;
			fin_i.getline(buf, 1024);

			if (contact){
				contact_table[i0][i1] = 1;
				contact_table[i1][i0] = 1;
				contact_list[i0].push_back(i1);
				contact_list[i1].push_back(i0);
			}
		}
		//		fout_interaction << "# " << sys.Shear_strain();
		//		fout_interaction << ' ' << cnt_interaction << endl;
	}
	
	
	
	void check_percolation(int i_start){
		vector <int> contact_cluster;
		queue <int> to_be_checked;
		contact_cluster.push_back(i_start);
		to_be_checked.push(i_start);
		while (!to_be_checked.empty()){
			int i_check = to_be_checked.front();
			to_be_checked.pop();
			for (int j=0; j<contact_list[i_check].size(); j++){
				int i_new = contact_list[i_check][j];
				bool new_member = true;
				for (int jj=0; jj<contact_cluster.size(); jj++){
					if (contact_cluster[jj] == i_new){
						new_member = false;
						break;
					}
				}
				if (new_member == true){
					contact_cluster.push_back(i_new);
					to_be_checked.push(i_new);
				}
			}
			contact_list[i_check].clear();
		}
		if (contact_cluster.size() > 1){
			if (contact_cluster.size() > max_cluster_size){
				max_cluster_size = (int)contact_cluster.size();
			}
		}
	}

};



#endif /* defined(__LF_DEM__System__) */

