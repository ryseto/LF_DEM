//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 7/9/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__System__
#define __LF_DEM__System__

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <queue>
#include <list>
#include <string>
#include "vec3d.h"
#include "Bond.h"
#include "Cluster.h"
#include <fstream>
using namespace std;
class Cluster;

class System{
private:
	int np;
	double volume_fraction;
	int np3;
	double lx;
	double ly;
	double lz;
	double lx_half;
	double ly_half;
	double lz_half;
	double strain;
	double shear_disp;
	double dimensionless_shear_rate;
	vector <vector<int>> bond_list;
	vector <Cluster *> cl;
	ifstream fin_p;
	ifstream fin_i;
	ofstream fout_yap;
	ofstream fout_data;
protected:
public:
	~System(){};
	vector<vec3d> position;
	vector<vec3d> velocity;
	vector<vec3d> ang_velocity;
	vector<double> radius;
	vector<Bond> bond;
		
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
		bond_list.resize(np);
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
	
	void openYapFile(string filename){
		fout_yap.open(filename.c_str());
	}
	void openDataFile(string filename){
		fout_data.open(filename.c_str());
	}
	
	void importConfiguration(){
		char buf[1024];
		fin_p >> buf >> strain >> shear_disp >> dimensionless_shear_rate;
		int i;
		for (int j=0; j < np; j ++){
			fin_p >> i >> radius[i];
			fin_p >> position[i].x >> position[i].y >> position[i].z;
			fin_p >> velocity[i].x >> velocity[i].y >> velocity[i].z;
			fin_p >> ang_velocity[i].x >> ang_velocity[i].y >> ang_velocity[i].z;
			fin_p.getline(buf, 1024);
		}
		double strain_;
		int num_interaction;
		fin_i >> buf >> strain_ >> num_interaction;
		cerr << buf << ' ' << strain_ << ' ' << num_interaction << endl;
		cerr << strain << ' ' << strain_ << ' ' << "num_i" << num_interaction<< endl;
		bond.resize(num_interaction);
		double dv_sum = 0;
		for (int j=0; j < num_interaction; j ++){
			int i0, i1;
			bool contact;
			fin_i >> i0 >> i1 >> contact;
			fin_i.getline(buf, 1024);
			vec3d d_omega = ang_velocity[i0]-ang_velocity[i1];
			if (d_omega.norm() < 0.2){
				vec3d dv = velocity[i1] - velocity[i0];
				vec3d dp = position[i1] - position[i0];
				vec3d omega = 0.5*(ang_velocity[i0]+ang_velocity[i1]);
				if (dp.z > 10){
					dv.x -= lz;
					dp.z -= lz;
					dp.x -= shear_disp;
				} else if (position[i1].z - position[i0].z < -10){
					dv.x += lz;
					dp.z += lz;
					dp.x += shear_disp;
				}
				if (dp.x > 10){
					dp.x -= lx;
				} else if (dp.x < -10){
					dp.x += lx;
				}
				vec3d deviation = dv-cross(omega, dp);
				dv_sum += dv.norm();
				if (deviation.norm() < 0.1){
					bond_list[i0].push_back(j);
					bond_list[i1].push_back(j);
					bond[j].set(i0, i1, true);
				}
			}
		}
		cerr << "dv_ave = "<< dv_sum/num_interaction << endl;
		//		fout_interaction << "# " << sys.Shear_strain();
		//		fout_interaction << ' ' << cnt_interaction << endl;
	}
	
	void buildCluster(){
		vector <bool> checked;
		checked.resize(np);
		for(int i = 0; i < np; i++){
			checked[i] = false;
		}
		int i = 0;
		int numcluster = 0;
		cerr << cl.size();

		cl.clear();
		while (i < np){
			if (checked[i] == false){
				cl.push_back(new Cluster(i));
				checked[i] = true;
				int mem;
				int j_cnt = 0;
				while (true){
					mem = cl[numcluster]->get_member(j_cnt++);
					if (mem == -1){
						break;
					}
					for (int j=0; j < bond_list[mem].size(); j++){
						int b = bond_list[mem][j];
						int next = bond[b].next(mem);
						if (next >= 0){
							cl[numcluster]->add(next);
							checked[next] = true;
						}
					}
				}
				numcluster++;
			}
			i++;
		}
		cerr << endl;
	}
	
	void clear(){
		for (int j=0; j < np; j ++){
			bond_list[j].clear();
		}
	}
	
	vec3d calcAngVelocity(Cluster *cluster, double &std_dev){
		vec3d sum_ang_velocity(0,0,0);
		vec3d sum_sq_ang_velocity(0,0,0);
		
		int c_size = (int)cluster->size();
		for (int i=0; i<c_size; i++){
			int mem = cluster->get_member(i);
			sum_ang_velocity += ang_velocity[mem];
			sum_sq_ang_velocity.x += ang_velocity[mem].x*ang_velocity[mem].x;
			sum_sq_ang_velocity.y += ang_velocity[mem].y*ang_velocity[mem].y;
			sum_sq_ang_velocity.z += ang_velocity[mem].z*ang_velocity[mem].z;
		}
		vec3d ave_ang_vel =  (1./c_size)*sum_ang_velocity;
		vec3d ave_ang_vel_sq = (1./c_size)*sum_sq_ang_velocity;
		vec3d tmp(sqrt(ave_ang_vel.x*ave_ang_vel.x - ave_ang_vel_sq.x),
				  sqrt(ave_ang_vel.y*ave_ang_vel.y - ave_ang_vel_sq.y),
				  sqrt(ave_ang_vel.z*ave_ang_vel.z - ave_ang_vel_sq.z));

		std_dev = tmp.norm();
		return ave_ang_vel;
	}
	
	void output_yap(){
		vector <bool> done;
		done.resize(np);
		for (int i=0; i< np;i++){
			done[i] = false;
		}
		int num_cl_particle = 0;
		int sum_cl_size = 0;
		int cnt_cl = 0;
		int max_size = 0;
		for (int i=0; i<cl.size(); i++){
			int cl_size = (int)(cl[i]->size());

			double std_dev_ang;
			vec3d ave_ang_vel = calcAngVelocity(cl[i], std_dev_ang);
			if (cl_size > 1){
				sum_cl_size += cl_size;
				cnt_cl++;
				if (max_size < cl_size){
					max_size = cl_size ;
				}
				for (int j=0; j< cl_size; j++){
					int mem = cl[i]->get_member(j);
					if (abs(ave_ang_vel.y - ang_velocity[mem].y) > 0.3){
						fout_yap << "@ 0" << endl;
						fout_yap << "r 0.2" << endl;
						fout_yap << "c ";
						fout_yap << position[mem].x << ' ' ;
						fout_yap << position[mem].y -0.1<< ' ' ;
						fout_yap << position[mem].z << ' ' ;
						fout_yap << "@ 10" << endl;
					}
					fout_yap << "@ 10" << endl;
					fout_yap << "r " << radius[mem] << endl;
					fout_yap << "c ";
					fout_yap << position[mem].x << ' ' ;
					fout_yap << position[mem].y << ' ' ;
					fout_yap << position[mem].z << endl;
					done[mem] = true;
					num_cl_particle ++;
				}
			}
			
		}
		cout << endl;
		
		fout_yap << "@ 2" << endl;
		for (int i=0; i<np; i++){
			if (done[i] == false){
				fout_yap << "r " << radius[i] << endl;
				fout_yap << "c ";
				fout_yap << position[i].x << ' ' ;
				fout_yap << position[i].y << ' ' ;
				fout_yap << position[i].z << endl;
			}
		}
		fout_yap << endl;
		fout_data << strain << ' ' << (1.*num_cl_particle)/np << ' ' << sum_cl_size*(1./cnt_cl) << ' ';
		fout_data << cnt_cl << ' ' << max_size << endl;
	}
	
	
};

#endif /* defined(__LF_DEM__System__) */






