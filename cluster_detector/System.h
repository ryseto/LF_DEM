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

	ofstream fout_yap;
	ofstream fout_data;
	ofstream fout_cluster;
protected:
public:
	~System(){};
	vector<vec3d> position;
	vector<double> orientation;
	vector<vec3d> velocity;
	vector<vec3d> ang_velocity;
	vector<double> radius;
	vector<Bond> bond;
	ifstream fin_p;
	ifstream fin_i;
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
		orientation.resize(np);
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
	}
	
	void openYapFile(string filename){
		fout_yap.open(filename.c_str());
	}
	void openDataFile(string filename){
		fout_data.open(filename.c_str());
	}
	void openClusterFile(string filename){
		fout_cluster.open(filename.c_str());
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
			fin_p >> buf >> buf >> buf;
			fin_p >> orientation[i];
			fin_p.getline(buf, 1024);
		//	cerr << ang_velocity[i].y << endl;
		}
		double strain_;
		int num_interaction;
		fin_i >> buf >> strain_ >> num_interaction;
		cerr << buf << ' ' << strain_ << ' ' << num_interaction << endl;
		cerr << strain << ' ' << strain_ << ' ' << "num_i" << num_interaction<< endl;
		bond.resize(num_interaction);
		double dv_sum = 0;
		static int cnt_strain = 0;
		for (int j=0; j < num_interaction; j ++){
			int i0, i1;
			bool contact;
			double nx, ny, nz, gap;
			double f_lub_n, f_lub_t;
			double f_con_n, f_con_t, f_ele;
			double s_xz, ns1, ns2;
			int contact_state;
			/* 1, 2: numbers of the interacting particles
			 * 3: 1=contact, 0=apart
			 * 4, 5, 6: normal vector
			 * 7: dimensionless gap = s - 2, s = 2r/(a1+a2)
			 * 8: normal     of lubrication force
			 * 9: tangential of lubrication force
			 * 10: normal part     of contact force
			 * 11: tangential part of contact force
			 * 12: normal colloidal force
			 * 13: Viscosity contribution of contact xF
			 * 14: N1 contribution of contact xF
			 * 15: N2 contribution of contact xF
			 * 16: friction state
			 *      0 = not frictional
			 *      1 = non-sliding
			 *      2 = sliding
			 */
			fin_i >> i0 >> i1 >> contact; // 1 2 3
			fin_i >> nx >> ny >> nz >> gap; //4 5 6 7
			fin_i >> f_lub_n >> f_lub_t; // 8 9
			fin_i >> f_con_n >> f_con_t; // 10 11
			fin_i >> f_ele; //12
			fin_i >> s_xz >> ns1 >> ns2; // 13 14 15
			fin_i >> contact_state; //16
			//cerr << contact_state << endl;
			if (contact_state == 1) {
				//vec3d dw = ang_velocity[i1]-ang_velocity[i0];
				double dw = ang_velocity[i1].y-ang_velocity[i0].y;
				cout << cnt_strain << ' '  << dw << endl;
				if (abs(dw) < 2){
					bond_list[i0].push_back(j);
					bond_list[i1].push_back(j);
					bond[j].set(i0, i1, contact_state);
				}
			}
		}
		cnt_strain ++;
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
		cerr << cl.size() << endl;
		cl.clear();
		
		while (i < np){
			if (checked[i] == false){
				cl.push_back(new Cluster(i, position[i]));
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
						if (bond[b].contact_state == 1){
							int next = bond[b].next(mem);
							if (next >= 0){
																
								vec3d dr = position[next]-position[mem];
								if (abs(dr.z) > 10 ){
									if (dr.z > 0){
										dr.z -= lz;
										dr.x -= shear_disp;
									} else {
										dr.z += lz;
										dr.x += shear_disp;
									}
								}
								if (abs(dr.x) > 10 ){
									if (dr.x > 0){
										dr.x -= lx;
									} else {
										dr.x += lx;
									}
								}
								if (abs(dr.y) > 10 ){
									if (dr.y > 0){
										dr.y -= ly;
									} else {
										dr.y += ly;
									}
								}
								
								cl[numcluster]->add(next, mem, dr, b);
								checked[next] = true;
							}
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
		int num_cl_particle = 0;
//		int sum_cl_size = 0;
//		int cnt_cl = 0;
//		int max_size = 0;
//		vec3d comp_axis(-1,0,1);
//		comp_axis.unitvector();
//		
//		vec3d elong_axis(1,1,0);
//		elong_axis.unitvector();
//		vector <double> cluster_length;
//		for (int i=0; i<cl.size(); i++){
//			int cl_size = (int)(cl[i]->size());
//			if (cl_size > 1){
//				vec3d cm(0,0,0);
//				for (int j=0; j< cl_size; j++){
//					cm += cl[i]->pos[j];
//				}
//				cm = cm /cl_size;
//				for (int j=0; j< cl[i]->size(); j++){
//					cl[i]->pos[j]-=cm;
//				}
//				double comp_max = 0;
//				double comp_min = 0;
//				for (int j=0; j< cl[i]->size(); j++){
//					double comp_length = dot(comp_axis, cl[i]->pos[j]);
//					if (comp_length > comp_max){
//						comp_max = comp_length;
//					}
//					if (comp_length < comp_min){
//						comp_min = comp_length;
//					}
//				}
//				cluster_length.push_back(comp_max-comp_min);
//			}
//		}
//
//		double ave_cluster_length = 0;
//		double max_cluster_length = 0;
//		if (cluster_length.size() > 0){
//			for (int k=0; k < cluster_length.size(); k++){
//				ave_cluster_length += cluster_length[k];
//				if (max_cluster_length < cluster_length[k]){
//					max_cluster_length = cluster_length[k];
//				}
//			}
//			ave_cluster_length = ave_cluster_length/cluster_length.size();
//		}
		
//		fout_yap << "@ 3" << endl;
//		
//		for (int i=0; i<cl.size(); i++){
//			int cl_size = (int)(cl[i]->size());
//			if (cl_size > 1){
//				sum_cl_size += cl_size;
//				cnt_cl++;
//				if (max_size < cl_size){
//					max_size = cl_size ;
//				}
//				vec3d cm(0,0,0);
//				for (int j=0; j< cl_size; j++){
//					cm += cl[i]->pos[j];
//				}
//				cm = cm /cl_size;
//				fout_yap << "@ 10" << endl;
//				for (int j=0; j< cl[i]->size(); j++){
//					fout_yap << "r " << radius[cl[i]->get_member(j)] << endl;
//					fout_yap << "c ";
//					fout_yap << cl[i]->pos[j].x << ' ' ;
//					fout_yap << cl[i]->pos[j].y << ' ' ;
//					fout_yap << cl[i]->pos[j].z << endl;
//					num_cl_particle ++;
//				}
//			}
//		}
//		fout_yap << endl;
//


		fout_yap << "@ 3" << endl;
		fout_yap << "y 3" << endl;
		for (int i=0; i<cl.size(); i++){
			int cl_size = (int)(cl[i]->size());
			if (cl_size > 1){
				for (int j=0; j< cl[i]->size()-1; j++){
					int b = cl[i]->cbond[j];
					int i0 = bond[b].par_num[0];
					int i1 = bond[b].par_num[1];
					vec3d dp = position[i0] - position[i1];
					if (dp.sq_norm() < 10){
						double dw = abs(ang_velocity[i0].y- ang_velocity[i1].y);
						fout_yap << "r " << 0.2*dw << endl;
						fout_yap << "s ";
						fout_yap << position[i0].x << ' ' ;
						fout_yap << position[i0].y - 0.1<< ' ' ;
						fout_yap << position[i0].z << ' ' ;
						fout_yap << position[i1].x << ' ' ;
						fout_yap << position[i1].y - 0.1<< ' ' ;
						fout_yap << position[i1].z << endl;
					}
					//	fout_yap << cl[i]->pos[j].z << endl;
					//num_cl_particle ++;
				}
			}
		}
		fout_yap << "@ 2" << endl;
		for (int i=0; i<cl.size(); i++){
			int cl_size = (int)(cl[i]->size());
			if (cl_size > 1){
				//				sum_cl_size += cl_size;
				//				cnt_cl++;
				//				if (max_size < cl_size){
				//					max_size = cl_size ;
				//				}
				for (int j=0; j< cl[i]->size(); j++){
					fout_yap << "r " << radius[cl[i]->get_member(j)] << endl;
					fout_yap << "c ";
					fout_yap << cl[i]->pos[j].x << ' ' ;
					fout_yap << cl[i]->pos[j].y << ' ' ;
					fout_yap << cl[i]->pos[j].z << endl;
					num_cl_particle ++;
				}
			}
		}
		fout_yap << "y 4" << endl;
		fout_yap << "@ 0" << endl;
		for (int i=0; i<cl.size(); i++){
			int cl_size = (int)(cl[i]->size());
			if (cl_size > 1){
				//				sum_cl_size += cl_size;
				//				cnt_cl++;
				//				if (max_size < cl_size){
				//					max_size = cl_size ;
				//				}
				for (int j=0; j< cl[i]->size(); j++){
					int ii = cl[i]->member[j];
					fout_yap << "l ";
					fout_yap << cl[i]->pos[j].x << ' ' ;
					fout_yap << cl[i]->pos[j].y -0.1 << ' ' ;
					fout_yap << cl[i]->pos[j].z << ' ' ;
					fout_yap << cl[i]->pos[j].x + sin(orientation[ii])<< ' ' ;
					fout_yap << cl[i]->pos[j].y -0.1<< ' ' ;
					fout_yap << cl[i]->pos[j].z + cos(orientation[ii])<< endl;
					num_cl_particle ++;
				}
			}
		}
		
		
		
		fout_yap << endl;
		
//		fout_data << strain << ' ';
//		fout_data << ave_cluster_length << ' ';
//		fout_data << cluster_length.size() << ' ';
//		fout_data << max_cluster_length << endl;

	}
	
	void output_clusters(){
		if (cl.size() > 1){
			int count_cluster = 0;
			for (int i=0; i<cl.size(); i++){
				int cl_size = (int)(cl[i]->size());
				if (cl_size > 1){
					count_cluster ++;
					//for (int j=0; j< cl_size; j++){
					//fout_cluster << cl[i]->get_member(j) << ' ';
					//
					//}
					//fout_cluster << "\n" ;
					fout_cluster << cl_size << ' ';
				}
			}
			if (count_cluster > 0){
				fout_cluster << "\n";
			}
		}
	}
		
	
};

#endif /* defined(__LF_DEM__System__) */






