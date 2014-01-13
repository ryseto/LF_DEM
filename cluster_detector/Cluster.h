//
//  Cluster.h
//  LF_DEM
//
//  Created by Ryohei Seto on 7/9/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Cluster__
#define __LF_DEM__Cluster__

#include <iostream>
#include <vector>
#include "vec3d.h"
//#include "System.h"
using namespace std;

class Cluster{
private:
	

protected:
public:
	Cluster(){;}
	Cluster(int i, vec3d pos_);
	
	~Cluster(){;}
	vector <int> member;
	vector <vec3d> pos;
	vector <int> cbond;
	void add(int i, int next, vec3d dr_, int b){
		vec3d pos_orig;
		for (int j=0; j<member.size(); j++){
			if (member[j]== i){
				return;
			}
			if (member[j] == next){
				pos_orig = pos[j];
			}
		}
		cbond.push_back(b);
		member.push_back(i);
		pos.push_back(pos_orig + dr_);
		//pos.push_back(dr_);
	}

	unsigned long size(){
		return member.size();
	}
	
	int get_member(int j){
		if (j < member.size()){
			return member[j];
		} else {
			return -1;
		}
	}


	
	
};

#endif /* defined(__LF_DEM__Cluster__) */
