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
//#include "System.h"
using namespace std;

class Cluster{
private:
	vector <int> member;
protected:
public:
	Cluster(){;}
	Cluster(int i);
	
	~Cluster(){;}
	
	void add(int i){
		for (int j=0; j<member.size(); j++){
			if (member[j]== i){
				return;
			}
		}
		member.push_back(i);
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
