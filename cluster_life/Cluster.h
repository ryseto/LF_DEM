//
//  Cluster.h
//  LF_DEM
//
//  Created by Ryohei Seto on 10/10/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Cluster__
#define __LF_DEM__Cluster__
#include <iostream>
#include <vector>
using namespace std;

class Cluster{
private:
	vector <int> member;
	
protected:
public:
	Cluster(){;}
	~Cluster(){;}
	void add(int i){
		member.push_back(i);
	}
	void printmember(){
		for (int i = 0; i<member.size(); i++){
			cerr << member[i] << ' ';
		}
		cerr << endl;
	}
	
};




#endif /* defined(__LF_DEM__Cluster__) */


