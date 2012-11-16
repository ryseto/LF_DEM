//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Interaction__
#define __LF_DEM__Interaction__
#include <iostream>

#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Interaction.h"
#include "System.h"
using namespace std;
class System;

class Interaction{
private:
	System *sys;
	vec3d r_vec;
	
protected:
	
public:
    /* For DEMsystem
     */
	Interaction(){};
	~Interaction(){};
	double r;
	bool active;
	int particle_num[2];
	vec3d nr_vec;
	
	void init(System *sys_);
	void create(int i, int j);
	void calcInteraction();

};


#endif /* defined(__LF_DEM__Interaction__) */
