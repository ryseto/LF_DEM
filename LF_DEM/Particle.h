//
//  Particle.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Particle__
#define __LF_DEM__Particle__
#include <iostream>
#include "vec3d.h"

class Particle{
private:
	
protected:
	
public:
    /* For DEMsystem
     */
	Particle();
	~Particle(){
	};
	vec3d pos;
	vec3d vel;
	vec3d avel;

};
#endif /* defined(__LF_DEM__Particle__) */

