//
//  MagneticParticle.h
//  LF_DEM
//
//  Created by Zhiyuan Zhao on 2020/9/14.
//  Copyright © 2020 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__MagneticParticle__
#define __LF_DEM__MagneticParticle__

#include <vector>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "RateDependence.h"

class MagneticParticle {
public:
    std::vector <vec3d> dipole_orient;
    vec3d dipole_orient_p0;
    vec3d dipole_orient_p1;
    
    MagneticParticle(){};
    MagneticParticle(unsigned np) :
    dipole_orient(np)
    {};
        
    void calcDipoleOrientation(int i, double angle) {
        dipole_orient[i].set(cos(angle),0,sin(angle));
    }
    
    void setPairDipole(int p0, int p1) {
        dipole_orient_p0 = dipole_orient[p0];
        dipole_orient_p1 = dipole_orient[p1];
    }
};

#endif /* __LF_DEM__MagneticParticle__ */
