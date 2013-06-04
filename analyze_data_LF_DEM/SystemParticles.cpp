//
//  SystemParticles.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 6/2/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#include "SystemParticles.h"

void
SystemParticles::eval_pair_correlation(){
	for (int i=0; i < np-1; i++){
		for (int j=i+1; j < np; j++){
			if (true || radius[i] == 1 && radius[j] == 1){
				vec3d dr = pos[i]-pos[j];
				if (abs(dr.y) > ly_half) {
					if (dr.y > 0){
						dr.y -= ly;
					} else {
						dr.y += ly;
					}
				}
				
				if ( abs(dr.y) < d ){
					if (abs(dr.z) > lz_half) {
						if (dr.z > 0) {
							dr.z -= lz;
							dr.x -= shear_disp;
						} else {
							dr.z += lz;
							dr.x += shear_disp;
						}
					}
					if (abs(dr.x) > lx_half) {
						if (dr.x > 0){
							dr.x -= lx;
						} else {
							dr.x += lx;
						}
					}
					if (dr.z < 0){
						dr.z = -dr.z;
						dr.x = -dr.x;
					}
					if (dr.z < zrange &&
						abs(dr.x) < xrange){
						int ix = (dr.x+xrange)/d;
						int iz = dr.z/d;

						g[ix][iz] += 1;
						
						cnt ++;
						
					}
				}
			}
		}
		
	}
	
	
}
