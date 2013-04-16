//
//  Stresslet.h
//  LF_DEM
//
//  Created by Ryohei Seto on 4/16/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef LF_DEM_Stresslet_h
#define LF_DEM_Stresslet_h
#include "vec3d.h"

class stresslet {
public:
	/* variables */
	double elm[6];
	
	inline stresslet(void){
		for (int i=0; i<6; i++){
			elm[i] = 0;
		}
	}
	
	inline stresslet(const vec3d & v1, const vec3d & v2){
		elm[0] = v1.x*v2.x;
		elm[1] = 0.5*(v1.x*v2.y+v1.y*v2.x);
		elm[2] = 0.5*(v1.x*v2.z+v1.z*v2.x);
		elm[3] = 0.5*(v1.y*v2.z+v1.z*v2.y);
		elm[4] = v1.y*v2.y;
		elm[5] = v1.z*v2.z;
	}
	
	inline stresslet(const double &_xx,
					 const double &_xy,
					 const double &_xz,
					 const double &_yz,
					 const double &_yy,
					 const double &_zz){
		elm[0] = _xx;
		elm[1] = _xy;
		elm[2] = _xz;
		elm[3] = _yz;
		elm[4] = _yy;
		elm[5] = _zz;
	}
	
	/* copy constructor */
	inline stresslet(const stresslet& other){
		for (int i=0; i<6; i++){
			elm[i] = other.elm[i];
		}
	}
	
	inline friend stresslet
	operator + (const stresslet &a1, const stresslet &a2){
		return stresslet(a1.elm[0]+a2.elm[0],
						 a1.elm[1]+a2.elm[1],
						 a1.elm[2]+a2.elm[2],
						 a1.elm[3]+a2.elm[3],
						 a1.elm[4]+a2.elm[4],
						 a1.elm[5]+a2.elm[5]);
	}
	
	inline friend stresslet
	operator + (const stresslet &s){
		return s;
	}
	
	/* subtraction */
	inline friend stresslet
	operator - (const stresslet &a1, const stresslet &a2){
		return stresslet(a1.elm[0]-a2.elm[0],
						 a1.elm[1]-a2.elm[1],
						 a1.elm[2]-a2.elm[2],
						 a1.elm[3]-a2.elm[3],
						 a1.elm[4]-a2.elm[4],
						 a1.elm[5]-a2.elm[5]);
	}
	
	inline friend stresslet
	operator - (const stresslet &s){
		return stresslet(s.elm[0],
						 s.elm[1],
						 s.elm[2],
						 s.elm[3],
						 s.elm[4],
						 s.elm[5]);
	}
	
	inline friend stresslet
	operator * (const double &d, const stresslet &s){
		return stresslet(d*s.elm[0],
						 d*s.elm[1],
						 d*s.elm[2],
						 d*s.elm[3],
						 d*s.elm[4],
						 d*s.elm[5]);
	}
	
	inline friend stresslet
	operator * (const stresslet &s, const double &d){
		return d*s;
	}
	
	inline friend stresslet
	operator * (const int &i, const stresslet &s){
		return stresslet(i*s.elm[0],
						 i*s.elm[1],
						 i*s.elm[2],
						 i*s.elm[3],
						 i*s.elm[4],
						 i*s.elm[5]);
	}
	
	inline friend stresslet
	operator * (const stresslet &s, const int &i){
		return stresslet(i*s.elm[0],
						 i*s.elm[1],
						 i*s.elm[2],
						 i*s.elm[3],
						 i*s.elm[4],
						 i*s.elm[5]);
	}
	
	inline stresslet&
	operator +=(const stresslet &s){
		for (int i=0; i<6; i++){
			elm[i] += s.elm[i];
		}
		return *this;
	}
	
	inline stresslet&
	operator -=(const stresslet &s){
		for (int i=0; i<6; i++){
			elm[i] -= s.elm[i];
		}
		return *this;
	}
	
	inline stresslet&
	operator *=(const double &d){
		for (int i=0; i<6; i++){
			elm[i] *= d;
		}
		return *this;
	}
	
	inline stresslet&
	operator *=(const int &i){
		for (int i=0; i<6; i++){
			elm[i] *= i;
		}
		return *this;
	}
	
	inline stresslet&
	operator /=(const double &d){
		double d_inv = 1/d;
		for (int i=0; i<6; i++){
			elm[i] *= d_inv;
		}
		return *this;
	}
	
	inline stresslet&
	operator /= (const int &i){
		double d_inv = 1./i;
		for (int i=0; i<6; i++){
			elm[i] *= d_inv;
		}
		return *this;
	}
	
	inline void
	set(const double &_xx, const double &_xy, const double &_xz,
		const double &_yz, const double &_yy, const double &_zz){
		elm[0] = _xx;
		elm[1] = _xy;
		elm[2] = _xz;
		elm[3] = _yz;
		elm[4] = _yy;
		elm[5] = _zz;
	}
	
	inline void
	reset(){
		for (int i=0; i<6; i++){
			elm[i] = 0;
		}
	}
	
	inline friend stresslet
	tensor_prod (const vec3d &v1, const vec3d &v2) {
		return stresslet(v1.x*v2.x,
						 0.5*(v1.x*v2.y+v1.y*v2.x),
						 0.5*(v1.x*v2.z+v1.z*v2.x),
						 0.5*(v1.y*v2.z+v1.z*v2.y),
						 v1.y*v2.y,
						 v1.z*v2.z);
	}
};
#endif
