//
//  StressTensor.h
//  LF_DEM
//
//  Created by Ryohei Seto on 4/16/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef LF_DEM_StressTensor_h
#define LF_DEM_StressTensor_h
#include "vec3d.h"

class StressTensor {
private:

public:
	/*
	 * (xx, xy, xz, yz, yy, zz)
	 */
	double elm[6];
	
	inline StressTensor(void){
		for (int i=0; i<6; i++){
			elm[i] = 0;
		}
	}
	
	void outputStressTensor(ofstream &fout){
		fout << elm[0] << ' ';
		fout << elm[1] << ' ';
		fout << elm[2] << ' ';
		fout << elm[3] << ' ';
		fout << elm[4] << ' ';
		fout << elm[5] << ' ';
	}
	
	inline StressTensor(const vec3d & v1, const vec3d & v2){
		elm[0] = v1.x*v2.x;
		elm[1] = 0.5*(v1.x*v2.y+v1.y*v2.x);
		elm[2] = 0.5*(v1.x*v2.z+v1.z*v2.x);
		elm[3] = 0.5*(v1.y*v2.z+v1.z*v2.y);
		elm[4] = v1.y*v2.y;
		elm[5] = v1.z*v2.z;
	}
	
	inline StressTensor(const double &_xx,
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
	inline StressTensor(const StressTensor& other){
		for (int i=0; i<6; i++){
			elm[i] = other.elm[i];
		}
	}
	
	inline friend StressTensor
	operator + (const StressTensor &a1, const StressTensor &a2){
		return StressTensor(a1.elm[0]+a2.elm[0],
							a1.elm[1]+a2.elm[1],
							a1.elm[2]+a2.elm[2],
							a1.elm[3]+a2.elm[3],
							a1.elm[4]+a2.elm[4],
							a1.elm[5]+a2.elm[5]);
	}
	
	inline friend StressTensor
	operator + (const StressTensor &s){
		return s;
	}
	
	/* subtraction */
	inline friend StressTensor
	operator - (const StressTensor &a1, const StressTensor &a2){
		return StressTensor(a1.elm[0]-a2.elm[0],
							a1.elm[1]-a2.elm[1],
							a1.elm[2]-a2.elm[2],
							a1.elm[3]-a2.elm[3],
							a1.elm[4]-a2.elm[4],
							a1.elm[5]-a2.elm[5]);
	}
	
	inline friend StressTensor
	operator - (const StressTensor &s){
		return StressTensor(s.elm[0],
							s.elm[1],
							s.elm[2],
							s.elm[3],
							s.elm[4],
							s.elm[5]);
	}
	
	inline friend StressTensor
	operator * (const double &d, const StressTensor &s){
		return StressTensor(d*s.elm[0],
							d*s.elm[1],
							d*s.elm[2],
							d*s.elm[3],
							d*s.elm[4],
							d*s.elm[5]);
	}
	
	inline friend StressTensor
	operator * (const StressTensor &s, const double &d){
		return d*s;
	}
	
	inline friend StressTensor
	operator * (const int &i, const StressTensor &s){
		return StressTensor(i*s.elm[0],
							i*s.elm[1],
							i*s.elm[2],
							i*s.elm[3],
							i*s.elm[4],
							i*s.elm[5]);
	}
	
	inline friend StressTensor
	operator * (const StressTensor &s, const int &i){
		return StressTensor(i*s.elm[0],
							i*s.elm[1],
							i*s.elm[2],
							i*s.elm[3],
							i*s.elm[4],
							i*s.elm[5]);
	}
	
	inline StressTensor&
	operator +=(const StressTensor &s){
		for (int i=0; i<6; i++){
			elm[i] += s.elm[i];
		}
		return *this;
	}
	
	inline StressTensor&
	operator -=(const StressTensor &s){
		for (int i=0; i<6; i++){
			elm[i] -= s.elm[i];
		}
		return *this;
	}
	
	inline StressTensor&
	operator *=(const double &d){
		for (int i=0; i<6; i++){
			elm[i] *= d;
		}
		return *this;
	}
	
	inline StressTensor&
	operator *=(const int &i){
		for (int j=0; j<6; j++){
			elm[j] *= i;
		}
		return *this;
	}
	
	inline StressTensor&
	operator /=(const double &d){
		double d_inv = 1/d;
		for (int i=0; i<6; i++){
			elm[i] *= d_inv;
		}
		return *this;
	}
	
	inline StressTensor&
	operator /= (const int &i){
		double d_inv = 1./i;
		for (int j=0; j<6; j++){
			elm[j] *= d_inv;
		}
		return *this;
	}
	
	/*
	 * Finite trace
	 */
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
	
	/*
	 * Traceless
	 */
	inline void
	set(const double &_xx, const double &_xy, const double &_xz,
		const double &_yz, const double &_yy){
		elm[0] = _xx;
		elm[1] = _xy;
		elm[2] = _xz;
		elm[3] = _yz;
		elm[4] = _yy;
		elm[5] = -(_xx+_yy);
	}
	
	inline void
	set(const vec3d &v1, const vec3d &v2){
		elm[0] = v1.x*v2.x;
		elm[1] = 0.5*(v1.x*v2.y+v1.y*v2.x);
		elm[2] = 0.5*(v1.x*v2.z+v1.z*v2.x);
		elm[3] = 0.5*(v1.y*v2.z+v1.z*v2.y);
		elm[4] = v1.y*v2.y;
		elm[5] = v1.z*v2.z;
	}
	
	inline void
	reset(){
		for (int i=0; i<6; i++){
			elm[i] = 0;
		}
	}
	
	inline friend StressTensor
	tensor_prod (const vec3d &v1, const vec3d &v2) {
		return StressTensor(v1.x*v2.x,
							0.5*(v1.x*v2.y+v1.y*v2.x),
							0.5*(v1.x*v2.z+v1.z*v2.x),
							0.5*(v1.y*v2.z+v1.z*v2.y),
							v1.y*v2.y,
							v1.z*v2.z);
	}
	
	/*	N1 = Sxx-Szz;
	 */
	double getNormalStress1(){
		return elm[0]-elm[5];
	}
	/*	N2 = Szz-Syy;
	 */
	double getNormalStress2(){
		return elm[5]-elm[4];
	}
	
	double getStressXZ(){
		return elm[2];
	}
	
	double getParticlePressure(){
		return -(1./3)*(elm[0]+elm[4]+elm[5]);
	}
	
};
#endif
