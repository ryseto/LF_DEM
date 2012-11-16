/*
 *  vec3d.h
 *  CCN_3D
 *
 *  Created by seto on 09/08/11.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef vec3d_h
#define vec3d_h 1
#include <iostream>
#include <iomanip>
#include <cmath>

class vec3d {
public:
	/* variables */
	double x, y, z;

	/* constructor/destructor */
	inline vec3d (void){}
	inline vec3d (const double &_x, 
								const double &_y, 
								const double &_z): x(_x), y(_y), z(_z) {}
	inline ~vec3d(void){}
		
	/* operators */
	inline vec3d& operator = (const vec3d& v){
		x = v.x, y = v.y, z = v.z;
		return *this;
	}
	
	inline friend bool operator==(const vec3d &v1, const vec3d &v2){
		if (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z)
			return true;
		return false;   
	}
	inline friend bool operator!=(const vec3d &v1, const vec3d &v2){
		if (v1.x != v2.x || v1.y != v2.y || v1.z != v2.z )
			return true;
		return false;
	}
	
	inline friend vec3d operator + (const vec3d &a1, const vec3d &a2){
		return vec3d( a1.x + a2.x, a1.y + a2.y, a1.z + a2.z);
	}
	inline friend vec3d operator + (const vec3d &v){
		return v;
	}
	
	/* subtraction */
	inline friend vec3d operator - (const vec3d &a1, const vec3d &a2){
		return vec3d( a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
	}
	inline friend vec3d operator - (const vec3d &v){
		return vec3d( -v.x, -v.y, -v.z);
	}
	
	/* multiplication */
	inline friend vec3d operator * (const double &d, const vec3d &v){
		return vec3d( d*v.x, d*v.y, d*v.z);
	}
	inline friend vec3d operator * (const vec3d &v, const double &d){
		return d*v;
	}
	inline friend vec3d operator * (const int &i, const vec3d &v){
		double d_tmp = (double)i; 
		return d_tmp*v;
	}
	inline friend vec3d operator * (const vec3d &v, const int &i){
		double d_tmp = (double)i; 
		return d_tmp*v;
	}
	
	/* scalar product */
	inline friend double dot(const vec3d &a1, const vec3d &a2){
	return a1.x*a2.x + a1.y*a2.y + a1.z*a2.z;
	}
	
	
	/* vector product */
	inline friend vec3d cross(const vec3d &v1, const vec3d &v2){
		return vec3d(v1.y*v2.z - v1.z*v2.y, 
                     v1.z*v2.x - v1.x*v2.z, 
                     v1.x*v2.y - v1.y*v2.x);
	}	

    /* vector product */
	inline friend vec3d cross_vec_array(const vec3d &v1, 
                                        const double *v2p){
		return vec3d(v1.y*(*(v2p+2)) - v1.z*(*(v2p+1)), 
                     v1.z*(*v2p) - v1.x*(*(v2p+2)), 
                     v1.x*(*(v2p+1)) - v1.y*(*v2p));
	}	

    
	/* division */
	inline friend vec3d operator / (const vec3d &v, const double &d){
		double d_tmp = 1.0/d;
		return d_tmp*v;
	}
	inline friend vec3d operator / (const vec3d &v, const int &i){
		double d_tmp = 1.0/i;
		return d_tmp*v;
	}

	// assign operator
	inline vec3d& operator +=(const vec3d &v){
		x += v.x, y += v.y, z += v.z;
		return *this;
	}
	inline vec3d& operator -=(const vec3d &v){
		x -= v.x, y -= v.y, z -= v.z;
		return *this;
	}
	inline vec3d& operator *=(const double & d){
		x *= d, y *= d, z *= d;
		return *this;
	}
	inline vec3d& operator *=(const int & i){
		double d_tmp = (double)i;
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
		return *this;
	}
	inline vec3d& operator /=(const double & d){
		double d_tmp = 1.0/d;
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
		return 	*this;
	}
	inline vec3d& operator /=(const int & i){
		double d_tmp = 1.0/i;
		x *= d_tmp, y *= d_tmp, z *= d_tmp;
		return *this;
	}

	/* utility */
	inline void set(const double &_x, const double &_y, const double &_z){x=_x, y=_y, z=_z;}
	inline void reset(){ x = 0, y = 0, z = 0;}
    inline void add(const double &_dx, const double &_dy, const double &_dz){
        x += _dx, y += + _dy, z += _dz;}
	inline void unitvector(){(*this) = (*this) / norm();} 
	inline void sign_reverse(){ (*this) = -(*this);}
	inline double sq_norm(){return x*x + y*y + z*z; }
	inline double sq_norm_xy(){return x*x + y*y;}
	inline double norm(){return sqrt(sq_norm());}

	inline friend double dist(const vec3d &a1, const vec3d &a2){
		return (a1 - a2).norm();}	
	inline friend double sq_dist(const vec3d &a1, const vec3d &a2){		
		return (a1 - a2).sq_norm() ;
	}
	inline void rotateInfinitesimal(const vec3d &dphi){
		/* dphi must be small vector. */
		(*this) += cross(dphi, *this);
	}
	inline vec3d product_rate_of_strain(double * E){
		vec3d product(E[0]*x + E[1]*y + E[2]*z ,
				  E[1]*x + E[4]*y + E[3]*z ,
				  E[2]*x + E[3]*y + (-E[0]-E[4])*z );
		
		return product;
	}
	void cerr(){
		std::cerr << x << ' '<< y << ' ' << z << std::endl;
	}	
};
#endif	
