/*
 *  vec3d.h
 *  LF_DEM
 *
 *  Created by seto on 09/08/11.
 *  Copyright 2009-2014 Ryohei Seto and Romain Mari. All rights reserved.
 *
 */

/**
 \class vec3d
 \brief 3d vector object.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef vec3d_h
#define vec3d_h 1
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

class vec3d {
public:
	/* variables */
	double x;
	double y;
	double z;
	
	/* constructor/destructor */
	vec3d (): x(0), y(0), z(0){}
	vec3d (double _x,
		   double _y,
		   double _z): x(_x), y(_y), z(_z) {}
	
	vec3d (double a): x(a), y(a), z(a) {}
	
	friend vec3d operator + (const vec3d& v)
	{
		return v;
	}
	
	friend vec3d operator - (const vec3d& v)
	{
		return vec3d(-v.x, -v.y, -v.z);
	}
	
	/* division */
	template <typename T>
	friend vec3d operator / (const vec3d& v,
							 const T& a)
	{
		return vec3d(v.x/a, v.y/a, v.z/a);
	}
	
	/* multiplication */
	template <typename T>
	friend vec3d operator * (const vec3d& v,
							 const T& a)
	{
		return a*v;
	}
	
	template <typename T>
	friend vec3d operator * (const T& a,
							 const vec3d& v)
	{
		return vec3d(a*v.x, a*v.y, a*v.z);
	}
	
	// assign operator
	vec3d& operator += (const vec3d& v)
	{
		x += v.x, y += v.y, z += v.z;
		return *this;
	}
	
	vec3d& operator -= (const vec3d& v)
	{
		x -= v.x, y -= v.y, z -= v.z;
		return *this;
	}
	
	template <typename T>
	vec3d& operator *= (const T& a)
	{
		x *= a, y *= a, z *= a;
		return *this;
	}
	
	template <typename T>
	vec3d& operator /= (const T& a)
	{
		x /= a, y /= a, z /= a;
		return 	*this;
	}
	
	bool is_nonzero()
	{
		if (x != 0 || z != 0 || y != 0) {
			return true;
		} else {
			return false;
		}
	}
	
	void set(double elm_x, double elm_y, double elm_z)
	{
		x = elm_x, y = elm_y, z = elm_z;
	}
	
	void reset()
	{
		x = 0, y = 0, z = 0;
	}
	
	void unitvector()
	{
		(*this) = (*this)/norm();
	}
	
	double sq_norm() const
	{
		return x*x+y*y+z*z;
	}
	
	double sq_norm_xy() const
	{
		return x*x+y*y;
	}
	
	double sq_norm_xz() const
	{
		return x*x+z*z;
	}
	
	double norm() const
	{
		return sqrt(sq_norm());
	}
	
	double norm_xz() const
	{
		return sqrt(sq_norm_xz());
	}
	
	void rotateInfinitesimal(const vec3d& dphi);
	
	void vertical_projection(const vec3d& v);
	
	vec3d product_rate_of_strain(double* E)
	{
		vec3d product(E[0]*x+E[1]*y+E[2]*z,
					  E[1]*x+E[4]*y+E[3]*z,
					  E[2]*x+E[3]*y+(-E[0]-E[4])*z);
		return product;
	}
	
	inline double angle_elevation_xz()
	{
		return atan2(z, x);
	}
	
	inline double angle_0_pi()
	{
		if (x > 0) {
			return atan2(x, z);
		} else {
			return atan2(-x, z);
		}
	}
	
	void periodicBoundaryBox(const double& lx,
							 const double& ly,
							 const double& lz)
	{
		if (x < 0) {
			x += lx;
		} else if (x > lx) {
			x -= lx;
		}
		if (y < 0) {
			y += ly;
		} else if (y > ly) {
			y -= ly;
		}
		if (z < 0) {
			z += lz;
		} else if (z > lz) {
			z -= lz;
		}
	}
	
	void cerr(std::string note)
	{
		std::cerr << note << ' ' << x << ' '<< y << ' ' << z << std::endl;
	}
	
	void cerr()
	{
		std::cerr << x << ' '<< y << ' ' << z << std::endl;
	}
};

/************ vec3d helping functions *******************/

/*** Symmetric binary operators ***/
inline bool operator == (const vec3d &v1,
						 const vec3d &v2)
{
	if (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z) {
		return true;
	}
	return false;
}

inline bool operator != (const vec3d& v1,
						 const vec3d& v2)
{
	if (v1.x != v2.x || v1.y != v2.y || v1.z != v2.z) {
		return true;
	}
	return false;
}

inline vec3d operator + (const vec3d& a1,
						 const vec3d& a2)
{
	return vec3d(a1.x+a2.x, a1.y+a2.y, a1.z+a2.z);
}

/* subtraction */
inline vec3d operator - (const vec3d& a1,
						 const vec3d& a2)
{
	return vec3d(a1.x-a2.x, a1.y-a2.y, a1.z-a2.z);
}

// output stream operator
inline std::ostream& operator << (std::ostream& out,
								  const vec3d& v)
{
	out << v.x << " " << v.y << " " << v.z;
	return out;
}

/******* Other, dist, dot, cross, etc **********/
inline double dist(const vec3d& a1, const vec3d& a2)
{
	return (a1-a2).norm();
}

inline double sq_dist(const vec3d& a1, const vec3d& a2)
{
	return (a1-a2).sq_norm();
}

/* scalar product */
inline double dot(const vec3d& a1,
				  const vec3d& a2)
{
	return a1.x*a2.x+a1.y*a2.y+a1.z*a2.z;
}

inline double dot(const vec3d* a1,
				  const vec3d& a2)
{
	return a1->x*a2.x+a1->y*a2.y+a1->z*a2.z;
}

inline double dot(const vec3d* a1,
				  const vec3d* a2)
{
	return a1->x*a2->x+a1->y*a2->y+a1->z*a2->z;
}

inline double dot(const vec3d& a1,
				  const vec3d* a2)
{
	return a1.x*a2->x+a1.y*a2->y+a1.z*a2->z;
}


/* vector product */
inline vec3d cross(const vec3d& v1,
				   const vec3d& v2)
{
	return vec3d(v1.y*v2.z - v1.z*v2.y,
				 v1.z*v2.x - v1.x*v2.z,
				 v1.x*v2.y - v1.y*v2.x);
}

inline vec3d cross(const vec3d* v1,
				   const vec3d& v2)
{
	return vec3d(v1->y*v2.z - v1->z*v2.y,
				 v1->z*v2.x - v1->x*v2.z,
				 v1->x*v2.y - v1->y*v2.x);
}

inline vec3d cross(const vec3d& v1,
				   const vec3d* v2)
{
	return vec3d(v1.y*v2->z - v1.z*v2->y,
				 v1.z*v2->x - v1.x*v2->z,
				 v1.x*v2->y - v1.y*v2->x);
}

/* vector product */
inline vec3d cross_vec_array(const vec3d& v1,
							 const double* v2p)
{
	return vec3d(v1.y*(*(v2p+2))-v1.z*(*(v2p+1)),
				 v1.z*(*v2p)-v1.x*(*(v2p+2)),
				 v1.x*(*(v2p+1))-v1.y*(*v2p));
}

inline void
vec3d::rotateInfinitesimal(const vec3d& dphi)
{
	/* dphi must be small vector. */
	(*this) += cross(dphi, *this);
}

inline void
vec3d::vertical_projection(const vec3d& v)
{
	(*this) -= dot(*this, v)*v;
}


// inline vec3d str2vec3d(const std::string& value)
// {
// 	std::string::size_type l1 = value.find("(", 0);
// 	if (l1 == std::string::npos) {
// 		exit(1);
// 	}
// 	std::string::size_type l2 = value.find(",", l1);
// 	if (l2 == std::string::npos) {
// 		exit(1);
// 	}
// 	std::string::size_type l3 = value.find(",", l2+1);
// 	if (l3 == std::string::npos) {
// 		exit(1);
// 	}
// 	std::string::size_type l4 = value.find(")", l3+1);
// 	if (l4 == std::string::npos) {
// 		exit(1);
// 	}
// 	double vx = atof(value.substr(l1+1, l2-l1-1).c_str());
// 	double vy = atof(value.substr(l2+1, l3-l2-1).c_str());
// 	double vz = atof(value.substr(l3+1, l4-l3-1).c_str());
// 	return vec3d(vx,vy,vz);
// }
#endif
