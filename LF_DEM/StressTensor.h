//
//  StressTensor.h
//  LF_DEM
//
//  Created by Ryohei Seto on 4/16/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

/**
 \class StressTensor
 \brief Stress tensor object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef LF_DEM_StressTensor_h
#define LF_DEM_StressTensor_h
#include "Matrix.h"
#include "vec3d.h"
#include <iostream>
#include <iomanip>
#include <vector>

class StressTensor {
public:
	/*
	 * (xx, xy, xz, yz, yy, zz)
	 */
	std::vector <double> elm;

	inline StressTensor(void)
	{
		elm.resize(6, 0);
	}

	inline StressTensor(double a)
	{
		elm.resize(6);
		for (int i=0; i<6; i++) {
			elm[i] = a;
		}
	}

	inline StressTensor(double _xx,
						double _xy,
						double _xz,
						double _yz,
						double _yy,
						double _zz)
	{
		elm.resize(6);
		elm[0] = _xx;
		elm[1] = _xy;
		elm[2] = _xz;
		elm[3] = _yz;
		elm[4] = _yy;
		elm[5] = _zz;
	}

	inline StressTensor(matrix m)
	{
		// 0(0,0) 1(0,1) 2(0,2)
		// 3(1,0) 4(1,1) 5(1,2)
		// 6(2,0) 7(2,1) 8(2,2)
		elm.resize(6);
		elm[0] = m.elm[0];
		elm[1] = m.elm[1];
		elm[2] = m.elm[2];
		elm[3] = m.elm[5];
		elm[4] = m.elm[4];
		elm[5] = m.elm[8];
	}

	// output stream operator
	inline friend std::ostream& operator << (std::ostream& out,
											 const StressTensor& st)
	{
		out << st.elm[0] << ' ';
		out << st.elm[1] << ' ';
		out << st.elm[2] << ' ';
		out << st.elm[3] << ' ';
		out << st.elm[4] << ' ';
		out << st.elm[5] << ' ';
		return out;
	}

	inline StressTensor(const vec3d& v1, const vec3d& v2)
	{
		elm.resize(6);
		elm[0] = v1.x*v2.x; //xx
		elm[1] = 0.5*(v1.x*v2.y+v1.y*v2.x); //xy
		elm[2] = 0.5*(v1.x*v2.z+v1.z*v2.x); //xz
		elm[3] = 0.5*(v1.y*v2.z+v1.z*v2.y); //yz
		elm[4] = v1.y*v2.y; //yy
		elm[5] = v1.z*v2.z; //zz
	}

	inline StressTensor(const vec3d& v)
	{
		elm.resize(6);
		elm[0] = v.x*v.x;
		elm[1] = v.x*v.y;
		elm[2] = v.x*v.z;
		elm[3] = v.y*v.z;
		elm[4] = v.y*v.y;
		elm[5] = v.z*v.z;
	}
	
	inline friend StressTensor operator + (const StressTensor& a1,
										   const StressTensor& a2)
	{
		return StressTensor(a1.elm[0]+a2.elm[0],
							a1.elm[1]+a2.elm[1],
							a1.elm[2]+a2.elm[2],
							a1.elm[3]+a2.elm[3],
							a1.elm[4]+a2.elm[4],
							a1.elm[5]+a2.elm[5]);
	}

	inline friend StressTensor operator + (const StressTensor& s)
	{
		return s;
	}

	/* subtraction */
	inline friend StressTensor operator - (const StressTensor& a1,
										   const StressTensor& a2)
	{
		return StressTensor(a1.elm[0]-a2.elm[0],
							a1.elm[1]-a2.elm[1],
							a1.elm[2]-a2.elm[2],
							a1.elm[3]-a2.elm[3],
							a1.elm[4]-a2.elm[4],
							a1.elm[5]-a2.elm[5]);
	}

	inline friend StressTensor operator - (const StressTensor& s)
	{
		return StressTensor(-s.elm[0],
							-s.elm[1],
							-s.elm[2],
							-s.elm[3],
							-s.elm[4],
							-s.elm[5]);
	}

	inline friend StressTensor operator * (double d,
										   const StressTensor& s)
	{
		return StressTensor(d*s.elm[0],
							d*s.elm[1],
							d*s.elm[2],
							d*s.elm[3],
							d*s.elm[4],
							d*s.elm[5]);
	}

	inline friend StressTensor operator * (const StressTensor& s,
										   double d)
	{
		return StressTensor(d*s.elm[0],
							d*s.elm[1],
							d*s.elm[2],
							d*s.elm[3],
							d*s.elm[4],
							d*s.elm[5]);
	}

	inline friend StressTensor operator * (int i,
										   const StressTensor& s)
	{
		return StressTensor(i*s.elm[0],
							i*s.elm[1],
							i*s.elm[2],
							i*s.elm[3],
							i*s.elm[4],
							i*s.elm[5]);
	}

	inline friend StressTensor operator * (const StressTensor& s,
										   int i)
	{
		return StressTensor(i*s.elm[0],
							i*s.elm[1],
							i*s.elm[2],
							i*s.elm[3],
							i*s.elm[4],
							i*s.elm[5]);
	}

	inline friend StressTensor operator / (double d,
										   const StressTensor& s)
	{
		return StressTensor(s.elm[0]/d,
							s.elm[1]/d,
							s.elm[2]/d,
							s.elm[3]/d,
							s.elm[4]/d,
							s.elm[5]/d);
	}

	inline friend StressTensor operator / (const StressTensor& s,
										   double d)
	{
		return StressTensor(s.elm[0]/d,
							s.elm[1]/d,
							s.elm[2]/d,
							s.elm[3]/d,
							s.elm[4]/d,
							s.elm[5]/d);
	}

	inline friend StressTensor operator / (int i,
										   const StressTensor& s)
	{
		return StressTensor(s.elm[0]/i,
							s.elm[1]/i,
							s.elm[2]/i,
							s.elm[3]/i,
							s.elm[4]/i,
							s.elm[5]/i);
	}

	inline friend StressTensor operator / (const StressTensor& s,
										   int i)
	{
		return StressTensor(s.elm[0]/i,
							s.elm[1]/i,
							s.elm[2]/i,
							s.elm[3]/i,
							s.elm[4]/i,
							s.elm[5]/i);
	}

	inline StressTensor& operator += (const StressTensor& s)
	{
		for (int i=0; i<6; i++) {
			elm[i] += s.elm[i];
		}
		return *this;
	}

	inline StressTensor& operator -= (const StressTensor& s)
	{
		for (int i=0; i<6; i++) {
			elm[i] -= s.elm[i];
		}
		return *this;
	}

	inline StressTensor& operator *= (double d)
	{
		for (int i=0; i<6; i++) {
			elm[i] *= d;
		}
		return *this;
	}

	inline StressTensor& operator *= (int i)
	{
		for (int j=0; j<6; j++) {
			elm[j] *= i;
		}
		return *this;
	}

	inline StressTensor& operator /= (double d)
	{
		double d_inv = 1.0/d;
		for (int i=0; i<6; i++) {
			elm[i] *= d_inv;
		}
		return *this;
	}

	inline StressTensor& operator /= (int i)
	{
		double d_inv = 1.0/i;
		for (int j=0; j<6; j++) {
			elm[j] *= d_inv;
		}
		return *this;
	}

	/*
	 * Finite trace
	 */
	inline void set(double _xx, double _xy, double _xz,
					double _yz, double _yy, double _zz)
	{
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
	inline void set(double _xx, double _xy, double _xz,
					double _yz, double _yy)
	{
		elm[0] = _xx;
		elm[1] = _xy;
		elm[2] = _xz;
		elm[3] = _yz;
		elm[4] = _yy;
		elm[5] = -(_xx+_yy);
	}

	inline void set(const vec3d& v1, const vec3d& v2)
	{
		/* XF Stress tensor:
		 * S = (r^{(j)}-r^{(i)}) F^{(i,j)}
		 */
		elm[0] = v1.x*v2.x;
		elm[1] = 0.5*(v1.x*v2.y+v1.y*v2.x);
		elm[2] = 0.5*(v1.x*v2.z+v1.z*v2.x);
		elm[3] = 0.5*(v1.y*v2.z+v1.z*v2.y);
		elm[4] = v1.y*v2.y;
		elm[5] = v1.z*v2.z;
	}

	inline void reset() {
		for (int i=0; i<6; i++) {
			elm[i] = 0;
		}
	}

	/*	N1 = Sxx-Szz;
	 */
	double getNormalStress1()
	{
		return elm[0]-elm[5];
	}

	/*	N2 = Szz-Syy;
	 */
	double getNormalStress2()
	{
		return elm[5]-elm[4];
	}

	double getParticlePressure()
	{
		return -(1./3)*(elm[0]+elm[4]+elm[5]);
	}

	matrix getMatrix()
	{
		matrix m(elm[0], elm[1], elm[2],
				 elm[1], elm[4], elm[3],
				 elm[2], elm[3], elm[5]);
		return m;
	}
	//	void cerr()
	//	{
	//		std::cerr << elm[0] << ' ' << elm[1] << ' '<< elm[2] << ' '<< elm[3] << ' '<< elm[4] << ' ' << elm[5] << std::endl;
	//	}
};

// Helper functions
inline double shearStressComponent(const StressTensor& s, double theta_shear)
{
	return cos(theta_shear)*s.elm[2]+sin(theta_shear)*s.elm[3];
}

#endif
