//
//  Matrix.h
//  LF_DEM
//
//  Created by Ryohei Seto on 2016/10/28.
//  Copyright © 2016 Ryohei Seto. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h
#include <vector>
#include "vec3d.h"

class Sym2Tensor;

class matrix {
public:
	/* variables */
	std::vector <double> elm;
	matrix (void)
	{
		elm.resize(9, 0);
	}
	
	matrix (double _e00, double _e01, double _e02,
			double _e10, double _e11, double _e12,
			double _e20, double _e21, double _e22)
	{
		elm.resize(9, 0);
		
		elm[0] = _e00;
		elm[1] = _e01;
		elm[2] = _e02;
		elm[3] = _e10;
		elm[4] = _e11;
		elm[5] = _e12;
		elm[6] = _e20;
		elm[7] = _e21;
		elm[8] = _e22;
	}
	
	double get(int i, int j)
	{
		return elm[3*i+j];
	}
	
	void set(int i, int j, double value)
	{
		// 0(0,0) 1(0,1) 2(0,2)
		// 3(1,0) 4(1,1) 5(1,2)
		// 6(2,0) 7(2,1) 8(2,2)
		elm[3*i+j] = value;
	}
	
	void setSymmetric(double e_00, double e_01, double e_02,
							 double e_12, double e_11, double e_22)
	{
		// 0(0,0) 1(0,1) 2(0,2)
		// 3(1,0) 4(1,1) 5(1,2)
		// 6(2,0) 7(2,1) 8(2,2)
		elm[0] = e_00;
		elm[1] = elm[3] = e_01;
		elm[2] = elm[6] = e_02;
		elm[4] = e_11;
		elm[5] = elm[7] = e_12;
		elm[8] = e_22;
	}
	
	void set_zero()
	{
		for (auto &elm_: elm) {
			elm_ = 0;
		}
	}
	
	void set_rotation(double angle, char axis)
	{
		if (axis == 'x') {
			elm[0] = 1, elm[1] = 0, elm[2] = 0;
			elm[3] = 0, elm[4] = cos(angle), elm[5] = -sin(angle);
			elm[6] = 0, elm[7] = sin(angle), elm[8] = cos(angle);
		} else if (axis == 'y') {
			elm[0] = cos(angle), elm[1] = 0, elm[2] = sin(angle);
			elm[3] = 0, elm[4] = 1, elm[5] = 0;
			elm[6] = -sin(angle), elm[7] = 0, elm[8] = cos(angle);
		} else if (axis == 'z') {
			elm[0] = cos(angle), elm[1] = -sin(angle), elm[2] = 0;
			elm[3] = sin(angle), elm[4] = cos(angle), elm[5] = 0;
			elm[6] = 0, elm[7] = 0, elm[8] = 1;
		} else {
			exit(1);
		}
	}
	
	matrix inverse()
	{
		// 00 11 22   048
		// 10 21 02   372
		// 20 01 12   615
		
		// 00 21 12   075
		// 20 11 02   742
		// 10 01 22   318
		
		double determinant = elm[0]*elm[4]*elm[8] + elm[3]*elm[7]*elm[2] + elm[6]*elm[1]*elm[5] \
		- elm[0]*elm[7]*elm[5] - elm[6]*elm[4]*elm[2] - elm[3]*elm[1]*elm[8];
		matrix m_inv;
		m_inv.elm[0] = (elm[4]*elm[8]-elm[5]*elm[7])/determinant;
		m_inv.elm[1] = (elm[2]*elm[7]-elm[1]*elm[8])/determinant;
		m_inv.elm[2] = (elm[1]*elm[5]-elm[2]*elm[4])/determinant;
		m_inv.elm[3] = (elm[6]*elm[5]-elm[3]*elm[8])/determinant;
		m_inv.elm[4] = (elm[0]*elm[8]-elm[6]*elm[2])/determinant;
		m_inv.elm[5] = (elm[2]*elm[3]-elm[0]*elm[5])/determinant;
		m_inv.elm[6] = (elm[3]*elm[7]-elm[4]*elm[6])/determinant;
		m_inv.elm[7] = (elm[1]*elm[6]-elm[0]*elm[7])/determinant;
		m_inv.elm[8] = (elm[0]*elm[4]-elm[3]*elm[1])/determinant;
		return m_inv;
	}
	
	matrix antiSymmetrise()
	{
		matrix m_unsym;
		for (int i=0;i<3;i++) {
			for (int j=0; j<3; j++) {
				m_unsym.elm[3*i+j] = (elm[3*i+j]-elm[3*j+i])/2;
			}
		}
		return m_unsym;
	}

	vec3d getLine(int i)
	{
		// 0(0,0) 1(0,1) 2(0,2)
		// 3(1,0) 4(1,1) 5(1,2)
		// 6(2,0) 7(2,1) 8(2,2)
		return vec3d(elm[3*i], elm[3*i+1], elm[3*i+2]);
	}
	
	vec3d getColumn(int j)
	{
		// 0(0,0) 1(0,1) 2(0,2)
		// 3(1,0) 4(1,1) 5(1,2)
		// 6(2,0) 7(2,1) 8(2,2)
		return vec3d(elm[j], elm[3+j], elm[6+j]);
	}
	
	friend matrix operator * (const matrix& m,
							  const double& a)
	{
		return a*m;
	}
	
	/* multiplication */
	friend matrix operator * (const double& a,
							  const matrix& m)
	{
		return matrix(a*m.elm[0], a*m.elm[1], a*m.elm[2],
					  a*m.elm[3], a*m.elm[4], a*m.elm[5],
					  a*m.elm[6], a*m.elm[7], a*m.elm[8]);
	}
	
	friend vec3d operator * (const matrix& m,
							 const vec3d& v)
	{
		return vec3d(m.elm[0]*v.x+m.elm[1]*v.y+m.elm[2]*v.z,
					 m.elm[3]*v.x+m.elm[4]*v.y+m.elm[5]*v.z,
					 m.elm[6]*v.x+m.elm[7]*v.y+m.elm[8]*v.z);
		
	}
	
	friend vec3d operator * (const vec3d& v,
							 const matrix& m)
	{
		return vec3d(m.elm[0]*v.x+m.elm[3]*v.y+m.elm[6]*v.z,
					 m.elm[1]*v.x+m.elm[4]*v.y+m.elm[7]*v.z,
					 m.elm[2]*v.x+m.elm[5]*v.y+m.elm[8]*v.z);
		
	}
	
	friend matrix operator * (const matrix& m1,
							  const matrix& m2)
	{
		matrix product;
		for (int i=0; i<3 ;i++) {
			for (int j=0; j<3 ;j++) {
				for (int k=0; k<3 ;k++) {
					product.elm[3*i+j] += m1.elm[3*i+k]*m2.elm[3*k+j];
				}
			}
		}
		return product;
	}
	
	friend matrix operator + (const matrix& m1,
							  const matrix& m2)
	{
		matrix sum;
		for (int i=0; i<9; i++) {
			sum.elm[i] = m1.elm[i]+m2.elm[i];
		}
		return sum;
	}
	
	friend matrix operator + (const double a,
							  const matrix& m)
	{
		return matrix(a+m.elm[0], m.elm[1], m.elm[2],
					  m.elm[3], a+m.elm[4], m.elm[5],
					  m.elm[6], m.elm[7], a+m.elm[8]);
	}
	matrix transpose() const
	{
		return matrix (elm[0], elm[3], elm[6],
					   elm[1], elm[4], elm[7],
					   elm[2], elm[5], elm[8]);
	}
	void print() const
	{
		std::cout << std::setw(10) << elm[0] << std::setw(10) << elm[1] << std::setw(10) <<  elm[2] << std::endl;
		std::cout << std::setw(10) << elm[3] << std::setw(10) << elm[4] << std::setw(10) <<  elm[5] << std::endl;
		std::cout << std::setw(10) << elm[6] << std::setw(10) << elm[7] << std::setw(10) <<  elm[8] << std::endl;
	}
};

std::pair<Sym2Tensor, vec3d> symmetryDecomposition(const matrix &m);
matrix vec2AntiSymFullMatrix(vec3d omega);
matrix sym2FullMatrix(Sym2Tensor sym);

#endif /* Matrix_h */
