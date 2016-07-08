//
//  LubricationFunctions.h
//  LF_DEM
//
//  Created by Ryohei Seto on 8/27/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef LF_DEM_LubricationFunctions_h
#define LF_DEM_LubricationFunctions_h

inline double func_g1_XA(double lamb)
{
	double lamb_p_1 = lamb+1;
	return 2*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g2_YA(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (4./15)*lamb*(2+lamb+2*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g2_YB(double lamb)
{
	double lamb_p_1 = lamb+1;
	return -(1./5)*lamb*(4+lamb)/(lamb_p_1*lamb_p_1);
}

/*
 * XC ---> no singular contribution
 */
inline double func_g2_YC(double lamb)
{
	return (2./5)*lamb/(lamb+1);
}

inline double func_g4_YC(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (4./5)*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g1_XG(double lamb)
{
	double lamb_p_1 = lamb+1;
	return 3*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g2_YG(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (1./10)*lamb*(4-lamb+7*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g2_YH(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (1./10)*lamb*(2-lamb)/(lamb_p_1*lamb_p_1);
}

inline double func_g5_YH(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (1./20)*lamb*lamb*(1+7*lamb)/(lamb_p_1*lamb_p_1);
}

inline double func_g1_XM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (6./5)*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g4_XM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (6./5)*lamb*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g2_YM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (6./25)*lamb*(1-lamb+4*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

inline double func_g5_YM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (3./50)*lamb*lamb*(7-10*lamb+7*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

#endif
