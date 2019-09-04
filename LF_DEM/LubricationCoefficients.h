//
//  LubricationFunctions.h
//  LF_DEM
//
//  Created by Ryohei Seto on 8/27/13.
//  Copyright (c) 2013 Ryohei Seto. All rights reserved.
//

#ifndef LF_DEM_LubricationCoefficients_h
#define LF_DEM_LubricationCoefficients_h


namespace Interactions
{
namespace Lub
{

struct LubricationCoefficients
{
	/* 
		Separation independent coefficients cX* and cY*.

		X_A = cXA*(1/h_ij)
		Y_A = cYA*log(1/h_ij)
		etc
	*/
	double cXA[4]; // ii ij ji jj
	double cYA[4];
	double cYB[4];
	double cYC[4];
	double cXG[4];
	double cYG[4];
	double cYH[4];
	double cXM[4];
	double cYM[4];
	LubricationCoefficients(double r0, double r1);
};

}
}

#endif
