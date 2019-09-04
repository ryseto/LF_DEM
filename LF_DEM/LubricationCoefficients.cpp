#include "LubricationCoefficients.h"

double func_g1_XA(double lamb)
{
	double lamb_p_1 = lamb+1;
	return 2*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g2_YA(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (4./15)*lamb*(2+lamb+2*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g2_YB(double lamb)
{
	double lamb_p_1 = lamb+1;
	return -(1./5)*lamb*(4+lamb)/(lamb_p_1*lamb_p_1);
}

/*
 * XC ---> no singular contribution
 */
double func_g2_YC(double lamb)
{
	return (2./5)*lamb/(lamb+1);
}

double func_g4_YC(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (4./5)*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g1_XG(double lamb)
{
	double lamb_p_1 = lamb+1;
	return 3*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g2_YG(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (1./10)*lamb*(4-lamb+7*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g2_YH(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (1./10)*lamb*(2-lamb)/(lamb_p_1*lamb_p_1);
}

double func_g5_YH(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (1./20)*lamb*lamb*(1+7*lamb)/(lamb_p_1*lamb_p_1);
}

double func_g1_XM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (6./5)*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g4_XM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (6./5)*lamb*lamb*lamb/(lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g2_YM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (6./25)*lamb*(1-lamb+4*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

double func_g5_YM(double lamb)
{
	double lamb_p_1 = lamb+1;
	return (3./50)*lamb*lamb*(7-10*lamb+7*lamb*lamb)/(lamb_p_1*lamb_p_1*lamb_p_1);
}

namespace Interactions {

namespace Lub {

LubricationCoefficients::LubricationCoefficients(double a0, double a1)
{
	auto ro = a0+a1;
	auto ro_12 = ro/2;
	auto lambda = a1/a0;
	auto invlambda = 1/lambda;
	auto lambda_square = lambda*lambda;
	auto lambda_cubic = lambda*lambda*lambda;
	auto lambda_p_1 = 1+lambda;
	auto lambda_p_1_square = lambda_p_1*lambda_p_1;
	auto lambda_p_1_cubic = lambda_p_1_square*lambda_p_1;
	auto a0a0_23 = a0*a0*(2.0/3);
	auto a1a1_23 = a1*a1*(2.0/3);
	auto roro_16 = ro*ro*(1.0/6);
	auto a0a0a0 = a0*a0*a0;
	auto a1a1a1 = a1*a1*a1;
	auto rororo = ro*ro*ro;
	auto a0a0a0_43 = a0a0a0*(4.0/3);
	auto a1a1a1_43 = a1a1a1*(4.0/3);
	auto rororo_16 = rororo*(1.0/6);
	auto a0a0a0_109 = a0a0a0*(10.0/9);
	auto a1a1a1_109 = a1a1a1*(10.0/9);
	auto rororo_536 = rororo*(5.0/36);
	/* XA
	 * X_{a,b}(l) = X_{b,a}(l) = X_{3-a,3-b}(1/l)
	 * X21(l) = X12(l)
	 * X22(l) = X11(1/l)
	 */
	auto g1_XA     = func_g1_XA(lambda);
	auto g1_inv_XA = func_g1_XA(invlambda);
	cXA[0] = a0*g1_XA;
	cXA[1] = ro_12*(-2/lambda_p_1)*g1_XA;
	cXA[2] = cXA[1];
	cXA[3] = a1*g1_inv_XA;
	/* YA
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	auto g2_YA     = func_g2_YA(lambda);
	auto g2_inv_YA = func_g2_YA(invlambda);
	cYA[0] = a0*g2_YA;
	cYA[1] = ro_12*(-2/lambda_p_1)*g2_YA;
	cYA[2] = cYA[1];
	cYA[3] = a1*g2_inv_YA;
	/* YB
	 * Y_{a,b}(l) = -Y_{3-a,3-b}(1/l)
	 * Y21(l) = -Y12(1/l)
	 * Y22(l) = -Y11(1/l)
	 */
	auto g2_YB     = func_g2_YB(lambda);
	auto g2_inv_YB = func_g2_YB(invlambda);
	cYB[0] = a0a0_23*g2_YB;
	cYB[1] = roro_16*(-4/lambda_p_1_square)*g2_YB;
	cYB[2] = roro_16*(4*lambda_square/lambda_p_1_square)*g2_inv_YB;
	cYB[3] = -a1a1_23*g2_inv_YB;
	/* YC
	 * Y_{a,b}(l) = Y_{b,a}(l) = Y_{3-a,3-b}(1/l})
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	auto g2_YC     = func_g2_YC(lambda);
	auto g2_inv_YC = func_g2_YC(invlambda);
	auto g4_YC     = func_g4_YC(lambda);
	cYC[0] = a0a0a0_43*g2_YC;
	cYC[1] = rororo_16*g4_YC;
	cYC[2] = cYC[1];
	cYC[3] = a1a1a1_43*g2_inv_YC;
	/* XG
	 * X_{a,b}(l) = -X_{3-a,3-b}(1/l)
	 * X21(l) = -X12(1/l)
	 * X22(l) = -X11(1/l)
	 */
	auto g1_XG     = func_g1_XG(lambda);
	auto g1_inv_XG = func_g1_XG(invlambda);
	cXG[0] = a0a0_23*g1_XG;
	cXG[1] = roro_16*(-4/lambda_p_1_square)*g1_XG;
	cXG[2] = roro_16*(4*lambda_square/lambda_p_1_square)*g1_inv_XG;
	cXG[3] = -a1a1_23*g1_inv_XG;
	/* YG
	 * Y_{a,b}(l) = -Y_{3-a,3-b}(1/l)
	 * Y21(l) = -Y12(1/l)
	 * Y22(l) = -Y11(1/l)
	 */
	auto g2_YG     = func_g2_YG(lambda);
	auto g2_inv_YG = func_g2_YG(invlambda);
	cYG[0] = a0a0_23*g2_YG;
	cYG[1] = roro_16*(-4/lambda_p_1_square)*g2_YG;
	cYG[2] = roro_16*(4*lambda_square/lambda_p_1_square)*g2_inv_YG;
	cYG[3] = -a1a1_23*g2_inv_YG;
	/* YH
	 * Y_{a,b}(l) = Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(1/l)
	 * Y22(l) = Y11(1/l)
	 */
	auto g2_YH     = func_g2_YH(lambda);
	auto g2_inv_YH = func_g2_YH(invlambda);
	auto g5_YH     = func_g5_YH(lambda);
	auto g5_inv_YH = func_g5_YH(invlambda);
	cYH[0] = a0a0a0_43*g2_YH;
	cYH[1] = rororo_16*(8/lambda_p_1_cubic)*g5_YH;
	cYH[2] = rororo_16*(8*lambda_cubic/lambda_p_1_cubic)*g5_inv_YH;
	cYH[3] = a1a1a1_43*g2_inv_YH;
	/* XM
	 * X_{a,b}(l) = X_{b,a}(l)= X_{3-a,3-b}(1/l)
	 * X21(l) = X12(l)
	 * X22(l) = X11(1/l)
	 */
	auto g1_XM     = func_g1_XM(lambda);
	auto g1_inv_XM = func_g1_XM(invlambda);
	auto g4_XM     = func_g4_XM(lambda);
	cXM[0] = a0a0a0_109*g1_XM;
	cXM[1] = rororo_536*(8/lambda_p_1_cubic)*g4_XM;
	cXM[2] = cXM[1];
	cXM[3] = a1a1a1_109*g1_inv_XM;
	/* YM
	 * Y_{a,b}(l) = Y_{b,a}(l)= Y_{3-a,3-b}(1/l)
	 * Y21(l) = Y12(l)
	 * Y22(l) = Y11(1/l)
	 */
	auto g2_YM     = func_g2_YM(lambda);
	auto g2_inv_YM = func_g2_YM(invlambda);
	auto g5_YM     = func_g5_YM(lambda);
	cYM[0] = a0a0a0_109*g2_YM;
	cYM[1] = rororo_536*(8/lambda_p_1_cubic)*g5_YM;
	cYM[2] = cYM[1];
	cYM[3] = a1a1a1_109*g2_inv_YM;
}

} // namespace Lub

} // namespace Interactions