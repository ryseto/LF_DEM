#ifndef __LF_DEM__Random__
#define __LF_DEM__Random__

#ifndef USE_DSFMT
#include "MersenneTwister.h"
#endif

#ifndef USE_DSFMT
#define GRANDOM ( r_gen->randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.
#endif
#ifdef USE_DSFMT
#define GRANDOM  ( sqrt( -2.0 * log( 1.0 - dsfmt_genrand_open_open(&r_gen) ) ) * cos(2.0 * 3.14159265358979323846264338328 * dsfmt_genrand_close_open(&r_gen) ) ) // RNG gaussian with mean 0. and variance 1.
#endif
using namespace std;

#ifdef USE_DSFMT
inline unsigned long
wagnerhash(time_t t, clock_t c)
{
	/**
		\brief Utility function to start up the DSFMT RNG with a nice seed.

	 From MersenneTwister v1.0 by Richard J. Wagner
	 comments below are from the original code.

	 Get a unsigned long from t and c
	 Better than unsigned long(x) in case x is floating point in [0,1]
	 Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
	*/

	static unsigned long differ = 0; // guarantee time-based seeds will change
	unsigned long h1 = 0;
	unsigned char *pp = (unsigned char *) &t;
	for (size_t i=0; i<sizeof(t); ++i){
		h1 *= UCHAR_MAX + 2U;
		h1 += pp[i];
	}
	unsigned long h2 = 0;
	pp = (unsigned char *) &c;
	for (size_t j=0; j<sizeof(c); ++j) {
		h2 *= UCHAR_MAX + 2U;
		h2 += pp[j];
	}
	return (h1 + differ++)^h2;
}
#endif

#endif