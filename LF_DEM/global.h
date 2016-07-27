#ifndef __LF_DEM__global__
#define __LF_DEM__global__

#include <string>
#include <algorithm>
#include <ostream>
#include <vector>
#include "vec3d.h"

#ifndef GIT_VERSION
/*
 * GIT_VERSION stores the output of `git describe --dirty --always`
 * run on the last compilation *made in a git repo*.
 *
 * When in a git repo:
 * - Makefile defines GIT_VERSION := $(shell git describe --dirty --always)
 * - In Xcode, VersionInfo.h is automatically generated
 *   by the following script in the Pre-Action of Build.
 *   -----------------------------------
 *   git=/usr/bin/git
 *   cd "${PROJECT_DIR}/${TARGETNAME}"
 *   version=`$git describe --dirty`
 *   echo "#define GIT_VERSION \"$version\"" > VersionInfo.h
 *   cp -p "${PROJECT_DIR}/tools/circularrheo.pl" ~/bin/
 *   cp -p "${PROJECT_DIR}/tools/analyzeCircularWideGap.pl" ~/bin/
 *   -----------------------------------
 *   cp ${BUILT_PRODUCTS_DIR}/${TARGETNAME} $HOME/bin/
 *   -----------------------------------
 *   and the code uses VersionInfo.h
 *
 * When NOT in a git repo, the code tries to import VersionInfo.h.
 * VersionInfo.h is created upon creation of a source tarball with `make tar` from a git repo.
 *
 */
#include "VersionInfo.h"
#endif

#ifdef USE_DSFMT
#include <limits.h>
#endif

inline void removeBlank(std::string& str)
{
	str.erase(remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
}

inline bool getSuffix(const std::string& str, std::string& value, std::string& suffix)
{
	size_t suffix_pos = str.find_first_of("abcdfghijklmnopqrstuvwxyz"); // omission of "e" is intended, to allow for scientific notation like "1e5h"
	value = str.substr(0, suffix_pos);
	if (suffix_pos != str.npos) {
		suffix = str.substr(suffix_pos, str.length());
		return true;
	} else {
		return false;
	}
}

inline void errorNoSuffix(std::string quantity)
{
	std::cerr << "Error : no unit scale (suffix) provided for " << quantity << std::endl; exit(1);
}

inline bool str2bool(const std::string& value)
{
	if (value == "true") {
		return true;
	} else if (value == "false") {
		return false;
	} else {
		std::cerr << "The value should be true or false" << std::endl;
		exit(1);
	}
}

inline vec3d str2vec3d(const std::string& value)
{
	std::string::size_type l1 = value.find("(", 0);
	if (l1 == std::string::npos) {
		exit(1);
	}
	std::string::size_type l2 = value.find(",", l1);
	if (l2 == std::string::npos) {
		exit(1);
	}
	std::string::size_type l3 = value.find(",", l2+1);
	if (l3 == std::string::npos) {
		exit(1);
	}
	std::string::size_type l4 = value.find(")", l3+1);
	if (l4 == std::string::npos) {
		exit(1);
	}
	double vx = atof(value.substr(l1+1, l2-l1-1).c_str());
	double vy = atof(value.substr(l2+1, l3-l2-1).c_str());
	double vz = atof(value.substr(l3+1, l4-l3-1).c_str());
	return vec3d(vx,vy,vz);
}

inline void Str2KeyValue(const std::string& str_parameter,
						 std::string& keyword,
						 std::string& value)
{
	std::string::size_type pos_equal = str_parameter.find("=");
	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
	return;
}

inline std::vector<std::string> splitString(const std::string& str){
	std::string stripped_str = str;
	std::size_t brk = 0;

	std::vector<std::string> elements;

	while ( stripped_str.length()>0 ) {
		brk = stripped_str.find(" ");
		std::string first_part;
		if (brk < std::string::npos) {
			brk += 1;
			first_part = stripped_str.substr(0, brk);
			stripped_str = stripped_str.substr(brk, std::string::npos);
		} else {
			first_part = stripped_str.substr(0, brk);
			stripped_str = "";
		}
		removeBlank(first_part);
		if ( first_part.length()>0 ) {
			elements.push_back(first_part);
		}
	}
	return elements;
}

#ifdef USE_DSFMT
inline unsigned long
hash(time_t t, clock_t c)
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

#endif /* defined(__LF_DEM__global__) */
