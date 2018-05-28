#ifndef __LF_DEM__global__
#define __LF_DEM__global__

#include <csignal>
#include <string>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include "vec3d.h"

#define RATE_INDEPENDENT 0
#define RATE_DEPENDENT 1
#define RATE_PROPORTIONAL 2
#define VELOCITY_STRESS 0
#define STRAIN_STRESS 1
#define XF_STRESS 2
#define BROWNIAN_STRESS 3

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

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

#ifdef SIGINT_CATCH
extern volatile sig_atomic_t sig_caught;

inline void sigint_handler(int signum)
{
	if (signum==SIGINT){
		sig_caught = signum;
		std::cerr << "Received a SIGINT" << std::endl;
	}
}
#endif

inline void removeBlank(std::string& str)
{
	str.erase(remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
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


inline std::vector<std::string> splitString(const std::string& str)
{
	std::string stripped_str = str;
	std::size_t brk = 0;

	std::vector<std::string> elements;

	while (stripped_str.length() > 0) {
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

namespace Parameters {
enum class ControlVariable : unsigned {
	rate,
	stress,
	force,
	// viscnb
};
}

inline void checkInFile(std::string filename)
{
	std::ifstream file_import(filename.c_str(), std::ios::in);
	if (!file_import) {
		std::ostringstream error_str;
		error_str  << " File '" << filename << "' not found." << std::endl;
		throw std::runtime_error(error_str.str());
	}
}
#endif /* defined(__LF_DEM__global__) */
