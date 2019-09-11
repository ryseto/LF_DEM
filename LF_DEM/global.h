#ifndef __LF_DEM__global__
#define __LF_DEM__global__

#include <csignal>
#include <string>
#include <iostream>

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

#endif /* defined(__LF_DEM__global__) */
