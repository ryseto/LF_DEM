#Makefile LF_DEM, Romain Mari

#=======================================#
#    Parameters to provide              #
#=======================================#


# OS_version: OSX or Linux or SI
OS_version = OSX
# install directory
install_dir = ~/bin/

DSFMT_RNG = no

# can define your own variables in a Makefile.in without having to commit every time, e.g. "OS_version = Linux"
mymake = MyMakefile.in
ifneq ("$(wildcard $(mymake))","")
include $(mymake)
endif




# Compiler
ifeq ($(OS_version),SI)
	CXX=icc
endif
ifeq ($(OS_version),OIST)
	CXX=icc
endif
ifeq ($(OS_version),Linux)
	CXX=g++
endif
ifeq ($(OS_version),OSX)
	CXX=clang++
endif

CXXFLAGS_optimized = -O3
CXXFLAGS_debug = -g -DDEV="1"
CXXFLAGS_dev = -O3 -DDEV="1"

CXXFLAGS_optimized += -std=c++11
CXXFLAGS_debug += -std=c++11
CXXFLAGS_dev += -std=c++11

GPP_Warnings = -Wall
ifeq ($(CXX),$(filter $(CXX),g++ clang++))
	CXXFLAGS_optimized += $(GPP_Warnings)
	CXXFLAGS_debug += -pg $(GPP_Warnings) -Wextra -Wshadow
	CXXFLAGS_dev += $(GPP_Warnings) -Wextra -Wshadow
endif

ICC_Warnings = -Wall
ifeq ($(CXX),icc)
	CXXFLAGS_optimized += $(ICC_Warnings)
	CXXFLAGS_debug += -pg $(ICC_Warnings)
	CXXFLAGS_dev += $(ICC_Warnings)
endif

is_git_repo := $(shell if [ -d .git ] || git rev-parse --git-dir > /dev/null 2>&1; then echo true;fi)
ifeq ($(is_git_repo),true)
	VERSION := $(shell git describe --dirty --always)
	CXXFLAGS_optimized += -DGIT_VERSION="\"$(VERSION)\""
	CXXFLAGS_debug += -DGIT_VERSION="\"$(VERSION)\""
	CXXFLAGS_dev += -DGIT_VERSION="\"$(VERSION)\""
endif

ifeq ($(OS_version),OIST)
	VERSION := $(shell git describe --dirty --always)
	CXXFLAGS_optimized += -DVERSION="\"$(VERSION)\""
	CXXFLAGS_debug += -DVERSION="\"$(VERSION)\""
endif

ifeq ($(OS_version),Linux)
	CXXFLAGS_debug += -pg
endif

ifeq ($(DSFMT_RNG),yes)
	CXXFLAGS_optimized += -DUSE_DSFMT="$DSFMT_RNG" -DDSFMT_MEXP=19937
	CXXFLAGS_debug += -DUSE_DSFMT="$DSFMT_RNG" -DDSFMT_MEXP=19937
	DSFMT_SRC = dSFMT-src-2.2.3/dSFMT.c
endif

# Enable use of Metis library ( yes / no )
UseMetis=no

#======== Include Paths =============
Cholmod_include_path_OSX     = /usr/local/include/
Cholmod_include_path_Linux   = /usr/include/suitesparse/
Cholmod_include_path_SI  = ~/usr/include/
Cholmod_include_path_OIST  = $(SUITESPARSE_ROOT)/include/

Cholmod_path = $(Cholmod_include_path_$(OS_version))
Metis_path = /opt/local/lib/

# Includes
Includes = -I $(Cholmod_path)

#======== Linking ==================
# Blas and Lapack
Blas_Linking_Flags = -lblas
Lapack_Linking_Flags = -llapack

ifeq ($(OS_version), SI)
	Blas_Linking_Flags = -mkl -lrt
	Lapack_Linking_Flags =
endif

ifeq ($(OS_version), OIST)
	Blas_Linking_Flags = -mkl -lrt
	Lapack_Linking_Flags =
endif

# Cholmod

Cholmod_Linking_Flags =
ifeq ($(OS_version), Linux)
Cholmod_Linking_Flags = -lcholmod -lamd -lcamd -lcolamd -lccolamd
endif

ifeq ($(OS_version), SI)
    Cholmod_Linking_Flags  =
endif

ifeq ($(OS_version), OIST)
	Cholmod_Linking_Flags  =
	ifeq ($(UseMetis),yes)
		Cholmod_Linking_Flags += -lmetis
	endif
endif

# extra linking flags for cholmod on Mac OSX
ifeq ($(OS_version),OSX)
	Cholmod_lib_path=$(Cholmod_path:include/=lib/)
	Cholmod_Linking_Flags += $(Cholmod_lib_path)/libcholmod.a $(Cholmod_lib_path)/libamd.a $(Cholmod_lib_path)/libcolamd.a $(Cholmod_lib_path)/libsuitesparseconfig.a
	ifeq ($(UseMetis),yes)
		Cholmod_Linking_Flags += $(Metis_path)/libmetis.dylib
	endif
endif

ifeq ($(OS_version), SI)
       Cholmod_lib_path=$(Cholmod_path:include/=lib/)
       Cholmod_Linking_Flags += $(Cholmod_lib_path)/libcholmod.a $(Cholmod_lib_path)/libamd.a $(Cholmod_lib_path)/libcolamd.a $(Cholmod_lib_path)/libsuitesparseconfig.a
       ifeq ($(UseMetis),yes)
#               Cholmod_Linking_Flags += $(Metis_path)/libmetis.dylib
			Cholmod_Linking_Flags += $(Cholmod_lib_path)/libccolamd.a $(Cholmod_lib_path)/libcamd.a $(Cholmod_lib_path)/libmetis.a
       endif
endif

ifeq ($(OS_version), OIST)
       Cholmod_lib_path=$(Cholmod_path:include/=lib/)
       Cholmod_Linking_Flags += $(Cholmod_lib_path)/libcholmod.a $(Cholmod_lib_path)/libamd.a $(Cholmod_lib_path)/libcolamd.a $(Cholmod_lib_path)/libsuitesparseconfig.a $(Cholmod_lib_path)/libcamd.a  $(Cholmod_lib_path)/libccolamd.a  $(Cholmod_lib_path)/libsuitesparseconfig.a
       ifeq ($(UseMetis),yes)
             Cholmod_Linking_Flags += $(Cholmod_lib_path)/libmetis.a
       endif
endif


OUT_CODE=LF_DEM
SRC=$(shell ls *.cpp)
LIB=$(shell ls *.h)
OBJ=$(SRC:.cpp=.o)
DSFMT_OBJ += $(DSFMT_SRC:.c=.o)
DATE=$(shell date +%D | sed "s./..g")

all:	CXXFLAGS= $(CXXFLAGS_optimized)
all:	$(OUT_CODE)

debug:	CXXFLAGS= $(CXXFLAGS_debug)
debug:	$(OUT_CODE)

dev:	CXXFLAGS= $(CXXFLAGS_dev)
dev:	$(OUT_CODE)

check:
	(cd check;	make test)

$(OUT_CODE): $(OBJ) $(DSFMT_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(DSFMT_OBJ)  -lstdc++ -lm $(Lapack_Linking_Flags) $(Blas_Linking_Flags) $(Cholmod_Linking_Flags)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(Includes) -o $@ -c $(@:.o=.cpp)

$(DSFMT_OBJ): $(DSFMT_SRC)
	$(CXX) $(CXXFLAGS) -o $@ -c $(@:.o=.c)

tar:
	echo "#define GIT_VERSION \"$(VERSION)\"" > VersionInfo.h
	tar -czvf LF_DEM_sources_$(DATE).tar.gz $(SRC) $(LIB) VersionInfo.h Makefile

clean:
	rm -f  $(OBJ)

doc:
	/Applications/Doxygen.app/Contents/Resources/doxygen Doxyfile

install:
	cp -p $(OUT_CODE) $(install_dir)
	(if [ "$(OS_version)" == "OIST" ]; then cp -p generateYaplotFile.pl $(install_dir); fi)
