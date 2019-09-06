#Makefile LF_DEM, Romain Mari

override_default_cholmod = false
UseMetis = no

# load Makefile inputs
makeconfig = config/Makefile_config.mk
gotconfig = 0
ifneq ("$(wildcard $(makeconfig))","")
	include $(makeconfig)
	gotconfig = 1
endif

# Compilation Flags
CXXFLAGS_optimized = -O3
CXXFLAGS_debug = -g -DDEV="1"
CXXFLAGS_dev = -g -O3 -DDEV="1"

GPP_Warnings = -Wall
ifeq ($(CXX),$(filter $(CXX),g++ clang++ gcc))
	CXXFLAGS_optimized += $(GPP_Warnings)
	CXXFLAGS_debug += -pg $(GPP_Warnings) -Wextra -Wshadow
	CXXFLAGS_dev += $(GPP_Warnings) -Wextra -Wshadow
endif
# ifeq ($(CXX),g++)
# 	CXXFLAGS_debug += -fprofile-arcs -ftest-coverage # profiling
# endif
ifeq ($(CXX),g++)
	CXXFLAGS_debug += -fno-inline
endif
ICC_Warnings = -Wall
ifeq ($(CXX),icpc)
	CXXFLAGS_optimized += $(ICC_Warnings)
	CXXFLAGS_debug += -pg $(ICC_Warnings)
	CXXFLAGS_dev += $(ICC_Warnings)
endif
ifeq ($(CXX),gcc)
	CXXFLAGS_debug += -pg
endif

CXXFLAGS_optimized += $(CXXFLAGS_EXTRA)

ifeq ($(DSFMT_RNG),yes)
	CXXFLAGS_optimized += -DUSE_DSFMT="$DSFMT_RNG" -DDSFMT_MEXP=19937
	CXXFLAGS_debug += -DUSE_DSFMT="$DSFMT_RNG" -DDSFMT_MEXP=19937
	DSFMT_SRC = dSFMT-src-2.2.3/dSFMT.c
endif

GSD_SRC = gsd.c

ifeq ($(UseMetis),yes)
	CXXFLAGS_optimized += -DUSE_METIS
	CXXFLAGS_debug += -DUSE_METIS
	CXXFLAGS_dev += -DUSE_METIS
endif

# Dealing with git versioning
is_git_repo := $(shell if [ -d .git ] || git rev-parse --git-dir > /dev/null 2>&1; then echo true;fi)
ifeq ($(is_git_repo),true)
	VERSION := $(shell git describe --dirty --always)
	CXXFLAGS_optimized += -DGIT_VERSION="\"$(VERSION)\""
	CXXFLAGS_debug += -DGIT_VERSION="\"$(VERSION)\""
	CXXFLAGS_dev += -DGIT_VERSION="\"$(VERSION)\""
endif


#======== Include Paths =============
ifeq ($(override_default_cholmod), false)
	Cholmod_path = -I $(SUITESPARSE_ROOT)/include/
	Cholmod_libpath = $(SUITESPARSE_ROOT)/lib/
ifneq ("$(wildcard $(Cholmod_libpath)libcholmod.so)","")
	Cholmod_Linking_Flags = $(Cholmod_libpath)libcholmod.so -Wl,-rpath=$(Cholmod_libpath)
else
	Cholmod_Linking_Flags = $(Cholmod_libpath)libcholmod.a $(Cholmod_libpath)libamd.a  $(Cholmod_libpath)libcolamd.a $(Cholmod_libpath)libsuitesparseconfig.a
endif
endif


# Includes
Includes = $(Cholmod_path)

OUT_CODE=LF_DEM
SRC=$(shell ls *.cpp)
LIB=$(shell ls *.h)
OBJ=$(SRC:.cpp=.o)
DSFMT_OBJ += $(DSFMT_SRC:.c=.o)
GSD_OBJ = $(GSD_SRC:.c=.o)
DATE=$(shell date +%D | sed "s./..g")

all:	CXXFLAGS= $(CXXFLAGS_optimized) -std=c++11
all:	CFLAGS= $(CXXFLAGS_optimized)
all:	$(OUT_CODE)

debug:	CXXFLAGS= $(CXXFLAGS_debug) -std=c++11
debug:	CFLAGS= $(CXXFLAGS_debug)
debug:	$(OUT_CODE)

dev:	CXXFLAGS= $(CXXFLAGS_dev) -std=c++11
dev:	CFLAGS= $(CXXFLAGS_dev)
dev:	$(OUT_CODE)

check:
	(cd check;	make test)

hasconfig:
	@(if [ $(gotconfig) -ne 1 ]; then \
		echo "\n no config/Makefile_config.mk \n \n";\
	fi)

$(OUT_CODE): $(OBJ) $(DSFMT_OBJ) $(GSD_OBJ) hasconfig
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(DSFMT_OBJ) $(GSD_OBJ) $(Cholmod_Linking_Flags) $(Cholmod_GPU) $(Lapack_Linking_Flags) $(Blas_Linking_Flags) $(Extra_Linking_Flags) -lstdc++ -lm

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(Includes) -c $(@:.o=.cpp)

$(GSD_OBJ): $(GSD_SRC)
	$(CC) $(CFLAGS) -o $@ -c $(@:.o=.c)

$(DSFMT_OBJ): $(DSFMT_SRC)
	$(CXX) $(CFLAGS) -o $@ -c $(@:.o=.c)

tar:
	echo "#define GIT_VERSION \"$(VERSION)\"" > VersionInfo.h
	tar -czvf LF_DEM_sources_$(DATE).tar.gz $(SRC) $(LIB) VersionInfo.h Makefile

clean:
	rm -f $(OBJ) $(DSFMT_OBJ) gsd.o

doc:
	doxygen Doxyfile

install:
	cp -p $(OUT_CODE) $(install_dir)
