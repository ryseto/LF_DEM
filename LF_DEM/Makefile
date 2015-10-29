#Makefile LF_DEM, Romain Mari


# load Makefile inputs
makeconfig = config/Makefile_config.mk
gotconfig = 0

ifneq ("$(wildcard $(makeconfig))","")
	include $(makeconfig)
	gotconfig = 1
endif

# Compilation Flags
CXXFLAGS_optimized = -O3 -std=c++11
CXXFLAGS_debug = -g -DDEV="1" -std=c++11
CXXFLAGS_dev = -O3 -DDEV="1" -std=c++11

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
ifeq ($(CXX),gcc)
	CXXFLAGS_debug += -pg
endif

CXXFLAGS_optimized += $(CXXFLAGS_EXTRA)

ifeq ($(DSFMT_RNG),yes)
	CXXFLAGS_optimized += -DUSE_DSFMT="$DSFMT_RNG" -DDSFMT_MEXP=19937
	CXXFLAGS_debug += -DUSE_DSFMT= "$DSFMT_RNG" -DDSFMT_MEXP=19937
	DSFMT_SRC = dSFMT-src-2.2.3/dSFMT.c
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
Cholmod_path = $(SUITESPARSE_ROOT)/include/
Cholmod_libpath = $(SUITESPARSE_ROOT)/lib/
Cholmod_Linking_Flags = $(Cholmod_libpath)libcholmod.a $(Cholmod_libpath)libamd.a  $(Cholmod_libpath)libcolamd.a $(Cholmod_libpath)libsuitesparseconfig.a
#Metis_path = /opt/local/lib/

# Includes
Includes = -I $(Cholmod_path)

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

hasconfig:
	@(if [ $(gotconfig) -ne 1 ]; then \
		echo "\n no config/Makefile_config.mk \n \n";\
	fi)
	
$(OUT_CODE): $(OBJ) $(DSFMT_OBJ) hasconfig
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(DSFMT_OBJ)  $(Cholmod_Linking_Flags) $(Lapack_Linking_Flags) $(Blas_Linking_Flags) $(Extra_Linking_Flags) -lstdc++ -lm

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(Includes) -c $(@:.o=.cpp)

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
