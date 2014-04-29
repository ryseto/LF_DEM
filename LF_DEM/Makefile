#Makefile LF_DEM, Romain Mari

#=======================================#
#    Parameters to provide              #
#=======================================#

# OS_version: OSX or Linux
OS_version=OSX
#SysEnv=Energy_seto
# Build : debug or optimized
#Build=optimized

# Compiler
ifeq ($(OS_version),Linux)
	CXX=g++
endif
ifeq ($(OS_version),OSX)
	CXX=clang
endif

# Enable use of Trilinos library ( yes / no )
UseTrilinos=no
# Enable use of Metis library ( yes / no )
UseMetis=no

#=======================================#
#=======================================#

#======== Compiler flags =============
CXXFLAGS_Linux_optimized=  -O3 -Wall  -Wextra -Wshadow
CXXFLAGS_Linux_debug=  -g -pg -Wall  -Wextra -Wshadow
CXXFLAGS_OSX_optimized=  -O3 -Wall  -Wextra -Wshadow
CXXFLAGS_OSX_debug=  -g -pg  -Wall  -Wextra -Wshadow

CXXFLAGS= $(CXXFLAGS_$(OS_version)_optimized)


#======== Include Paths =============
Cholmod_include_path_OSX     = /usr/local/include/

ifeq ($(SysEnv),Energy_seto)
	Cholmod_include_path_OSX     = ~/usr/local/include/
endif


Cholmod_include_path_Linux   = /usr/include/suitesparse/
ifeq ($(UseTrilinos),yes)
	Trilinos_include_path_Linux  = /home/romain/romain/recherche/trilinos-11.0.3-Source/
	Trilinos_include_path_OSX    = /usr/local/include/
endif


Trilinos_path = $(Trilinos_include_path_$(OS_version))
Cholmod_path = $(Cholmod_include_path_$(OS_version))
Metis_path = /opt/local/lib/

# Includes
Includes = -I $(Cholmod_path) 
ifeq ($(UseTrilinos),yes)
	Includes += -I $(Trilinos_path)
endif

#======== Linking ==================
# Blas and Lapack
Blas_Linking_Flags = -lblas 
Lapack_Linking_Flags = -llapack

# Cholmod
Cholmod_Linking_Flags = -lcholmod -lamd -lcamd -lcolamd -lccolamd
ifeq ($(SysEnv),Energy_seto)
	Cholmod_Linking_Flags     =
endif

ifeq ($(UseMetis),yes)
	Cholmod_Linking_Flags += -lmetis 
endif

# extra linking flags for cholmod on Mac OSX
ifeq ($(OS_version),OSX)
	Cholmod_lib_path=$(Cholmod_path:include/=lib/)
	Cholmod_Linking_Flags += $(Cholmod_lib_path)/libcholmod.a $(Cholmod_lib_path)/libamd.a $(Cholmod_lib_path)/libcolamd.a $(Cholmod_lib_path)/libsuitesparseconfig.a
	ifeq ($(UseMetis),yes)
		Cholmod_Linking_Flags += $(Metis_path)/libmetis.dylib
	endif
endif

# Trilinos
ifeq ($(UseTrilinos),yes)
	Trilinos_Linking_Flags =  -lbelosepetra -lbelos -lifpack -lepetra -lkokkosnodetsqr -lteuchos
endif



OUT_CODE=LF_DEM
SRC= $(shell ls *.cpp)
LIB=  $(shell ls *.h)  
OBJ=$(SRC:.cpp=.o)

DATE=$(shell date +%D | sed "s./..g")

all:	$(OUT_CODE)

debug:	CXXFLAGS= $(CXXFLAGS_$(OS_version)_debug)
debug:	$(OUT_CODE)

$(OUT_CODE):$(OBJ)
	$(CXX) $(CXXFLAGS_$(OS_version)_$(Build)) -o $@ $(OBJ)  -lstdc++ -lm $(Lapack_Linking_Flags) $(Blas_Linking_Flags) $(Cholmod_Linking_Flags) $(Trilinos_Linking_Flags)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(Includes) -o $@ -c $(@:.o=.cpp) 

tar:
	tar cvfz LF_DEM_sources_$(DATE).tar.gz $(SRC) $(LIB)  LF_DEM_sources_$(DATE).ps Makefile input.*

clean:
	rm -f  $(OBJ)

ps:
	enscript --pretty-print --color --fancy-header -r -2 -p LF_DEM_sources_$(DATE).ps $(LIB) $(SRC)

