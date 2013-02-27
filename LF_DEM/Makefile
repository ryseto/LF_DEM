#Makefile LF_DEM, Romain Mari

# Compiler and options
CXX=g++

# OS_version: OSX or Linux
OS_version=OSX
# Build : debug or optimized
Build=debug


#======== Include Paths =============
Trilinos_include_path_OSX    = /usr/local/include/
Cholmod_include_path_OSX     = /usr/local/include/
Trilinos_include_path_Linux  = /home/romain/romain/recherche/trilinos-11.0.3-Source/
Cholmod_include_path_Linux   = /usr/include/suitesparse/

Trilinos_path = $(Trilinos_include_path_$(OS_version))
Cholmod_path = $(Cholmod_include_path_$(OS_version))

# Includes
Includes = -I $(Cholmod_path) -I $(Trilinos_path)


#======== Linking ==================
Blas_Linking_Flags = -lblas 
Lapack_Linking_Flags = -llapack
Cholmod_Linking_Flags = -lcholmod
# extra linking flags for cholmod on Mac OSX
ifeq ($(OS_version),OSX)
	Cholmod_lib_path=$(Cholmod_path:include/=lib/)
	Cholmod_Linking_Flags += $(Cholmod_lib_path)/libcholmod.a $(Cholmod_lib_path)/libamd.a $(Cholmod_lib_path)/libcolamd.a $(Cholmod_lib_path)/libsuitesparseconfig.a
endif
Trilinos_Linking_Flags =  -lbelosepetra -lbelos -lifpack -lepetra -lkokkosnodetsqr -lteuchos


CXXFLAGS_Linux_optimized=  -O3 -Wall  -Wextra -Wshadow # linux
CXXFLAGS_Linux_debug=  -g -pg -Wall  -Wextra -Wshadow # linux
CXXFLAGS_OSX_optimized=  -O3 -Wall  -Wextra -Wshadow # mac OSX
CXXFLAGS_OSX_debug=  -g -pg

CXXFLAGS= $(CXXFLAGS_$(OS_version)_$(Build))


OUT_CODE=LF_DEM
SRC= $(shell ls *.cpp)
LIB=  $(shell ls *.h)  
OBJ=$(SRC:.cpp=.o)

DATE=$(shell date +%D | sed "s./..g")

all:	$(OUT_CODE)

$(OUT_CODE):$(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)  -lstdc++ -lm $(Lapack_Linking_Flags) $(Blas_Linking_Flags) $(Cholmod_Linking_Flags) $(Trilinos_Linking_Flags)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(Includes) -o $@ -c $(@:.o=.cpp) 

tar:
	tar cvfz LF_DEM_sources_$(DATE).tar.gz $(SRC) $(LIB)  LF_DEM_sources_$(DATE).ps Makefile input.*

clean:
	rm -f  $(OBJ)

ps:
	enscript --pretty-print --color --fancy-header -r -2 -p LF_DEM_sources_$(DATE).ps $(LIB) $(SRC)

