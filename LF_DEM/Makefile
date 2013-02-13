#Makefile LF_DEM, Romain Mari

# Compiler and options
CXX=g++

OS_version=OSX # OSX or Linux
Opt=Debug # debug or optimized


Trilinos_path=/usr/local
Cholmod_path=/usr/local

# Includes
Includes = -I $(Cholmod_path)/include/ -I $(Trilinos_path)/include/ 

# Linking
Blas_Linking_Flags = -lblas 
Lapack_Linking_Flags = -llapack
Cholmod_Linking_Flags = -lcholmod
# extra linking flags for mac OS X
Cholmod_Linking_Flags += $(Cholmod_path)/lib/libcholmod.a $(Cholmod_path)/lib/libamd.a $(Cholmod_path)/lib/libcolamd.a $(Cholmod_path)/lib/libsuitesparseconfig.a
Trilinos_Linking_Flags =  -lbelosepetra -lbelos -lifpack -lepetra -lkokkosnodetsqr -lteuchos


CXXFLAGS= -O3

# tests do not work !!!???!
ifeq ($(OS_version),OSX)
ifeq ($(Opt), optimized)
		CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow # mac OS
else
		CXXFLAGS=  -g -pg # mac OS debug
endif
endif
ifeq ($(OS_version), Linux)
ifeq ($(Opt), optimized)
		CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow  # linux
else
		CXXFLAGS=  -g  -pg -Wall  -Wextra -Wshadow # linux debug
endif
endif



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

