#Makefile Romain Mari 06/2010
CXX=g++ # -cxxlib-icc-abr

INCLUDEPATHS= -I /usr/local/include/ -I /Users/LevichFellow/recherche/trilinos-11.0.3-Build/include/
#CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow -I /usr/local/include/ # mac OS
#CXXFLAGS=  -g -pg -Wall  -Wextra -Wshadow  $(INCLUDEPATHS) # mac OS debug
CXXFLAGS=  -g -pg $(INCLUDEPATHS) # mac OS debug no warning 
#CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow -I /usr/include/suitesparse/ # linux
#CXXFLAGS=  -g  -pg -Wall  -Wextra -Wshadow -I /usr/include/suitesparse/

LA_LIBRARIES = -llapack -lblas 

Cholmod_LIBRARIES = -lcholmod
Cholmod_LIBRARIES += /usr/local/Lib/libcholmod.a /usr/local/Lib/libamd.a /usr/local/Lib/libcolamd.a /usr/local/lib/libsuitesparseconfig.a # mac os

Trilinos_LIBRARY_DIRS= -L/Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib 
#Trilinos_LIBRARIES=  -lbelosepetra -lbelos -lml -lifpack -lamesos -lgaleri-xpetra -lgaleri -laztecoo -lisorropia -lepetraext -ltpetraext -ltpetrainout -ltpetra -ltriutils -lzoltan -lsimpi -lepetra -lkokkosdisttsqr -lkokkosnodetsqr -lkokkoslinalg -lkokkosnodeapi -lkokkos -lkokkosdisttsqr -lkokkosnodetsqr -lkokkoslinalg -lkokkosnodeapi -lkokkos -ltpi -lteuchos
Trilinos_LIBRARIES=  /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libbelosepetra.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libbelos.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libml.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libifpack.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libamesos.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libgaleri-xpetra.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libgaleri.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libaztecoo.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libisorropia.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libepetraext.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libtpetraext.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libtpetrainout.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libtpetra.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libtriutils.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libzoltan.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libsimpi.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libepetra.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkosdisttsqr.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkosnodetsqr.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkoslinalg.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkosnodeapi.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkos.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkosdisttsqr.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkosnodetsqr.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkoslinalg.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkosnodeapi.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libkokkos.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libtpi.dylib /Users/LevichFellow/recherche/trilinos-11.0.3-Build/lib/libteuchos.dylib


OUT_CODE=LF_DEM
SRC= $(shell ls *.cpp)
LIB=  $(shell ls *.h)  
OBJ=$(SRC:.cpp=.o)

DATE=$(shell date +%D | sed "s./..g")

all:	$(OUT_CODE)

$(OUT_CODE):$(OBJ) 
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ -lstdc++ -lm $(LA_LIBRARIES) $(Cholmod_LIBRARIES) $(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES)


$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ -c $(@:.o=.cpp) 

tar:
	tar cvfz LF_DEM_sources_$(DATE).tar.gz $(SRC) $(LIB)  LF_DEM_sources_$(DATE).ps Makefile input.*

clean:
	rm -f  $(OBJ)

ps:
	enscript --pretty-print --color --fancy-header -r -2 -p LF_DEM_sources_$(DATE).ps $(LIB) $(SRC)

