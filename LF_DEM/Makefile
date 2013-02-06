#Makefile LF_DEM, Romain Mari
CXX=g++ # -cxxlib-icc-abr
trilinos_path=/usr/local
INCLUDEPATHS= -I /usr/local/include/ -I $(trilinos_path)/include/ 
#CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow -I /usr/local/include/ # mac OS
#CXXFLAGS=  -g -pg -Wall  -Wextra -Wshadow  $(INCLUDEPATHS) # mac OS debug
CXXFLAGS=  -g -pg # # mac OS debug no warning 
#CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow -I /usr/include/suitesparse/ # linux
#CXXFLAGS=  -g  -pg -Wall  -Wextra -Wshadow -I /usr/include/suitesparse/

LA_LIBRARIES = -llapack -lblas 

Cholmod_LIBRARIES = -lcholmod
Cholmod_LIBRARIES += /usr/local/lib/libcholmod.a /usr/local/lib/libamd.a /usr/local/lib/libcolamd.a /usr/local/lib/libsuitesparseconfig.a   # mac os

Trilinos_LIBRARY_DIRS= -L$(trilinos_path)/lib 
Trilinos_LIBRARIES=  -lbelosepetra -lbelos -lml -lifpack -lamesos -lgaleri-xpetra -lgaleri -laztecoo -lepetraext -ltpetraext -ltpetrainout -ltpetra -ltriutils -lepetra -lkokkosdisttsqr -lkokkosnodetsqr -lkokkoslinalg -lkokkosnodeapi -lkokkos -lkokkosdisttsqr -lkokkosnodetsqr -lkokkoslinalg -lkokkosnodeapi -lkokkos -ltpi -lteuchos


OUT_CODE=LF_DEM
SRC= $(shell ls *.cpp)
LIB=  $(shell ls *.h)  
OBJ=$(SRC:.cpp=.o)

DATE=$(shell date +%D | sed "s./..g")

all:	$(OUT_CODE)

$(OUT_CODE):$(OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)  -lstdc++ -lm $(LA_LIBRARIES) $(Cholmod_LIBRARIES) $(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES)


$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDEPATHS) -o $@ -c $(@:.o=.cpp) 

tar:
	tar cvfz LF_DEM_sources_$(DATE).tar.gz $(SRC) $(LIB)  LF_DEM_sources_$(DATE).ps Makefile input.*

clean:
	rm -f  $(OBJ)

ps:
	enscript --pretty-print --color --fancy-header -r -2 -p LF_DEM_sources_$(DATE).ps $(LIB) $(SRC)

