#Makefile Romain Mari 06/2010
CXX=g++ # -cxxlib-icc-abr

#CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow -I /usr/local/include/ # mac OS
CXXFLAGS=  -g -pg -Wall  -Wextra -Wshadow -I /usr/local/include/ # mac OS debug
#CXXFLAGS=  -O3 -Wall  -Wextra -Wshadow -I /usr/include/suitesparse/ # linux
#CXXFLAGS=  -g  -pg -Wall  -Wextra -Wshadow -I /usr/include/suitesparse/
LINKLIBS = -llapack -lblas -lcholmod
LINKLIBS += /usr/local/Lib/libcholmod.a /usr/local/Lib/libamd.a /usr/local/Lib/libcolamd.a /usr/local/lib/libsuitesparseconfig.a # mac os


OUT_CODE=LF_DEM
SRC= $(shell ls *.cpp)
LIB=  $(shell ls *.h)  

OBJ=$(SRC:.cpp=.o)

DATE=$(shell date +%D | sed "s./..g")


all:	$(OUT_CODE)

$(OUT_CODE):$(OBJ) 
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ -lstdc++ -lm $(LINKLIBS)


$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ -c $(@:.o=.cpp)

tar: ps
	tar cvfz LF_DEM_sources_$(DATE).tar.gz $(SRC) $(LIB)  LF_DEM_sources_$(DATE).ps Makefile input.*

clean:
	rm -f  $(OBJ)

ps:
	enscript --pretty-print --color --fancy-header -r -2 -p LF_DEM_sources_$(DATE).ps $(LIB) $(SRC)

