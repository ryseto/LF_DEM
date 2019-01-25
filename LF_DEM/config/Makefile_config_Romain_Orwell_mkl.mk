#=======================================#
#    Parameters to provide              #
#=======================================#


# Use DSFMT instead of MT as a RNG ( yes / no )
DSFMT_RNG = yes
# Enable use of Metis library ( yes / no )
UseMetis = no

# the directory where LF_DEM will be copied on `make install`
install_dir = ~/usr/bin/

# C++ compiler
CXX = clang++
CC = clang

# Libraries
#
# SuiteSparse library install folder
SUITESPARSE_ROOT = /usr/local

# Extra flags to the compiler, if needed (e.g. optimization flags)
CXXFLAGS_EXTRA =


#======== Linking ==================

# Blas and Lapack (you may want to modify the BLAS, as -lblas might point to non optimized BLAS)
MKLROOT = /opt/intel/mkl/
Blas_Linking_Flags = -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 \
		   -Wl,--start-group \
               $(MKLROOT)/lib/intel64/libmkl_intel_lp64.so \
               $(MKLROOT)/lib/intel64/libmkl_core.so \
               $(MKLROOT)/lib/intel64/libmkl_gnu_thread.so \
               -Wl,--end-group -fopenmp=libomp -lpthread -lm -ldl
Lapack_Linking_Flags =

# Extra linking here, if needed
Extra_Linking_Flags =
# for OpenMP with g++
# Extra_Linking_Flags = -fopenmp
