#=======================================#
#    Parameters to provide              #
#=======================================#

# Use DSFMT instead of MT as a RNG ( yes / no )
DSFMT_RNG = yes
# Enable use of Metis library ( yes / no )
UseMetis = yes

# the directory where LF_DEM will be copied on `make install`
install_dir = ~/bin/

# C++ compiler
CXX=clang
CC=clang

# Libraries
#
# SuiteSparse library install folder
SUITESPARSE_ROOT = /home/rjm238/usr/

# Extra flags to the compiler, if needed (e.g. optimization flags)
CXXFLAGS_EXTRA = -march=native # -flto


#======== Linking ==================
# Blas and Lapack (you may want to modify the BLAS, as -lblas might point to non optimized BLAS)
MKLROOT = /alt/applic/intel/intel13.1.0.146/mkl
Blas_Linking_Flags = -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 \
		   -Wl,--start-group \
               $(MKLROOT)/lib/intel64/libmkl_intel_lp64.so \
               $(MKLROOT)/lib/intel64/libmkl_core.so \
               $(MKLROOT)/lib/intel64/libmkl_gnu_thread.so \
               -Wl,--end-group -lpthread -lm -ldl
Lapack_Linking_Flags = 

# Extra linking here, if needed
# Extra_Linking_Flags =
# for OpenMP with g++
Extra_Linking_Flags = -fopenmp
