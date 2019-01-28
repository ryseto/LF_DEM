#=======================================#
#    Parameters to provide              #
#=======================================#


# Use DSFMT instead of MT as a RNG ( yes / no )
DSFMT_RNG = no
# Enable use of Metis library ( yes / no )
UseMetis = yes

# the directory where LF_DEM will be copied on `make install`
install_dir = ~/bin/

# C++ compiler
CXX = icpc
CC = icc

# Libraries
#
# SuiteSparse library install folder
SUITESPARSE_ROOT = /applis/site/stow/gcc_4.4.6/suitesparse_4.2.1/

# Extra flags to the compiler, if needed (e.g. optimization flags)
CXXFLAGS_EXTRA = -xAVX


#======== Linking ==================

# Blas and Lapack (you may want to modify the BLAS, as -lblas might point to non optimized BLAS)
Blas_Linking_Flags = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl 
#Lapack_Linking_Flags = -llapack
# Intel MKL
# Blas_Linking_Flags = -mkl -lrt
# Lapack_Linking_Flags =

# Extra linking here, if needed
Extra_Linking_Flags = -lrt
# for OpenMP with g++
# Extra_Linking_Flags = -fopenmp
