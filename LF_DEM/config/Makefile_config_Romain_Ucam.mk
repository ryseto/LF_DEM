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
CXX=g++
CC=gcc

# Libraries
#
# SuiteSparse library install folder
SUITESPARSE_ROOT = /alt/applic/user-maint/rjm238/suitesparse/

# Extra flags to the compiler, if needed (e.g. optimization flags)
CXXFLAGS_EXTRA =


#======== Linking ==================
# Blas and Lapack (you may want to modify the BLAS, as -lblas might point to non optimized BLAS)
Blas_Linking_Flags = /usr/lib/libopenblas.so.0
Lapack_Linking_Flags = -llapack

# Extra linking here, if needed
# Extra_Linking_Flags =
# for OpenMP with g++
Extra_Linking_Flags = -fopenmp
