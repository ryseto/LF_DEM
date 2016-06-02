#=======================================#
#    Parameters to provide              #
#=======================================#

# Use DSFMT instead of MT as a RNG ( yes / no )
DSFMT_RNG = no
# Enable use of Metis library ( yes / no )
UseMetis = no

# the directory where LF_DEM will be copied on `make install`
install_dir = ~/bin/

# C++ compiler
CXX=g++

# Libraries
#
# SuiteSparse library install folder
SUITESPARSE_ROOT = /home/rjm238/usr/

# Extra flags to the compiler, if needed (e.g. optimization flags)
CXXFLAGS_EXTRA = -march=native
# -flto


#======== Linking ==================
# Blas and Lapack (you may want to modify the BLAS, as -lblas might point to non optimized BLAS)
Blas_Linking_Flags = /usr/lib/libopenblas.so.0
Lapack_Linking_Flags = -llapack

# Extra linking here, if needed
# Extra_Linking_Flags =
# for OpenMP with g++
Extra_Linking_Flags = -fopenmp
