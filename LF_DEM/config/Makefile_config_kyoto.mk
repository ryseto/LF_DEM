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
CXX = icpc
CC = icc

# Libraries
#
# SuiteSparse library install folder
SUITESPARSE_ROOT = /home/seto/SuiteSparse/

# Extra flags to the compiler, if needed (e.g. optimization flags)
CXXFLAGS_EXTRA = 
CFLAGS_EXTRA = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 \
       -Wall -DWITHOUT_MPI -D__ArrayExtensions -DNDEBUG \
       -fomit-frame-pointer -DHAVE_SSE2=1 -DDSFMT_MEXP=19937 \
       -I$(HDF_DIR)/include -I/usr/local/include

#======== Linking ==================

# Blas and Lapack (you may want to modify the BLAS, as -lblas might point to non optimized BLAS)
Blas_Linking_Flags = -mkl -lrt
Lapack_Linking_Flags = 
# Intel MKL
# Blas_Linking_Flags = -mkl -lrt
# Lapack_Linking_Flags =

# Extra linking here, if needed
Extra_Linking_Flags =
# for OpenMP with g++
# Extra_Linking_Flags = -fopenmp
