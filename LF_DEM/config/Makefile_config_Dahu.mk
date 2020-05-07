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
CXX = g++
CC = gcc

CUDA_HOME = /usr/

NIX_ROOT = /nix/var/nix/profiles/per-user/rmari/lfdem/

# Libraries
#
# SuiteSparse library install folder
#SUITESPARSE_ROOT = /home/rmari/lib/
SUITESPARSE_ROOT = ~/lib/

# Extra flags to the compiler, if needed (e.g. optimization flags)
CXXFLAGS_EXTRA = -march=native #-DUSE_GPU=1
Cholmod_GPU_libpath = $(SUITESPARSE_ROOT)/lib/
Cholmod_GPU = #$(Cholmod_GPU_libpath)libSuiteSparse_GPURuntime.so $(Cholmod_GPU_libpath)libGPUQREngine.so

#======== Linking ==================

# Blas and Lapack (you may want to modify the BLAS, as -lblas might point to non optimized BLAS)
#Blas_Linking_Flags = -lblas
#Lapack_Linking_Flags = -llapack
# Intel MKL
Blas_Linking_Flags = -Wl,--no-as-needed  ${NIX_ROOT}/lib/libmkl_intel_lp64.so  ${NIX_ROOT}/lib/libmkl_intel_thread.so  ${NIX_ROOT}/lib/libmkl_core.so ${NIX_ROOT}/lib/libiomp5.so -lpthread -lm -ldl
Lapack_Linking_Flags =

# Extra linking here, if needed
Extra_Linking_Flags = #$(CUDA_HOME)/lib/x86_64-linux-gnu/libcudart.so
# for OpenMP with g++
# Extra_Linking_Flags = -fopenmp
