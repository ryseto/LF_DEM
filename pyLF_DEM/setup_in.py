#
# pyLF_DEM
#
###################################################################
#
# inputs to setup.py
#
# The purpose of this file is to define environment specific
# variables that will be used by distutils setup.py.
# These are essentially compilation flags for the different
# librairies used by LF_DEM
#
###################################################################

cc = "g++"
cxx = "g++"

LFDEM_ROOT = "../"
SUITESPARSE_ROOT = "/usr/local"
MKLROOT = "/opt/intel/mkl/"
Blas_Linking_Flags = ["-Wl,--no-as-needed", "-L${MKLROOT}/lib/intel64",
                      "-Wl,--start-group",
                      MKLROOT+"/lib/intel64/libmkl_intel_lp64.so",
                      MKLROOT+"/lib/intel64/libmkl_core.so",
                      MKLROOT+"/lib/intel64/libmkl_gnu_thread.so",
                      MKLROOT+"/lib/intel64/libmkl_avx.so",
                      MKLROOT+"/lib/intel64/libmkl_def.so",
                      "-Wl,--end-group",
                      "-lpthread",
                      "-lm",
                      "-ldl"]
Lapack_Linking_Flags = []

SuiteSparse_Linking_Flags = [SUITESPARSE_ROOT+"/lib/libcholmod.so",
                             '-Wl,-rpath='+SUITESPARSE_ROOT+"/lib/"]

Blas_Include_Dirs = []
Lapack_Include_Dirs = []
SuiteSparse_Include_Dirs = [SUITESPARSE_ROOT+"/include/"]

Extra_Compile_Args = []
Extra_Linking_Flags = ["-fopenmp"]

Library_Dirs = ["/usr/lib/"]
