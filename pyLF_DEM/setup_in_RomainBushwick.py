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

import re

cc = "g++"
cxx = "g++"

LFDEM_ROOT = "../"

# find all variables defines in Makefile_config.mk
mkconfig = {}
with open(LFDEM_ROOT+"LF_DEM/config/Makefile_config.mk") as f:
    for line in f:
        line = line.strip()
        if len(line) and line[0] != '#':
            k, val = line.split(' =', maxsplit=1)
            mkconfig[k] = val.lstrip().rstrip()

# replace all the ${...} by their value
for k in mkconfig:
    for m in re.finditer('\$\{\w*\}', mkconfig[k]):
         mkconfig[k] = mkconfig[k].replace(m.group(0), mkconfig[m.group(0)[2:-1]])
         

SUITESPARSE_ROOT = mkconfig['SUITESPARSE_ROOT']
Blas_Linking_Flags = mkconfig['Blas_Linking_Flags'].split()
Lapack_Linking_Flags = [mkconfig['Lapack_Linking_Flags']]
SuiteSparse_Linking_Flags = [SUITESPARSE_ROOT+"/lib/libcholmod.so",
                             '-Wl,-rpath='+SUITESPARSE_ROOT+"/lib/"]
Blas_Include_Dirs = []
Lapack_Include_Dirs = []
SuiteSparse_Include_Dirs = [SUITESPARSE_ROOT+"/include/"]

Extra_Compile_Args = []
Extra_Linking_Flags = []

Library_Dirs = ["/usr/lib/"]
