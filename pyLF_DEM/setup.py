#!/usr/bin/env python

"""
setup.py file for pyLF_DEM
"""

from distutils.core import setup, Extension
import os
import subprocess

import re

LFDEM_ROOT = "../"

# find all variables defines in Makefile_config.mk
print("Parsing LF_DEM/config/Makefile_config.mk...")
mkconfig = {}
with open(LFDEM_ROOT+"LF_DEM/config/Makefile_config.mk") as f:
    for line in f:
        line = line.strip()
        if len(line) and line[0] != '#':
            print(line)
            k, val = line.split('=', maxsplit=1)
            k = k.strip()
            mkconfig[k] = val.lstrip().rstrip()

# replace all the ${...} by their value
for k in mkconfig:
    for m in re.finditer('\$\{\w*\}', mkconfig[k]):
         mkconfig[k] = mkconfig[k].replace(m.group(0), mkconfig[m.group(0)[2:-1]])
    for m in re.finditer('\$\(\w*\)', mkconfig[k]):
         mkconfig[k] = mkconfig[k].replace(m.group(0), mkconfig[m.group(0)[2:-1]])

os.environ["CC"] = mkconfig['CC']
os.environ["CXX"] = mkconfig['CXX']


extra_link_args = mkconfig['Blas_Linking_Flags'].split()
extra_link_args += [mkconfig['Lapack_Linking_Flags']]

SUITESPARSE_ROOT = mkconfig['SUITESPARSE_ROOT']
extra_link_args += [SUITESPARSE_ROOT+"/lib/libcholmod.so",
                             '-Wl,-rpath='+SUITESPARSE_ROOT+"/lib/"]
extra_link_args = [a for a in extra_link_args if len(a)]

include_dirs = [SUITESPARSE_ROOT+"/include/"]

extra_compile_args = ["-std=c++11", "-O3", "-fpermissive"] # -fpermissive for gsd.c

pylf_dem_version = subprocess.check_output(
                   ["git", "describe", "--dirty"])
pylf_dem_version = pylf_dem_version.decode("utf-8")[:-1]

extra_compile_args += ["-DGIT_VERSION=\""+pylf_dem_version+"\""]

srcs = [x for x in os.listdir(LFDEM_ROOT+"/LF_DEM/") if x.find(".cpp")+4==len(x)]
srcs.append(LFDEM_ROOT+"/LF_DEM/gsd.c")
for i in range(len(srcs)):
    srcs[i] = LFDEM_ROOT+"/LF_DEM/"+srcs[i]
srcs.append('pyLFDEM_wrap.cxx')


pyLFDEM_module = Extension('_pyLFDEM',
                           sources=srcs,
                           language="c++",
                           include_dirs=include_dirs,
                           extra_compile_args=extra_compile_args,
                           extra_link_args=extra_link_args
                          )

setup(name='pyLFDEM',
      version='0.3',
      author="Romain Mari",
      description="""LF_DEM Simulation object""",
      ext_modules=[pyLFDEM_module],
      py_modules=["_pyLFDEM"],
      )
