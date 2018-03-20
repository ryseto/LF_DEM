#!/usr/bin/env python

"""
setup.py file for pyLF_DEM
"""

from distutils.core import setup, Extension
import os
import subprocess
import setup_in

os.environ["CC"] = setup_in.cc
os.environ["CXX"] = setup_in.cxx

extra_link_args = setup_in.Blas_Linking_Flags + setup_in.Lapack_Linking_Flags\
                  + setup_in.SuiteSparse_Linking_Flags + setup_in.Extra_Linking_Flags
include_dirs = setup_in.Blas_Include_Dirs + setup_in.Lapack_Include_Dirs\
               + setup_in.SuiteSparse_Include_Dirs
extra_compile_args = setup_in.Extra_Compile_Args+["-std=c++11", "-O3"]

pylf_dem_version = subprocess.check_output(
                   ["git", "describe", "--dirty"])
pylf_dem_version = pylf_dem_version.decode("utf-8")[:-1]

extra_compile_args += ["-DGIT_VERSION=\""+pylf_dem_version+"\""]
library_dirs = setup_in.Library_Dirs

srcs = [x for x in os.listdir(setup_in.LFDEM_ROOT+"/LF_DEM/") if x.find(".cpp")+4==len(x)]

for i in range(len(srcs)):
    srcs[i] = setup_in.LFDEM_ROOT+"/LF_DEM/"+srcs[i]
srcs.append('pyLFDEM_wrap.cxx')


pyLFDEM_module = Extension('_pyLFDEM',
                           sources=srcs,
                           language="c++",
                           include_dirs=include_dirs,
                           extra_compile_args=extra_compile_args,
                           library_dirs=library_dirs,
                           extra_link_args=extra_link_args
                          )

setup(name='pyLFDEM',
      version='0.3',
      author="Romain Mari",
      description="""LF_DEM Simulation object""",
      ext_modules=[pyLFDEM_module],
      py_modules=["_pyLFDEM"],
      )
