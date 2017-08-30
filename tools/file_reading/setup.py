# Cython compile instructions

from distutils.core import setup, Extension
from Cython.Build import cythonize

clfdem_file_module = Extension('clfdem_file',
                               sources=['clfdem_file.pyx'],
                               language="c++",
                               extra_compile_args=["-std=c++11"])

setup(
  ext_modules=cythonize(clfdem_file_module),
)
