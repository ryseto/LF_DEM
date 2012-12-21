import sys
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("pair_correlation", ["pair_correlation.pyx"])]

setup(
  name = 'Computes spherically resolved pair correlation function for LF_DEM data',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
