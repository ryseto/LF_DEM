import sys
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("Spherical_coordinate_histogram", ["Spherical_coordinate_histogram.pyx"], include_dirs=[numpy.get_include()]) ]

setup(
  name = 'Spherical coordinate container',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
