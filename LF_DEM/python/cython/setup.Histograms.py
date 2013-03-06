import sys
import numpy
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("Histograms", ["Histograms.pyx"], include_dirs=[numpy.get_include()]) ]

setup(
  name = 'Abstract Histogram',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
