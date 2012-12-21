import sys
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("cart2sph", ["cart2sph.pyx"])]

setup(
  name = 'cartesian to spherical coordinates',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
